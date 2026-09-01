#!/usr/bin/env Rscript
# ============================================================================
# Oklahoma Induced Seismicity Analysis — Bivariate ETAS with Sensitivity
#
# Treatment: OCC directive AOI_20150318 (wastewater injection reduction)
# State space: Oklahoma state boundary
#
# Core fits (primary tessellation, default grid_1.0R):
#   A. Naive bivariate ETAS     — joint MLE, homogeneous background
#   B. SEM bivariate ETAS       — adaptive_SEM, homogeneous background
#   C. Naive bivariate + KDE bg — publication naive pair
#   D. SEM bivariate + KDE bg   — publication SEM pair
#
# Partition sensitivity (excludes the primary tessellation):
#   - County boundaries
#   - Remaining grids: 2D, 5D (and 1D if county is primary)
#   - AOI region (AOI polygon = treated, rest = control)
#
# Usage:
#   Rscript oklahoma_analysis.R
#   Rscript oklahoma_analysis.R --test
# ============================================================================

suppressPackageStartupMessages({
  library(spatstat)
  library(sf)
  library(tigris)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
TEST_MODE <- "--test" %in% args
QUICK_CHECK <- "--quick-check" %in% args

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[1]),
                        winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
REPO_DIR <- if (basename(SCRIPT_DIR) == "oklahoma" &&
                basename(dirname(SCRIPT_DIR)) == "inst") {
  normalizePath(dirname(dirname(SCRIPT_DIR)), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
# R --file= can mangle spaces (Google Drive "My Drive" → "My~+~Drive").
repo_dir_env <- trimws(Sys.getenv("OK_REPO_DIR", ""))
if (nzchar(repo_dir_env) && file.exists(file.path(repo_dir_env, "DESCRIPTION"))) {
  REPO_DIR <- normalizePath(repo_dir_env, winslash = "/", mustWork = FALSE)
  SCRIPT_DIR <- file.path(REPO_DIR, "inst", "oklahoma")
} else if (!file.exists(file.path(REPO_DIR, "DESCRIPTION"))) {
  stop("Could not locate PPDisentangle DESCRIPTION. Set OK_REPO_DIR to the package root.")
}

# Ensure this script uses the CURRENT repository source (not a stale installed build).
# This is critical for SEM relabelling logic diagnostics and recent algorithm changes.
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required. Install with install.packages('pkgload').")
}
tryCatch({
  pkgload::load_all(REPO_DIR, quiet = TRUE, export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)
  cat("Loaded PPDisentangle from local source via pkgload::load_all().\n")
}, error = function(e) {
  stop("Failed to load local PPDisentangle source via pkgload::load_all(): ", e$message)
})
# Canonical output path (outside the git repo):
#   ../PPDisentangle-output/oklahoma/
OUT_DIR  <- PPDisentangle::pp_output_path("oklahoma", repo_root = REPO_DIR)
PLOT_DIR <- file.path(OUT_DIR, "plots")
report_out_dir_raw <- trimws(Sys.getenv("OK_REPORT_OUTPUT_DIR", ""))
REPORT_OUT_DIR <- if (nzchar(report_out_dir_raw)) {
  normalizePath(report_out_dir_raw, winslash = "/", mustWork = FALSE)
} else {
  OUT_DIR
}
for (d in unique(c(OUT_DIR, PLOT_DIR, REPORT_OUT_DIR))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
remove_legacy_report_outputs <- function() {
  legacy_paths <- c(
    file.path(SCRIPT_DIR, "oklahoma_report.html"),
    file.path(SCRIPT_DIR, "oklahoma_report.pdf"),
    file.path(SCRIPT_DIR, "oklahoma_report.tex")
  )
  legacy_glob <- list.files(
    SCRIPT_DIR,
    pattern = "^oklahoma_report_[0-9]{8}_[0-9]{6}\\.(html|pdf)$",
    full.names = TRUE
  )
  for (p in c(legacy_paths, legacy_glob)) {
    if (file.exists(p)) unlink(p)
  }
  legacy_files_dir <- file.path(SCRIPT_DIR, "oklahoma_report_files")
  if (dir.exists(legacy_files_dir)) unlink(legacy_files_dir, recursive = TRUE)
}
path_rel_to <- function(path, start) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  start <- normalizePath(start, winslash = "/", mustWork = FALSE)
  path_seg <- strsplit(path, "/", fixed = TRUE)[[1L]]
  start_seg <- strsplit(start, "/", fixed = TRUE)[[1L]]
  path_seg <- path_seg[nzchar(path_seg)]
  start_seg <- start_seg[nzchar(start_seg)]
  n <- min(length(path_seg), length(start_seg))
  i <- 1L
  while (i <= n && identical(path_seg[[i]], start_seg[[i]])) i <- i + 1L
  rel_parts <- c(rep("..", length(start_seg) - i + 1L), path_seg[i:length(path_seg)])
  if (length(rel_parts) < 1L) "." else paste(rel_parts, collapse = "/")
}
analysis_plots <- list()
kde_info <- NULL
store_analysis_plot <- function(name, plot_obj) {
  if (is.null(plot_obj)) return(invisible(NULL))
  analysis_plots[[name]] <<- plot_obj
  invisible(NULL)
}
SLURM_JOB_ID_RAW <- trimws(Sys.getenv("SLURM_JOB_ID", ""))
FILE_TAG <- if (nzchar(SLURM_JOB_ID_RAW)) paste0("_job", SLURM_JOB_ID_RAW) else ""
add_file_tag <- function(filename) {
  if (!nzchar(FILE_TAG)) return(filename)
  dot <- regexpr("\\.[^.]+$", filename, perl = TRUE)
  if (dot[1] > 0L) {
    paste0(substr(filename, 1L, dot[1] - 1L), FILE_TAG, substr(filename, dot[1], nchar(filename)))
  } else {
    paste0(filename, FILE_TAG)
  }
}

# ---- Configuration ----
DATA_DIR   <- file.path(SCRIPT_DIR, "oklahoma_induced_seismicity_data_regional20150318")
# Magnitude cutoff (catalog + ETAS reference m0). Prefer matching metadata$min_mag.
ETAS_M0    <- suppressWarnings(as.numeric(Sys.getenv("OK_ETAS_M0", "2.5")))
if (!is.finite(ETAS_M0) || ETAS_M0 <= 0) ETAS_M0 <- 2.5
# Gutenberg-Richter beta: estimated from pre-treatment *control* magnitudes after
# data load unless OK_BETA_GR is set explicitly.
BETA_GR_ENV <- suppressWarnings(as.numeric(Sys.getenv("OK_BETA_GR", "")))
BETA_GR    <- if (is.finite(BETA_GR_ENV) && BETA_GR_ENV > 0) BETA_GR_ENV else NA_real_
BETA_GR_SOURCE <- if (is.finite(BETA_GR)) "env:OK_BETA_GR" else "pending_pre_control_estimate"
CRS_PROJ   <- 5070

VANILLA_MAXIT <- if (QUICK_CHECK) 120 else if (TEST_MODE) 500 else 5000
env_vanilla_maxit <- suppressWarnings(as.integer(Sys.getenv("OK_VANILLA_MAXIT", "")))
if (!is.na(env_vanilla_maxit) && env_vanilla_maxit > 0L) {
  VANILLA_MAXIT <- env_vanilla_maxit
}
# Spatial scale is magnitude-independent: d(m) = D (gamma fixed at 0 in all fits).
GAMMA_FIXED <- 0
# Interior branching-ratio starts (~0.5-0.85 after GR scaling) to avoid
# Nelder-Mead beginning on the hard eta=0.98 cliff.
VANILLA_STARTS <- list(
  list(mu = 1.0, A = 0.20, alpha_m = 0.8, c = 0.05, p = 2.2,
       D = 5.0, gamma = GAMMA_FIXED, q = 1.6),
  list(mu = 2.0, A = 0.10, alpha_m = 0.5, c = 0.1,  p = 2.2,
       D = 2.0, gamma = GAMMA_FIXED, q = 2.0),
  list(mu = 1.5, A = 0.35, alpha_m = 0.8, c = 0.05, p = 2.2,
       D = 5.0, gamma = GAMMA_FIXED, q = 1.6),
  list(mu = 1.0, A = 0.55, alpha_m = 0.5, c = 0.08, p = 2.2,
       D = 3.0, gamma = GAMMA_FIXED, q = 1.8)
)
# Merge gamma=0 into any fixed_params list (all models exclude gamma).
with_gamma_fixed <- function(fixed_params = NULL) {
  out <- if (is.null(fixed_params)) list() else as.list(fixed_params)
  out[["gamma"]] <- GAMMA_FIXED
  out
}

SEM_N_LABELLINGS  <- if (QUICK_CHECK) 3 else if (TEST_MODE) 5  else 20
env_sem_n_labellings <- suppressWarnings(as.integer(Sys.getenv("OK_SEM_N_LABELLINGS", "")))
if (!is.na(env_sem_n_labellings) && env_sem_n_labellings >= 0L) {
  SEM_N_LABELLINGS <- env_sem_n_labellings
}
SEM_N_ITER        <- if (QUICK_CHECK) 1 else if (TEST_MODE) 2 else 10
env_sem_n_iter <- suppressWarnings(as.integer(Sys.getenv("OK_SEM_N_ITER", "")))
if (!is.na(env_sem_n_iter) && env_sem_n_iter > 0L) {
  SEM_N_ITER <- env_sem_n_iter
}
SEM_INNER_ITER    <- if (QUICK_CHECK) 2 else if (TEST_MODE) 5 else 1000
sem_inner_raw <- Sys.getenv("OK_SEM_INNER_ITER", "")
if (!nzchar(sem_inner_raw)) sem_inner_raw <- Sys.getenv("PP_SEM_INNER", "")
env_sem_inner_iter <- suppressWarnings(as.integer(sem_inner_raw))
if (!is.na(env_sem_inner_iter) && env_sem_inner_iter > 0L) {
  SEM_INNER_ITER <- env_sem_inner_iter
}
SENS_SEM_INNER_ITER <- suppressWarnings(as.integer(Sys.getenv("OK_SENS_SEM_INNER_ITER", as.character(SEM_INNER_ITER))))
if (!is.finite(SENS_SEM_INNER_ITER) || is.na(SENS_SEM_INNER_ITER) || SENS_SEM_INNER_ITER < 1L) {
  SENS_SEM_INNER_ITER <- SEM_INNER_ITER
}
SEM_INNER_PROPS   <- if (QUICK_CHECK) 3 else if (TEST_MODE) 5  else 20
env_sem_inner_props <- suppressWarnings(as.integer(Sys.getenv("OK_SEM_INNER_PROPS", "")))
if (is.finite(env_sem_inner_props) && !is.na(env_sem_inner_props) && env_sem_inner_props >= 1L) {
  SEM_INNER_PROPS <- env_sem_inner_props
}
SEM_CHANGE_FACTOR <- 0.01
SEM_CHANGE_FACTOR_MIN_MULT <- suppressWarnings(as.numeric(Sys.getenv("OK_SEM_CHANGE_FACTOR_MIN_MULT", "1.0")))
if (!is.finite(SEM_CHANGE_FACTOR_MIN_MULT) || is.na(SEM_CHANGE_FACTOR_MIN_MULT) || SEM_CHANGE_FACTOR_MIN_MULT <= 0) {
  SEM_CHANGE_FACTOR_MIN_MULT <- 1.0
}
SEM_CHANGE_FACTOR_MAX_MULT <- suppressWarnings(as.numeric(Sys.getenv("OK_SEM_CHANGE_FACTOR_MAX_MULT", "5.0")))
if (!is.finite(SEM_CHANGE_FACTOR_MAX_MULT) || is.na(SEM_CHANGE_FACTOR_MAX_MULT) || SEM_CHANGE_FACTOR_MAX_MULT < SEM_CHANGE_FACTOR_MIN_MULT) {
  SEM_CHANGE_FACTOR_MAX_MULT <- max(5.0, SEM_CHANGE_FACTOR_MIN_MULT)
}
SEM_OPTIM_METHOD <- tolower(trimws(Sys.getenv("OK_SEM_OPTIM_METHOD", "sample_weighted")))
if (!SEM_OPTIM_METHOD %in% c("max", "sample_weighted")) SEM_OPTIM_METHOD <- "sample_weighted"
SEM_MONOTONE_COMPLETE_LL <- tolower(Sys.getenv("OK_SEM_MONOTONE_COMPLETE_LL", "false")) %in%
  c("1", "true", "yes", "y")
SEM_START_FROM_C <- tolower(Sys.getenv("OK_SEM_START_FROM_C", "false")) %in%
  c("1", "true", "yes", "y")
SEM_BIV_N_THREADS <- suppressWarnings(as.integer(Sys.getenv("OK_SEM_BIV_N_THREADS", "1")))
if (!is.finite(SEM_BIV_N_THREADS) || is.na(SEM_BIV_N_THREADS) || SEM_BIV_N_THREADS < 1L) {
  SEM_BIV_N_THREADS <- 1L
}
SEM_SINGLE_FLIP_FROM_ITER <- suppressWarnings(as.numeric(Sys.getenv(
  "OK_SEM_SINGLE_FLIP_FROM_ITER", ""
)))
if (!is.finite(SEM_SINGLE_FLIP_FROM_ITER) || is.na(SEM_SINGLE_FLIP_FROM_ITER) ||
    SEM_SINGLE_FLIP_FROM_ITER < 1) {
  SEM_SINGLE_FLIP_FROM_ITER <- Inf
} else {
  SEM_SINGLE_FLIP_FROM_ITER <- as.integer(SEM_SINGLE_FLIP_FROM_ITER)
}
SEM_SELECTION_TEMPERATURE <- suppressWarnings(as.numeric(Sys.getenv("OK_SEM_SELECTION_TEMPERATURE", "0.2")))
if (!is.finite(SEM_SELECTION_TEMPERATURE) || is.na(SEM_SELECTION_TEMPERATURE) || SEM_SELECTION_TEMPERATURE <= 0) {
  SEM_SELECTION_TEMPERATURE <- 0.2
}
SEM_MAX_RELABEL_STEP_FRAC <- suppressWarnings(as.numeric(Sys.getenv("OK_SEM_MAX_RELABEL_STEP_FRAC", "0.05")))
if (!is.finite(SEM_MAX_RELABEL_STEP_FRAC) || is.na(SEM_MAX_RELABEL_STEP_FRAC) || SEM_MAX_RELABEL_STEP_FRAC <= 0) {
  SEM_MAX_RELABEL_STEP_FRAC <- 0.05
}
SEM_MAX_RELABEL_STEP_FRAC <- min(SEM_MAX_RELABEL_STEP_FRAC, 1.0)
SEM_FORCE_PARAM_UPDATE_FLIP_FRAC <- suppressWarnings(as.numeric(Sys.getenv("OK_SEM_FORCE_PARAM_UPDATE_FLIP_FRAC", "0.05")))
if (!is.finite(SEM_FORCE_PARAM_UPDATE_FLIP_FRAC) || is.na(SEM_FORCE_PARAM_UPDATE_FLIP_FRAC) || SEM_FORCE_PARAM_UPDATE_FLIP_FRAC <= 0) {
  SEM_FORCE_PARAM_UPDATE_FLIP_FRAC <- 0.05
}
SEM_FORCE_PARAM_UPDATE_FLIP_FRAC <- min(SEM_FORCE_PARAM_UPDATE_FLIP_FRAC, 1.0)
RUN_SEM_PILOT <- tolower(trimws(Sys.getenv("OK_RUN_SEM_PILOT", "false"))) %in% c("1", "true", "yes", "y")
SEM_PILOT_INNER_ITER <- suppressWarnings(as.integer(Sys.getenv("OK_SEM_PILOT_INNER_ITER", "100")))
if (!is.finite(SEM_PILOT_INNER_ITER) || is.na(SEM_PILOT_INNER_ITER) || SEM_PILOT_INNER_ITER < 1L) {
  SEM_PILOT_INNER_ITER <- 100L
}
# Pilot is intentionally configured to saturate available cores (one fit/core).
SEM_PILOT_CORES <- suppressWarnings(as.integer(Sys.getenv("OK_SEM_PILOT_CORES", "")))
if (!is.finite(SEM_PILOT_CORES) || is.na(SEM_PILOT_CORES) || SEM_PILOT_CORES < 1L) {
  SEM_PILOT_CORES <- NA_integer_
}
SEM_PILOT_MAX_COMBOS <- suppressWarnings(as.integer(Sys.getenv("OK_SEM_PILOT_MAX_COMBOS", "24")))
if (!is.finite(SEM_PILOT_MAX_COMBOS) || is.na(SEM_PILOT_MAX_COMBOS) || SEM_PILOT_MAX_COMBOS < 1L) {
  SEM_PILOT_MAX_COMBOS <- 24L
}
SEM_STAGNATION_TRIGGER_EVERY <- 50
SEM_TEMPORAL_WEIGHT <- suppressWarnings(as.numeric(Sys.getenv("OK_SEM_TEMPORAL_WEIGHT", "0")))
if (!is.finite(SEM_TEMPORAL_WEIGHT) || is.na(SEM_TEMPORAL_WEIGHT)) SEM_TEMPORAL_WEIGHT <- 0
SEM_TEMPORAL_WEIGHT <- min(max(SEM_TEMPORAL_WEIGHT, 0), 1)
SEM_TEMPORAL_SCALE_DAYS <- 15
SEM_T_TRUNC_RAW <- trimws(Sys.getenv("OK_SEM_T_TRUNC_DAYS", "90"))
if (!nzchar(SEM_T_TRUNC_RAW)) {
  stop("OK_SEM_T_TRUNC_DAYS is empty; auto t_trunc is disabled. Set a positive number of days (typically 90).")
}
SEM_T_TRUNC_DAYS <- suppressWarnings(as.numeric(SEM_T_TRUNC_RAW))
if (!is.finite(SEM_T_TRUNC_DAYS) || is.na(SEM_T_TRUNC_DAYS) || SEM_T_TRUNC_DAYS <= 0) {
  stop(sprintf(
    "OK_SEM_T_TRUNC_DAYS must be a positive number of days (got %s). Auto t_trunc is disabled.",
    SEM_T_TRUNC_RAW
  ))
}
SEM_T_TRUNC_SOURCE <- if (nzchar(trimws(Sys.getenv("OK_SEM_T_TRUNC_DAYS", "")))) "env" else "default_90"
# Kept only for saved-config / report compatibility; auto t_trunc is disabled.
SEM_T_TRUNC_REL <- suppressWarnings(as.numeric(Sys.getenv("OK_SEM_T_TRUNC_REL", "0.05")))
if (!is.finite(SEM_T_TRUNC_REL) || is.na(SEM_T_TRUNC_REL) || SEM_T_TRUNC_REL <= 0 || SEM_T_TRUNC_REL >= 1) {
  SEM_T_TRUNC_REL <- 0.05
}
SEM_PARAM_UPDATE  <- if (QUICK_CHECK) 10 else if (TEST_MODE) 10 else 25
SEM_OUTER_MAXIT       <- if (QUICK_CHECK) 40 else if (TEST_MODE) 200 else 5000
SEM_OUTER_MAXIT_BIV   <- 5000
env_sem_outer_maxit <- suppressWarnings(as.integer(Sys.getenv("OK_SEM_OUTER_MAXIT", "")))
if (!is.na(env_sem_outer_maxit) && env_sem_outer_maxit > 0L) {
  SEM_OUTER_MAXIT <- env_sem_outer_maxit
}
env_sem_outer_maxit_biv <- suppressWarnings(as.integer(Sys.getenv("OK_SEM_OUTER_MAXIT_BIV", "")))
if (!is.na(env_sem_outer_maxit_biv) && env_sem_outer_maxit_biv > 0L) {
  SEM_OUTER_MAXIT_BIV <- env_sem_outer_maxit_biv
}
SEM_WARMSTART_FIXED <- tolower(Sys.getenv("OK_SEM_WARMSTART_FIXED", "false")) %in% c("1", "true", "yes", "y")
# Optional speed mode for proof-of-concept runs.
RUN_SENSITIVITY <- tolower(Sys.getenv("OK_RUN_SENSITIVITY", "true")) %in% c("1", "true", "yes", "y")
# Partition grids/AOI are independent of KDE bandwidth sensitivity.
RUN_PARTITION_SENSITIVITY <- tolower(Sys.getenv(
  "OK_RUN_PARTITION_SENSITIVITY",
  if (RUN_SENSITIVITY) "true" else "false"
)) %in% c("1", "true", "yes", "y")
# Absolute-km KDE bandwidth grid for Fit C/D (empty => Diggle multipliers when RUN_SENSITIVITY).
KDE_BW_SENS_KM_RAW <- trimws(Sys.getenv("OK_KDE_BW_SENS_KM", ""))
KDE_BW_SENS_KM <- if (nzchar(KDE_BW_SENS_KM_RAW)) {
  suppressWarnings(as.numeric(unlist(strsplit(KDE_BW_SENS_KM_RAW, "[,;|\\s]+"))))
} else {
  numeric(0)
}
KDE_BW_SENS_KM <- sort(unique(KDE_BW_SENS_KM[is.finite(KDE_BW_SENS_KM) & KDE_BW_SENS_KM > 0]))
# Fit C/D temporal-truncation sensitivity (independent of bandwidth/partition sens).
T_TRUNC_SENS_DAYS_RAW <- trimws(Sys.getenv("OK_T_TRUNC_SENS_DAYS", "1,5,7,10,14,21"))
T_TRUNC_SENS_DAYS <- suppressWarnings(as.numeric(unlist(strsplit(T_TRUNC_SENS_DAYS_RAW, "[,;|\\s]+"))))
T_TRUNC_SENS_DAYS <- sort(unique(T_TRUNC_SENS_DAYS[is.finite(T_TRUNC_SENS_DAYS) & T_TRUNC_SENS_DAYS > 0]))
RUN_T_TRUNC_SENSITIVITY <- tolower(Sys.getenv(
  "OK_RUN_T_TRUNC_SENSITIVITY",
  if (length(T_TRUNC_SENS_DAYS) > 0L) "true" else "false"
)) %in% c("1", "true", "yes", "y")
if (!RUN_T_TRUNC_SENSITIVITY) T_TRUNC_SENS_DAYS <- numeric(0)
RUN_FIT_VARIABILITY <- tolower(Sys.getenv("OK_RUN_FIT_VARIABILITY", "false")) %in% c("1", "true", "yes", "y")
FIT_VARIABILITY_REPS <- suppressWarnings(as.integer(Sys.getenv("OK_FIT_VARIABILITY_REPS", "")))
if (!is.finite(FIT_VARIABILITY_REPS) || is.na(FIT_VARIABILITY_REPS) || FIT_VARIABILITY_REPS < 1L) {
  FIT_VARIABILITY_REPS <- NA_integer_
}
FIT_VARIABILITY_CORES_RAW <- Sys.getenv("OK_FIT_VARIABILITY_CORES", "")
FIT_VARIABILITY_PATCH_FILE <- trimws(Sys.getenv("OK_FIT_VARIABILITY_PATCH_FILE", ""))
FIT_VARIABILITY_ONLY <- tolower(trimws(Sys.getenv("OK_FIT_VARIABILITY_ONLY", "false"))) %in%
  c("1", "true", "yes", "y")
BOOTSTRAP_PATCH_FILE <- trimws(Sys.getenv("OK_BOOTSTRAP_PATCH_FILE", ""))
BOOTSTRAP_ONLY <- tolower(trimws(Sys.getenv("OK_BOOTSTRAP_ONLY", "false"))) %in%
  c("1", "true", "yes", "y")
T_TRUNC_SENS_PATCH_FILE <- trimws(Sys.getenv("OK_T_TRUNC_SENS_PATCH_FILE", ""))
T_TRUNC_SENS_ONLY <- tolower(trimws(Sys.getenv("OK_T_TRUNC_SENS_ONLY", "false"))) %in%
  c("1", "true", "yes", "y")
# True C/D-only scope: skip homogeneous A/B and univariate G–J fits/ATEs.
CD_ONLY <- tolower(trimws(Sys.getenv("OK_CD_ONLY", "false"))) %in%
  c("1", "true", "yes", "y")
AB_ONLY <- tolower(trimws(Sys.getenv("OK_AB_ONLY", "false"))) %in%
  c("1", "true", "yes", "y")
# Lean univariate+KDE (public G/H): skip A–F and homogeneous I/J.
UNIV_KDE_ONLY <- tolower(trimws(Sys.getenv("OK_UNIV_KDE_ONLY", "false"))) %in%
  c("1", "true", "yes", "y")
SKIP_ATE <- tolower(trimws(Sys.getenv("OK_SKIP_ATE", "false"))) %in%
  c("1", "true", "yes", "y")
# Diagnostic control snapshots do not initialize SEM. Skip them in C/D-only
# and univ-kde-only by default; set OK_SKIP_CONTROL_SNAPSHOTS=0 to keep them.
SKIP_CONTROL_SNAPSHOTS <- tolower(trimws(Sys.getenv(
  "OK_SKIP_CONTROL_SNAPSHOTS",
  if (isTRUE(CD_ONLY) || isTRUE(UNIV_KDE_ONLY)) "true" else "false"
))) %in% c("1", "true", "yes", "y")
if (sum(c(isTRUE(AB_ONLY), isTRUE(CD_ONLY), isTRUE(UNIV_KDE_ONLY))) > 1L) {
  stop("OK_AB_ONLY, OK_CD_ONLY, and OK_UNIV_KDE_ONLY are mutually exclusive.")
}
SMOKE_SEM_D_SEEDS <- suppressWarnings(as.integer(Sys.getenv("OK_SMOKE_SEM_D_SEEDS", "0")))
if (!is.finite(SMOKE_SEM_D_SEEDS) || is.na(SMOKE_SEM_D_SEEDS) || SMOKE_SEM_D_SEEDS < 0L) {
  SMOKE_SEM_D_SEEDS <- 0L
}
SMOKE_SEM_D_TRUNC <- suppressWarnings(as.numeric(Sys.getenv("OK_SMOKE_SEM_D_TRUNC", "3")))
if (!is.finite(SMOKE_SEM_D_TRUNC) || is.na(SMOKE_SEM_D_TRUNC) || SMOKE_SEM_D_TRUNC <= 0) {
  SMOKE_SEM_D_TRUNC <- 3
}
# Multi-seed Fit D smoke reuses the trunc-only skip path (no main fits / ATE / bootstrap).
if (SMOKE_SEM_D_SEEDS > 0L) {
  T_TRUNC_SENS_ONLY <- TRUE
}
KDE_VARIANT_MODE <- tolower(trimws(Sys.getenv("OK_KDE_VARIANT_MODE", "triple")))
if (!KDE_VARIANT_MODE %in% c("single", "triple")) KDE_VARIANT_MODE <- "triple"
RUN_KDE_PROFILE_SWEEP <- identical(KDE_VARIANT_MODE, "triple")
# Background KDE bandwidth. Default is isotropic Scott (bw.scott.iso).
# "scott" uses spatstat::bw.scott (anisotropic); "scott-iso" uses bw.scott.iso.
KDE_BW_METHOD_RAW <- tolower(trimws(Sys.getenv("OK_KDE_BW_METHOD", "scott-iso")))
KDE_BW_METHOD <- switch(
  KDE_BW_METHOD_RAW,
  "scott" =, "bw.scott" =, "scott_aniso" = "scott",
  "scott-iso" =, "scott_iso" =, "bw.scott.iso" = "scott-iso",
  "diggle" =, "bw.diggle" = "diggle",
  "diggle2" =, "2diggle" =, "2*diggle" =, "digglex2" = "diggle2",
  NA_character_
)
if (!nzchar(KDE_BW_METHOD_RAW)) KDE_BW_METHOD <- "scott-iso"
if (is.na(KDE_BW_METHOD)) {
  stop("OK_KDE_BW_METHOD must be one of: scott-iso (default), scott, diggle, diggle2; got: ",
       KDE_BW_METHOD_RAW)
}
OK_VERBOSE <- tolower(Sys.getenv("OK_VERBOSE", "false")) %in% c("1", "true", "yes", "y")
DF_VERBOSE <- tolower(Sys.getenv("OK_DF_VERBOSE", "false")) %in% c("1", "true", "yes", "y")
SEM_WORKER_LOGS <- tolower(Sys.getenv("OK_SEM_WORKER_LOGS", "true")) %in% c("1", "true", "yes", "y")
SEM_WORKER_LOG_VERBOSE <- tolower(Sys.getenv("OK_SEM_WORKER_LOG_VERBOSE", "false")) %in% c("1", "true", "yes", "y")
SEM_WORKER_LOG_SPLIT <- tolower(Sys.getenv("OK_SEM_WORKER_LOG_SPLIT", "false")) %in% c("1", "true", "yes", "y")

STRUCT_DEFAULTS <- list(c = 0.05, p = 1.2, D = 5.0, gamma = GAMMA_FIXED, q = 1.5)
# Structural terms are fixed downstream after first-half pre-treatment calibration.
# gamma is always fixed at 0 (not estimated); other structural terms may be profile-specific.
FIXED_STRUCTURAL <- NULL

ATE_N_SIMS    <- if (QUICK_CHECK) 5 else if (TEST_MODE) 20 else 40
env_ate_n_sims <- suppressWarnings(as.integer(Sys.getenv("OK_ATE_N_SIMS", "")))
if (!is.na(env_ate_n_sims) && env_ate_n_sims > 0L) {
  ATE_N_SIMS <- env_ate_n_sims
}
ATE_WINDOW_DAYS <- 100

RUN_BOOTSTRAP_ATE <- tolower(Sys.getenv("OK_RUN_BOOTSTRAP_ATE", "false")) %in% c("1", "true", "yes", "y")
BOOT_N_REPS <- suppressWarnings(as.integer(Sys.getenv("OK_BOOT_N_REPS", "0")))
if (!is.finite(BOOT_N_REPS) || is.na(BOOT_N_REPS) || BOOT_N_REPS < 0L) BOOT_N_REPS <- 0L
# Default to bootstrapping the all-free KDE pair (C and D; primary publication pair).
BOOT_TARGETS_RAW <- toupper(Sys.getenv("OK_BOOT_TARGETS", "C,D"))
BOOT_TARGETS <- unique(trimws(unlist(strsplit(BOOT_TARGETS_RAW, ","))))
BOOT_TARGETS <- BOOT_TARGETS[BOOT_TARGETS %in% c("C", "D")]
if (length(BOOT_TARGETS) < 1) BOOT_TARGETS <- c("C", "D")
# Sensitivity targets are wired directly from bootstrap targets.
SENS_TARGETS <- unique(BOOT_TARGETS)
if (!setequal(SENS_TARGETS, c("C", "D"))) {
  stop(sprintf(
    "Sensitivity/bootstrap target alignment requires C,D; got: %s",
    paste(SENS_TARGETS, collapse = ",")
  ))
}
# Parametric-bootstrap uncertainty requires refitting each simulated catalogue.
# "none" is allowed only when bootstrap is disabled (e.g. cd-primary / quick).
BOOT_REFIT_SCOPE <- tolower(trimws(Sys.getenv("OK_BOOT_REFIT_SCOPE", "full")))
if (!BOOT_REFIT_SCOPE %in% c("none", "partial", "full")) {
  stop("OK_BOOT_REFIT_SCOPE must be 'none', 'partial', or 'full'; got: ", BOOT_REFIT_SCOPE)
}
if (identical(BOOT_REFIT_SCOPE, "partial")) {
  warning("OK_BOOT_REFIT_SCOPE='partial' is not a distinct valid law; running full refits.")
  BOOT_REFIT_SCOPE <- "full"
}
BOOT_SEM_INNER_ITER <- suppressWarnings(as.integer(Sys.getenv("OK_BOOT_SEM_INNER_ITER", "100")))
if (!is.finite(BOOT_SEM_INNER_ITER) || is.na(BOOT_SEM_INNER_ITER) || BOOT_SEM_INNER_ITER < 1L) {
  BOOT_SEM_INNER_ITER <- 100L
}
BOOT_OUTER_CORES_RAW <- Sys.getenv("OK_BOOT_OUTER_CORES", "")
BOOT_SEED <- suppressWarnings(as.integer(Sys.getenv("OK_BOOT_SEED", "")))
OK_GLOBAL_SEED <- suppressWarnings(as.integer(Sys.getenv("OK_GLOBAL_SEED", "1")))
if ((!is.finite(BOOT_SEED) || is.na(BOOT_SEED)) && is.finite(OK_GLOBAL_SEED) && !is.na(OK_GLOBAL_SEED)) {
  BOOT_SEED <- as.integer(OK_GLOBAL_SEED)
}
OK_IDENTICAL_RANDOMNESS <- tolower(Sys.getenv("OK_IDENTICAL_RANDOMNESS", "false")) %in% c("1", "true", "yes", "y")
OK_BOOT_IDENTICAL_RANDOMNESS <- tolower(Sys.getenv("OK_BOOT_IDENTICAL_RANDOMNESS", "false")) %in% c("1", "true", "yes", "y")
OK_BOOT_GUARD_DEGENERATE <- tolower(Sys.getenv("OK_BOOT_GUARD_DEGENERATE", "true")) %in% c("1", "true", "yes", "y")
OK_ATE_USE_CRN <- tolower(Sys.getenv("OK_ATE_USE_CRN", "true")) %in% c("1", "true", "yes", "y")
OK_ATE_CRN_PAIR <- tolower(Sys.getenv("OK_ATE_CRN_PAIR", "true")) %in% c("1", "true", "yes", "y")
OK_ATE_CONDITIONAL_ON_PRE <- tolower(Sys.getenv("OK_ATE_CONDITIONAL_ON_PRE", "true")) %in% c("1", "true", "yes", "y")
OK_ATE_CRN_BASE <- suppressWarnings(as.integer(Sys.getenv("OK_ATE_CRN_BASE", "")))
# Full bivariate forward sim for ATE (default: all-control vs all-treated).
# OK_ATE_CONTRAST=all_or_nothing|observed
# Set OK_ATE_BIVARIATE=false only for legacy marginal/univariate ATE sims.
OK_ATE_BIVARIATE <- tolower(Sys.getenv("OK_ATE_BIVARIATE", "true")) %in% c("1", "true", "yes", "y")
OK_ATE_CONTRAST <- tolower(trimws(Sys.getenv("OK_ATE_CONTRAST", "all_or_nothing")))
if (!OK_ATE_CONTRAST %in% c("all_or_nothing", "observed")) {
  stop("OK_ATE_CONTRAST must be 'all_or_nothing' or 'observed', got: ", OK_ATE_CONTRAST)
}
# Optional scenario slug for parallel comparison jobs (e.g. univ_aon, biv_obs).
OK_ATE_SCENARIO <- tolower(gsub("[^a-z0-9_]+", "_", trimws(Sys.getenv("OK_ATE_SCENARIO", ""))))
OK_ATE_SCENARIO <- gsub("^_+|_+$", "", OK_ATE_SCENARIO)
OK_ATE_METHOD_LABEL <- paste0(
  if (isTRUE(OK_ATE_BIVARIATE)) "bivariate" else "univariate",
  "_", OK_ATE_CONTRAST
)
source(file.path(SCRIPT_DIR, "oklahoma_geometry.R"), local = FALSE)
PRIMARY_PARTITION <- oklahoma_normalize_primary_partition(
  Sys.getenv("OK_PRIMARY_PARTITION", "grid_1.0R"),
  quick_check = QUICK_CHECK
)
source(file.path(SCRIPT_DIR, "ate_bivariate.R"), local = FALSE)
cat(sprintf("ATE evaluation: bivariate=%s contrast=%s\n",
            OK_ATE_BIVARIATE, OK_ATE_CONTRAST))
ETAS_BRANCHING_MAX <- suppressWarnings(as.numeric(Sys.getenv("OK_ETAS_BRANCHING_MAX", "0.98")))
if (!is.finite(ETAS_BRANCHING_MAX) || is.na(ETAS_BRANCHING_MAX) ||
    ETAS_BRANCHING_MAX <= 0 || ETAS_BRANCHING_MAX >= 1) {
  ETAS_BRANCHING_MAX <- 0.98
}
# Finite t_trunc renormalizes Omori, so p > 0 is admissible. q > 1 remains
# required for the infinite-plane spatial kernel. Cascade stability is the
# fitted GR branching ratio / spectral radius (cap 0.98). Finite first-moment
# restrictions p > 2 and q > 1.5 are not imposed.
ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS <- FALSE
ETAS_P_LOWER_BOUND <- suppressWarnings(as.numeric(Sys.getenv(
  "OK_ETAS_P_LOWER_BOUND", "0"
)))
ETAS_Q_LOWER_BOUND <- suppressWarnings(as.numeric(Sys.getenv(
  "OK_ETAS_Q_LOWER_BOUND", "1"
)))
if (!is.finite(ETAS_P_LOWER_BOUND) || ETAS_P_LOWER_BOUND < 0) {
  ETAS_P_LOWER_BOUND <- 0
}
if (!is.finite(ETAS_Q_LOWER_BOUND) || ETAS_Q_LOWER_BOUND < 1) {
  ETAS_Q_LOWER_BOUND <- 1
}
BOOT_BRANCHING_MAX <- suppressWarnings(as.numeric(Sys.getenv(
  "OK_BOOT_BRANCHING_MAX", as.character(ETAS_BRANCHING_MAX)
)))
if (!is.finite(BOOT_BRANCHING_MAX) || is.na(BOOT_BRANCHING_MAX) ||
    BOOT_BRANCHING_MAX <= 0 || BOOT_BRANCHING_MAX >= 1) {
  BOOT_BRANCHING_MAX <- ETAS_BRANCHING_MAX
}
BOOT_MAX_PRE_EVENTS <- suppressWarnings(as.integer(Sys.getenv("OK_BOOT_MAX_PRE_EVENTS", "25000")))
if (!is.finite(BOOT_MAX_PRE_EVENTS) || is.na(BOOT_MAX_PRE_EVENTS) || BOOT_MAX_PRE_EVENTS < 100L) BOOT_MAX_PRE_EVENTS <- 25000L
BOOT_MAX_POST_EVENTS_PER_PROC <- suppressWarnings(as.integer(Sys.getenv("OK_BOOT_MAX_POST_EVENTS_PER_PROC", "30000")))
if (!is.finite(BOOT_MAX_POST_EVENTS_PER_PROC) || is.na(BOOT_MAX_POST_EVENTS_PER_PROC) || BOOT_MAX_POST_EVENTS_PER_PROC < 100L) BOOT_MAX_POST_EVENTS_PER_PROC <- 30000L
BOOT_MAX_TOTAL_EVENTS <- suppressWarnings(as.integer(Sys.getenv("OK_BOOT_MAX_TOTAL_EVENTS", "90000")))
if (!is.finite(BOOT_MAX_TOTAL_EVENTS) || is.na(BOOT_MAX_TOTAL_EVENTS) || BOOT_MAX_TOTAL_EVENTS < 1000L) BOOT_MAX_TOTAL_EVENTS <- 90000L
REPORT_FORMATS_RAW <- tolower(trimws(Sys.getenv("OK_REPORT_FORMATS", "html,pdf")))
REPORT_FORMATS <- unique(trimws(unlist(strsplit(REPORT_FORMATS_RAW, ","))))
REPORT_FORMATS <- REPORT_FORMATS[REPORT_FORMATS %in% c("html", "pdf")]
MEMORY_SAFE <- tolower(Sys.getenv("OK_MEMORY_SAFE", "true")) %in% c("1", "true", "yes", "y")
TRIM_SENS_OBJECTS <- tolower(Sys.getenv("OK_TRIM_SENS_OBJECTS", if (MEMORY_SAFE) "true" else "false")) %in% c("1", "true", "yes", "y")
SENS_CORES_RAW <- Sys.getenv("OK_SENS_CORES", "")
ATE_SIM_CORES_RAW <- Sys.getenv("OK_ATE_SIM_CORES", "")

etas_names <- c("mu", "A", "alpha_m", "c", "p", "D", "gamma", "q")

# Local-first parallelization:
# - defaults to all local cores minus one
# - can be overridden via OK_CORES
default_local_cores <- max(1L, parallel::detectCores() - 1L)
OK_CORES_RAW <- Sys.getenv("OK_CORES", "")
N_CORES <- as.integer(ifelse(nzchar(OK_CORES_RAW), OK_CORES_RAW, default_local_cores))
N_CORES <- max(1L, min(N_CORES, parallel::detectCores()))
if (MEMORY_SAFE && !TEST_MODE && !QUICK_CHECK && !nzchar(OK_CORES_RAW)) {
  N_CORES <- min(N_CORES, 8L)
}
if (!is.na(SEM_PILOT_CORES)) {
  SEM_PILOT_CORES <- max(1L, min(SEM_PILOT_CORES, N_CORES))
}
# Bootstrap replicates are mostly single-threaded (SEM + nested ATE at n_cores=1),
# so default to one outer worker per allocated core. Override with OK_BOOT_OUTER_CORES
# or OK_BOOT_OUTER_CAP_MEMSAFE if a host is memory-constrained.
BOOT_OUTER_DEFAULT <- N_CORES
BOOT_OUTER_CORES <- suppressWarnings(as.integer(ifelse(nzchar(BOOT_OUTER_CORES_RAW), BOOT_OUTER_CORES_RAW, as.character(BOOT_OUTER_DEFAULT))))
if (!is.finite(BOOT_OUTER_CORES) || is.na(BOOT_OUTER_CORES) || BOOT_OUTER_CORES < 1L) BOOT_OUTER_CORES <- 1L
BOOT_OUTER_CORES <- max(1L, min(BOOT_OUTER_CORES, N_CORES))
if (RUN_BOOTSTRAP_ATE && BOOT_N_REPS > 0L && BOOT_OUTER_CORES > 1L) {
  boot_outer_cap_raw <- trimws(Sys.getenv("OK_BOOT_OUTER_CAP_MEMSAFE", ""))
  if (nzchar(boot_outer_cap_raw)) {
    boot_outer_cap <- suppressWarnings(as.integer(boot_outer_cap_raw))
    if (is.finite(boot_outer_cap) && !is.na(boot_outer_cap) && boot_outer_cap >= 1L &&
        BOOT_OUTER_CORES > boot_outer_cap) {
      cat(sprintf(
        "Bootstrap outer-core cap: reducing BOOT_OUTER_CORES from %d to %d (OK_BOOT_OUTER_CAP_MEMSAFE)\n",
        BOOT_OUTER_CORES, boot_outer_cap
      ))
      BOOT_OUTER_CORES <- boot_outer_cap
    }
  }
}
SENS_CORES_DEFAULT <- if (MEMORY_SAFE) min(2L, N_CORES) else N_CORES
SENS_CORES <- suppressWarnings(as.integer(ifelse(nzchar(SENS_CORES_RAW), SENS_CORES_RAW, as.character(SENS_CORES_DEFAULT))))
if (!is.finite(SENS_CORES) || is.na(SENS_CORES) || SENS_CORES < 1L) SENS_CORES <- 1L
SENS_CORES <- max(1L, min(SENS_CORES, N_CORES))
# Overwritten after the sensitivity job list is built so leftover cores go to M-step.
SENS_BIV_N_THREADS <- SEM_BIV_N_THREADS
# Use all allocated cores for ATE forward sims by default (sims are independent;
# CRN seeds are per-replicate so parallel is safe). Override with OK_ATE_SIM_CORES.
ATE_SIM_CORES_DEFAULT <- N_CORES
ATE_SIM_CORES <- suppressWarnings(as.integer(ifelse(nzchar(ATE_SIM_CORES_RAW), ATE_SIM_CORES_RAW, as.character(ATE_SIM_CORES_DEFAULT))))
if (!is.finite(ATE_SIM_CORES) || is.na(ATE_SIM_CORES) || ATE_SIM_CORES < 1L) ATE_SIM_CORES <- N_CORES
ATE_SIM_CORES <- max(1L, min(ATE_SIM_CORES, N_CORES))
FIT_VARIABILITY_CORES <- suppressWarnings(as.integer(ifelse(
  nzchar(FIT_VARIABILITY_CORES_RAW), FIT_VARIABILITY_CORES_RAW, as.character(N_CORES)
)))
if (!is.finite(FIT_VARIABILITY_CORES) || is.na(FIT_VARIABILITY_CORES) || FIT_VARIABILITY_CORES < 1L) {
  FIT_VARIABILITY_CORES <- N_CORES
}
FIT_VARIABILITY_CORES <- max(1L, min(FIT_VARIABILITY_CORES, N_CORES))
if (!is.finite(FIT_VARIABILITY_REPS) || is.na(FIT_VARIABILITY_REPS) || FIT_VARIABILITY_REPS < 1L) {
  FIT_VARIABILITY_REPS <- FIT_VARIABILITY_CORES
}
if (isTRUE(FIT_VARIABILITY_ONLY)) {
  RUN_FIT_VARIABILITY <- TRUE
  if (!nzchar(FIT_VARIABILITY_PATCH_FILE)) {
    stop("OK_FIT_VARIABILITY_ONLY=1 requires OK_FIT_VARIABILITY_PATCH_FILE (path to existing oklahoma_results*.rds to patch).")
  }
  patch_check <- normalizePath(FIT_VARIABILITY_PATCH_FILE, winslash = "/", mustWork = FALSE)
  if (!file.exists(patch_check)) {
    stop("OK_FIT_VARIABILITY_ONLY: patch file does not exist: ", patch_check)
  }
  chk <- tryCatch(readRDS(patch_check), error = function(e) NULL)
  if (is.null(chk)) {
    stop("OK_FIT_VARIABILITY_ONLY: could not readRDS patch file: ", patch_check)
  }
}
exclusive_only_modes <- sum(c(
  isTRUE(BOOTSTRAP_ONLY),
  isTRUE(FIT_VARIABILITY_ONLY),
  isTRUE(T_TRUNC_SENS_ONLY)
))
if (exclusive_only_modes > 1L) {
  stop("OK_BOOTSTRAP_ONLY, OK_FIT_VARIABILITY_ONLY, and OK_T_TRUNC_SENS_ONLY are mutually exclusive.")
}
if (isTRUE(T_TRUNC_SENS_ONLY)) {
  RUN_T_TRUNC_SENSITIVITY <- TRUE
  if (length(T_TRUNC_SENS_DAYS) < 1L) {
    T_TRUNC_SENS_DAYS <- c(1, 5, 7, 10, 14, 21)
  }
  RUN_SENSITIVITY <- FALSE
  RUN_FIT_VARIABILITY <- FALSE
  RUN_BOOTSTRAP_ATE <- FALSE
  BOOT_N_REPS <- 0L
  if (!nzchar(T_TRUNC_SENS_PATCH_FILE)) {
    stop("OK_T_TRUNC_SENS_ONLY=1 requires OK_T_TRUNC_SENS_PATCH_FILE (path to existing oklahoma_results*.rds to patch).")
  }
  trunc_patch_check <- normalizePath(T_TRUNC_SENS_PATCH_FILE, winslash = "/", mustWork = FALSE)
  if (!file.exists(trunc_patch_check)) {
    stop("OK_T_TRUNC_SENS_ONLY: patch file does not exist: ", trunc_patch_check)
  }
  trunc_chk <- tryCatch(readRDS(trunc_patch_check), error = function(e) NULL)
  if (is.null(trunc_chk)) {
    stop("OK_T_TRUNC_SENS_ONLY: could not readRDS patch file: ", trunc_patch_check)
  }
  # Faster early KDE setup; only Fit D trunc re-fits are needed.
  KDE_VARIANT_MODE <- "single"
  RUN_KDE_PROFILE_SWEEP <- FALSE
  rm(trunc_chk)
}
if (isTRUE(BOOTSTRAP_ONLY)) {
  RUN_BOOTSTRAP_ATE <- TRUE
  RUN_SENSITIVITY <- FALSE
  RUN_FIT_VARIABILITY <- FALSE
  if (!nzchar(BOOTSTRAP_PATCH_FILE)) {
    stop("OK_BOOTSTRAP_ONLY=1 requires OK_BOOTSTRAP_PATCH_FILE (path to existing for_paper/oklahoma_results*.rds).")
  }
  boot_patch_check <- normalizePath(BOOTSTRAP_PATCH_FILE, winslash = "/", mustWork = FALSE)
  if (!file.exists(boot_patch_check)) {
    stop("OK_BOOTSTRAP_ONLY: patch file does not exist: ", boot_patch_check)
  }
  boot_chk <- tryCatch(readRDS(boot_patch_check), error = function(e) NULL)
  if (is.null(boot_chk)) {
    stop("OK_BOOTSTRAP_ONLY: could not readRDS patch file: ", boot_patch_check)
  }
  # Primary all-free pair is fits_named$C/D; accept legacy E/F or fitE/fitF aliases.
  has_C <- !is.null(boot_chk$fits_named$C$params) ||
    !is.null(boot_chk$fits_named$E$params) ||
    !is.null(boot_chk$fitE$params)
  has_D <- !is.null(boot_chk$fits_named$D$params) ||
    !is.null(boot_chk$fits_named$F$params) ||
    !is.null(boot_chk$fitF$params)
  if (!isTRUE(has_C) || !isTRUE(has_D)) {
    stop("OK_BOOTSTRAP_ONLY: patch RDS must contain fits_named$C/D$params (or legacy E/F / fitE/fitF).")
  }
  if (BOOT_N_REPS < 1L) {
    stop("OK_BOOTSTRAP_ONLY requires OK_BOOT_N_REPS >= 1.")
  }
  # Faster early KDE setup; main C/D (all_free) params come from the patch.
  KDE_VARIANT_MODE <- "single"
  RUN_KDE_PROFILE_SWEEP <- FALSE
  rm(boot_chk)
}
# After exclusive-mode overrides: active bootstrap cannot use refit_scope=none.
if (isTRUE(RUN_BOOTSTRAP_ATE) && BOOT_N_REPS > 0L && identical(BOOT_REFIT_SCOPE, "none")) {
  warning("OK_BOOT_REFIT_SCOPE='none' is incompatible with active bootstrap; using 'full'.")
  BOOT_REFIT_SCOPE <- "full"
}
PARALLEL_BACKEND <- tolower(trimws(Sys.getenv(
  "OK_PARALLEL_BACKEND",
  if (MEMORY_SAFE) "psock" else "fork"
)))
if (!PARALLEL_BACKEND %in% c("fork", "psock", "sequential")) PARALLEL_BACKEND <- "psock"

if (is.finite(OK_GLOBAL_SEED) && !is.na(OK_GLOBAL_SEED)) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(OK_GLOBAL_SEED)
}

derive_run_seed <- function(base_seed, label = "", offset = 0L) {
  seed_base <- suppressWarnings(as.integer(base_seed))
  if (!is.finite(seed_base) || is.na(seed_base)) return(NA_integer_)
  label_int <- utf8ToInt(as.character(label))
  if (length(label_int) < 1L) label_int <- 0L
  # Deterministic per-run seed using only explicit inputs (no PID/time entropy).
  seed_val <- seed_base +
    sum(label_int * seq_along(label_int)) +
    as.integer(offset)
  seed_val <- abs(as.integer(seed_val %% .Machine$integer.max))
  if (!is.finite(seed_val) || is.na(seed_val) || seed_val < 1L) seed_val <- 1L
  seed_val
}
RNG_STREAM_CALL_COUNTER <- 0L

parse_numeric_vector_env <- function(env_name, default_vals, lower = -Inf) {
  raw <- trimws(Sys.getenv(env_name, ""))
  if (!nzchar(raw)) return(as.numeric(default_vals))
  vals <- suppressWarnings(as.numeric(trimws(unlist(strsplit(raw, ",")))))
  vals <- vals[is.finite(vals) & !is.na(vals) & vals >= lower]
  if (length(vals) < 1L) return(as.numeric(default_vals))
  unique(as.numeric(vals))
}

run_parallel <- function(X, FUN, cores, label = "job") {
  n <- length(X)
  cores_use <- max(1L, min(as.integer(cores), as.integer(n)))
  t0 <- proc.time()[["elapsed"]]
  RNG_STREAM_CALL_COUNTER <<- RNG_STREAM_CALL_COUNTER + 1L
  run_seed <- derive_run_seed(OK_GLOBAL_SEED, label = label, offset = n + RNG_STREAM_CALL_COUNTER)
  cat(sprintf("  [parallel:%s] start: n=%d cores=%d backend=%s\n",
              label, n, cores_use, PARALLEL_BACKEND))
  if (n <= 1L || cores_use <= 1L || identical(PARALLEL_BACKEND, "sequential")) {
    if (is.finite(run_seed) && !is.na(run_seed)) set.seed(run_seed)
    out <- lapply(X, FUN)
    cat(sprintf("  [parallel:%s] done in %.1fs (sequential)\n",
                label, proc.time()[["elapsed"]] - t0))
    return(out)
  }
  if (identical(PARALLEL_BACKEND, "fork")) {
    if (is.finite(run_seed) && !is.na(run_seed)) set.seed(run_seed)
    out <- parallel::mclapply(X, FUN, mc.cores = cores_use, mc.set.seed = TRUE)
    cat(sprintf("  [parallel:%s] done in %.1fs (fork)\n",
                label, proc.time()[["elapsed"]] - t0))
    return(out)
  }
  # Stream worker stdout/stderr into the main job log for progress visibility.
  cl <- parallel::makePSOCKcluster(cores_use, outfile = "")
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  if (is.finite(run_seed) && !is.na(run_seed)) {
    parallel::clusterSetRNGStream(cl, iseed = run_seed)
  }
  if (exists("REPO_DIR", envir = .GlobalEnv, inherits = FALSE)) {
    parallel::clusterExport(cl, varlist = c("REPO_DIR"), envir = .GlobalEnv)
  }
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(spatstat)
      library(sf)
      library(tigris)
      library(data.table)
      library(dplyr)
      library(ggplot2)
      library(parallel)
    })
    if (exists("REPO_DIR", inherits = TRUE) && requireNamespace("pkgload", quietly = TRUE)) {
      try(pkgload::load_all(REPO_DIR, quiet = TRUE, export_all = TRUE, helpers = FALSE, attach_testthat = FALSE),
          silent = TRUE)
    }
    NULL
  })
  export_mode <- tolower(trimws(Sys.getenv(
    "OK_PSOCK_EXPORT_MODE",
    if (MEMORY_SAFE) "minimal" else "all"
  )))
  if (!export_mode %in% c("minimal", "all")) export_mode <- "minimal"
  verify_worker_symbols <- function(symbols) {
    if (length(symbols) < 1L) return(TRUE)
    checks <- tryCatch(
      parallel::clusterCall(
        cl,
        function(nms) unname(vapply(nms, exists, logical(1), inherits = TRUE)),
        nms = symbols
      ),
      error = function(e) NULL
    )
    if (is.null(checks)) return(FALSE)
    all(vapply(checks, function(x) length(x) == length(symbols) && all(x), logical(1)))
  }
  safe_worker_globals <- function(nms) {
    # Socket clusters and their guard environments cannot be serialized into
    # nested PSOCK workers. Bootstrap ATE calls run sequentially inside each
    # outer worker and therefore do not need the reusable main-process pool.
    setdiff(unique(nms), c("ate_cl_reuse", "ate_cl_reuse_guard"))
  }
  if (identical(export_mode, "all")) {
    worker_globals <- safe_worker_globals(ls(envir = .GlobalEnv, all.names = TRUE))
  } else {
    # Include transitive globals referenced by helper functions. Inspecting only
    # FUN misses values used inside helpers (notably bootstrap simulation state)
    # and can turn every PSOCK replicate into a caught worker-side failure.
    collect_global_dependencies <- function(root_fun) {
      queue <- list(root_fun)
      needed <- character(0)
      analysed_functions <- character(0)
      while (length(queue) > 0L) {
        current <- queue[[1L]]
        queue <- queue[-1L]
        globs <- tryCatch(
          codetools::findGlobals(current, merge = FALSE),
          error = function(e) list(variables = character(0), functions = character(0))
        )
        refs <- unique(c(globs$variables, globs$functions))
        refs <- refs[nzchar(refs)]
        refs <- refs[vapply(
          refs,
          function(nm) exists(nm, envir = .GlobalEnv, inherits = FALSE),
          logical(1)
        )]
        new_refs <- setdiff(refs, needed)
        needed <- unique(c(needed, new_refs))
        helper_names <- new_refs[vapply(
          new_refs,
          function(nm) is.function(get(nm, envir = .GlobalEnv, inherits = FALSE)),
          logical(1)
        )]
        helper_names <- setdiff(helper_names, analysed_functions)
        if (length(helper_names) > 0L) {
          queue <- c(
            queue,
            lapply(helper_names, get, envir = .GlobalEnv, inherits = FALSE)
          )
          analysed_functions <- unique(c(analysed_functions, helper_names))
        }
      }
      needed
    }
    worker_globals <- safe_worker_globals(collect_global_dependencies(FUN))
  }
  if (length(worker_globals) > 0L) {
    parallel::clusterExport(cl, varlist = worker_globals, envir = .GlobalEnv)
  }
  if (!verify_worker_symbols(worker_globals) && identical(export_mode, "minimal")) {
    warning(sprintf(
      "PSOCK preflight missing symbols for %s under minimal export; retrying with full global export.",
      label
    ))
    worker_globals <- safe_worker_globals(ls(envir = .GlobalEnv, all.names = TRUE))
    if (length(worker_globals) > 0L) {
      parallel::clusterExport(cl, varlist = worker_globals, envir = .GlobalEnv)
    }
    if (!verify_worker_symbols(worker_globals)) {
      warning(sprintf(
        "PSOCK preflight still missing symbols for %s after full export.",
        label
      ))
    }
  }
  out <- tryCatch(
    parallel::parLapply(cl, X, FUN),
    error = function(e) {
      warning(sprintf("PSOCK failed for %s; falling back to sequential: %s", label, e$message))
      lapply(X, FUN)
    }
  )
  cat(sprintf("  [parallel:%s] done in %.1fs (psock)\n",
              label, proc.time()[["elapsed"]] - t0))
  out
}

run_parallel_on_cluster <- function(cl, X, FUN, label = "job") {
  n <- length(X)
  t0 <- proc.time()[["elapsed"]]
  cl_size <- tryCatch(length(cl), error = function(e) NA_integer_)
  cat(sprintf("  [parallel:%s] start: n=%d cores=%d backend=psock(reuse)\n",
              label, n, ifelse(is.finite(cl_size), cl_size, 0L)))
  if (n <= 1L) {
    out <- lapply(X, FUN)
    cat(sprintf("  [parallel:%s] done in %.1fs (sequential)\n",
                label, proc.time()[["elapsed"]] - t0))
    return(out)
  }
  if (is.finite(OK_GLOBAL_SEED) && !is.na(OK_GLOBAL_SEED)) {
    # Advance reused-cluster RNG streams per call so repeated jobs do not replay.
    RNG_STREAM_CALL_COUNTER <<- RNG_STREAM_CALL_COUNTER + 1L
    seed_step <- derive_run_seed(OK_GLOBAL_SEED, label = label, offset = n + RNG_STREAM_CALL_COUNTER)
    parallel::clusterSetRNGStream(cl, iseed = seed_step)
  }
  out <- tryCatch(
    parallel::parLapply(cl, X, FUN),
    error = function(e) {
      warning(sprintf("PSOCK(reuse) failed for %s; falling back to sequential: %s", label, e$message))
      lapply(X, FUN)
    }
  )
  cat(sprintf("  [parallel:%s] done in %.1fs (psock(reuse))\n",
              label, proc.time()[["elapsed"]] - t0))
  out
}

pp_mode_env <- trimws(Sys.getenv("PP_MODE", ""))
mode_label <- if (nzchar(pp_mode_env)) {
  toupper(pp_mode_env)
} else if (QUICK_CHECK) {
  "QUICK_CHECK"
} else if (TEST_MODE) {
  "TEST"
} else {
  "FULL"
}
cat("=== Oklahoma County-Based ETAS Analysis ===\n")
cat(sprintf("Mode: %s | SEM iters: %d | Change factor: %.3f | Cores: %d\n",
            mode_label, SEM_N_ITER, SEM_CHANGE_FACTOR,
            N_CORES))
cat(sprintf("SEM selection: method=%s | temperature=%.3f | monotone=%s | start_from_C=%s | single_flip_from_iter=%s | biv_threads=%d\n",
            SEM_OPTIM_METHOD, SEM_SELECTION_TEMPERATURE,
            SEM_MONOTONE_COMPLETE_LL, SEM_START_FROM_C,
            if (is.finite(SEM_SINGLE_FLIP_FROM_ITER)) {
              as.character(SEM_SINGLE_FLIP_FROM_ITER)
            } else {
              "off"
            },
            SEM_BIV_N_THREADS))
cat(sprintf("SEM change-factor bounds: min_mult=%.3f | max_mult=%.3f\n",
            SEM_CHANGE_FACTOR_MIN_MULT, SEM_CHANGE_FACTOR_MAX_MULT))
cat(sprintf("SEM relabel controls: max_step_frac=%.3f | force_param_update_flip_frac=%.3f\n",
            SEM_MAX_RELABEL_STEP_FRAC, SEM_FORCE_PARAM_UPDATE_FLIP_FRAC))
if (RUN_SEM_PILOT) {
  sem_pilot_cores_planned <- if (!is.na(SEM_PILOT_CORES)) SEM_PILOT_CORES else N_CORES
  cat(sprintf("SEM pilot: enabled | inner_iter=%d | cores=%d | max_combos=%d\n",
              SEM_PILOT_INNER_ITER, sem_pilot_cores_planned, SEM_PILOT_MAX_COMBOS))
} else {
  cat("SEM pilot: disabled\n")
}
cat(sprintf("Memory safe: %s | Sens cores: %d | ATE sim cores: %d | Boot outer cores: %d | Trim sensitivity objects: %s\n",
            MEMORY_SAFE, SENS_CORES, ATE_SIM_CORES, BOOT_OUTER_CORES, TRIM_SENS_OBJECTS))
cat(sprintf("ATE CRN: use=%s pair=%s conditional_on_pre=%s base=%s\n",
            OK_ATE_USE_CRN, OK_ATE_CRN_PAIR, OK_ATE_CONDITIONAL_ON_PRE,
            ifelse(is.finite(OK_ATE_CRN_BASE) && !is.na(OK_ATE_CRN_BASE), as.character(OK_ATE_CRN_BASE), "auto")))
cat(sprintf("Parallel backend: %s\n", PARALLEL_BACKEND))
cat(sprintf(
  "ETAS p/q constraints: finite moments=%s | normalization bounds p>%.3f q>%.3f | branching cap=%.3f\n",
  ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
  ETAS_P_LOWER_BOUND, ETAS_Q_LOWER_BOUND, ETAS_BRANCHING_MAX
))
cat(sprintf("SEM inner iters: main=%d, sensitivity=%d, bootstrap=%d\n",
            SEM_INNER_ITER, SENS_SEM_INNER_ITER, BOOT_SEM_INNER_ITER))
cat(sprintf("Targets (shared): sensitivity=%s | bootstrap=%s\n",
            paste(SENS_TARGETS, collapse = ","),
            paste(BOOT_TARGETS, collapse = ",")))
cat(sprintf("Fit variability stage: run=%s | reps=%d | cores=%d\n",
            RUN_FIT_VARIABILITY, FIT_VARIABILITY_REPS, FIT_VARIABILITY_CORES))
if (isTRUE(FIT_VARIABILITY_ONLY)) {
  cat(sprintf("Fit variability ONLY mode: will skip main county fits; patch into %s\n",
              FIT_VARIABILITY_PATCH_FILE))
}
if (isTRUE(BOOTSTRAP_ONLY)) {
  cat(sprintf(
    "Bootstrap ONLY mode: skip main county fits; hydrate C/D (all_free) from %s; bivariate=%s contrast=%s reps=%d scope=%s\n",
    BOOTSTRAP_PATCH_FILE, OK_ATE_BIVARIATE, OK_ATE_CONTRAST, BOOT_N_REPS, BOOT_REFIT_SCOPE
  ))
}
if (isTRUE(T_TRUNC_SENS_ONLY)) {
  if (SMOKE_SEM_D_SEEDS > 0L) {
    cat(sprintf(
      "Fit D SEM smoke: %d seeds at t_trunc=%.4g; SEM inner=%d; no ATE; patch source %s\n",
      SMOKE_SEM_D_SEEDS, SMOKE_SEM_D_TRUNC, SEM_INNER_ITER, T_TRUNC_SENS_PATCH_FILE
    ))
  } else {
    cat(sprintf(
      "t_trunc sensitivity ONLY mode: skip main county fits/ATEs/bootstrap; SEM inner=%d; grid=%s; patch into %s\n",
      SEM_INNER_ITER, paste(signif(T_TRUNC_SENS_DAYS, 4), collapse = ","), T_TRUNC_SENS_PATCH_FILE
    ))
  }
}
cat(sprintf("SEM warm-start fixed adaptive step: %s\n", SEM_WARMSTART_FIXED))
cat(sprintf("SEM proposal/event t_trunc (days): %s (source=%s)\n",
            as.character(signif(SEM_T_TRUNC_DAYS, 4)), SEM_T_TRUNC_SOURCE))
cat(sprintf("SEM temporal relabel weight: %.3f\n", SEM_TEMPORAL_WEIGHT))
cat(sprintf("B/D SEM verbose tracing: %s\n", DF_VERBOSE))
cat(sprintf("Verbose optimizer/SEM tracing: %s\n", OK_VERBOSE))
cat(sprintf("SEM worker logs enabled: %s | worker logs verbose: %s | split-to-main-log: %s\n",
            SEM_WORKER_LOGS, SEM_WORKER_LOG_VERBOSE, SEM_WORKER_LOG_SPLIT))
cat(sprintf("KDE variant mode: %s | bandwidth method: %s\n",
            KDE_VARIANT_MODE, KDE_BW_METHOD))
cat(sprintf("Primary tessellation: %s\n", PRIMARY_PARTITION))
if (length(KDE_BW_SENS_KM) > 0L) {
  cat(sprintf("KDE bandwidth sensitivity (km): %s | SEM inner=%d\n",
              paste(signif(KDE_BW_SENS_KM, 4), collapse = ", "),
              as.integer(SENS_SEM_INNER_ITER)))
}
cat(sprintf("Partition sensitivity: %s\n",
            if (isTRUE(RUN_PARTITION_SENSITIVITY)) "TRUE" else "FALSE"))
cat(sprintf("CD-only scope (skip A/B + univ): %s\n", if (isTRUE(CD_ONLY)) "TRUE" else "FALSE"))
cat(sprintf("AB-only scope (skip C/D + univ): %s\n", if (isTRUE(AB_ONLY)) "TRUE" else "FALSE"))
cat(sprintf("Univ-KDE-only scope (skip A–F + homog I/J): %s\n",
            if (isTRUE(UNIV_KDE_ONLY)) "TRUE" else "FALSE"))
cat(sprintf("Skip ATE evaluation: %s\n", if (isTRUE(SKIP_ATE)) "TRUE" else "FALSE"))
cat(sprintf("Skip control snapshots: %s\n", if (isTRUE(SKIP_CONTROL_SNAPSHOTS)) "TRUE" else "FALSE"))

analysis_start_time <- Sys.time()
analysis_start_elapsed <- proc.time()[["elapsed"]]
timing_rows <- list()
timing_row_idx <- 0L
add_timing_row <- function(stage, elapsed_sec, status = "ok", detail = NA_character_) {
  timing_row_idx <<- timing_row_idx + 1L
  timing_rows[[timing_row_idx]] <<- data.frame(
    order = timing_row_idx,
    stage = as.character(stage),
    elapsed_sec = as.numeric(elapsed_sec),
    elapsed_min = as.numeric(elapsed_sec) / 60,
    status = as.character(status),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}
mem_snapshot <- function() {
  gc_mb <- tryCatch(sum(gc(verbose = FALSE)[, 2], na.rm = TRUE), error = function(e) NA_real_)
  rss_mb <- tryCatch({
    status_path <- "/proc/self/status"
    if (!file.exists(status_path)) return(NA_real_)
    vmrss <- grep("^VmRSS:", readLines(status_path, warn = FALSE), value = TRUE)
    if (length(vmrss) < 1L) return(NA_real_)
    kb <- suppressWarnings(as.numeric(gsub("[^0-9]", "", vmrss[1])))
    if (!is.finite(kb) || is.na(kb)) NA_real_ else kb / 1024
  }, error = function(e) NA_real_)
  sprintf("rss=%.1fMB gc_heap=%.1fMB",
          ifelse(is.finite(rss_mb), rss_mb, NA_real_),
          ifelse(is.finite(gc_mb), gc_mb, NA_real_))
}
sanitize_log_tag <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) x <- "sem_fit"
  x
}
SEM_WORKER_LOG_DIR <- file.path(OUT_DIR, "worker_logs")
if (SEM_WORKER_LOGS && !dir.exists(SEM_WORKER_LOG_DIR)) {
  dir.create(SEM_WORKER_LOG_DIR, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================================
# 1. Data loading
# ============================================================================
cat("\n--- Step 1: Data loading ---\n")

ensure_oklahoma_catalog <- function(data_dir, etas_m0) {
  meta_path <- file.path(data_dir, "metadata.json")
  events_path <- file.path(data_dir, "events_all.csv")
  need_build <- !dir.exists(data_dir) || !file.exists(meta_path) || !file.exists(events_path)
  meta_min <- NA_real_
  if (!need_build) {
    meta_try <- tryCatch(
      jsonlite::fromJSON(readLines(meta_path)),
      error = function(e) NULL
    )
    if (!is.null(meta_try) && !is.null(meta_try$design$min_mag)) {
      meta_min <- suppressWarnings(as.numeric(meta_try$design$min_mag))
    }
    # Rebuild if catalog cutoff is coarser than the analysis m0 (e.g. old 2.5 vs new 2.0).
    if (!is.finite(meta_min) || meta_min > etas_m0 + 1e-9) {
      need_build <- TRUE
    }
  }
  if (need_build) {
    cat(sprintf(
      "Catalog missing or min_mag (%.3f) > ETAS_M0 (%.3f) — rebuilding with OK_FORCE_DOWNLOAD...\n",
      if (is.finite(meta_min)) meta_min else NA_real_, etas_m0
    ))
    Sys.setenv(OK_FORCE_DOWNLOAD = "true")
    source(file.path(SCRIPT_DIR, "Oklahoma_data_and_viz.R"), local = TRUE)
  }
  invisible(TRUE)
}
ensure_oklahoma_catalog(DATA_DIR, ETAS_M0)

ev_all <- fread(file.path(DATA_DIR, "events_all.csv"))
meta   <- jsonlite::fromJSON(readLines(file.path(DATA_DIR, "metadata.json")))

t_star_utc <- as.POSIXct(meta$design$t_star_utc, tz = "UTC")
ev_all[, time_utc := as.POSIXct(time_utc, tz = "UTC")]
ev_all[, t_days := as.numeric(difftime(time_utc, t_star_utc, units = "days"))]
ev_all[, x_km := x_m / 1000]
ev_all[, y_km := y_m / 1000]
if (!"mag" %in% names(ev_all)) stop("No 'mag' column in events data")
ev_all <- ev_all[mag >= ETAS_M0]

# BETA_GR is estimated only after exact Oklahoma support has been constructed
# and regional catalogue rows outside that support have been removed.

post_end_days <- as.numeric(difftime(
  as.POSIXct(meta$design$post_end_utc, tz = "UTC"), t_star_utc, units = "days"))

# ============================================================================
# 2. Geometry, AOI-exterior structural init, primary tessellation
# ============================================================================
cat("\n--- Step 2: Building study geometry ---\n")

options(tigris_use_cache = TRUE)
counties_sf <- counties(state = "OK", cb = TRUE, year = 2022)
counties_sf <- st_transform(counties_sf, CRS_PROJ)
counties_sf <- st_make_valid(counties_sf)

ok_boundary <- st_union(counties_sf)
ok_boundary <- st_make_valid(ok_boundary)

# The modeled support is the Oklahoma county union. Its bounding rectangle is
# substantially larger and must not be used for fitting or simulation.
win_km <- oklahoma_sf_to_owin_km(ok_boundary)

# Convert county polygons to spatstat owin tiles (preserve rings/multipart).
county_owins <- oklahoma_sf_features_to_owins_km(counties_sf, name_col = "NAME")
valid_idx <- !vapply(county_owins, is.null, logical(1))
if (sum(valid_idx) < 50) {
  cat("  Warning: only", sum(valid_idx), "of", nrow(counties_sf), "counties converted.\n")
}

county_owins_valid <- county_owins[valid_idx]
counties_sf_valid  <- counties_sf[valid_idx, ]
county_partition <- tess(tiles = county_owins_valid, window = win_km)

cat(sprintf("  Counties in tessellation: %d / %d\n",
            county_partition$n, nrow(counties_sf)))

# Assign treatment: county centroid inside OCC AOI (sensitivity / maps).
aoi_path <- file.path(DATA_DIR, "occ_aoi_layer_2.geojson")
aoi_sf <- st_read(aoi_path, quiet = TRUE)
aoi_sf <- st_transform(aoi_sf, CRS_PROJ)
aoi_sf <- st_make_valid(aoi_sf)
aoi_union <- st_union(aoi_sf)

county_centroids <- st_centroid(counties_sf_valid)
inside_aoi <- lengths(st_within(county_centroids, aoi_union)) > 0

county_processes <- ifelse(inside_aoi, "treated", "control")
names(county_processes) <- counties_sf_valid$NAME
county_treated_idx <- county_processes == "treated"
county_treated_names <- names(county_processes)[county_treated_idx]
county_control_ss <- as.owin(county_partition[!county_treated_idx])
county_treated_ss <- as.owin(county_partition[county_treated_idx])
county_state_spaces <- list(control = county_control_ss, treated = county_treated_ss)
county_scheme <- list(
  partition = county_partition,
  processes = county_processes,
  state_spaces = county_state_spaces,
  label = "county",
  treated_idx = county_treated_idx
)

cat(sprintf("  Treated counties: %d (%s)\n",
            sum(county_treated_idx),
            paste(county_treated_names, collapse = ", ")))
cat(sprintf("  Control counties: %d\n", sum(!county_treated_idx)))

aoi_owin <- tryCatch({
  oklahoma_sf_to_owin_km(aoi_union)
}, error = function(e) {
  cat("  AOI owin construction failed, using treated-county union as fallback\n")
  county_treated_ss
})
aoi_treat_ss <- tryCatch(
  intersect.owin(aoi_owin, win_km),
  error = function(e) county_treated_ss
)
aoi_ctrl_ss <- tryCatch(
  setminus.owin(win_km, aoi_treat_ss),
  error = function(e) county_control_ss
)

# Keep events inside the Oklahoma window; primary labels are assigned later.
keep_in_win <- function(df) {
  ok <- is.finite(df$x) & is.finite(df$y) & is.finite(df$t) &
    spatstat.geom::inside.owin(df$x, df$y, win_km)
  df[ok, , drop = FALSE]
}
attach_event_cols <- function(df) {
  df$W <- 1.0
  df$n <- nrow(df)
  df$background <- TRUE
  df$location_process <- NA_character_
  df
}
pp_pre  <- attach_event_cols(keep_in_win(as.data.frame(ev_all[t_days < 0,
  .(x = x_km, y = y_km, t = t_days, mag = mag)])))
pp_post <- attach_event_cols(keep_in_win(as.data.frame(ev_all[t_days >= 0 & t_days <= post_end_days,
  .(x = x_km, y = y_km, t = t_days, mag = mag)])))
cat(sprintf(
  "  Window support: pre=%d post=%d events inside Oklahoma\n",
  nrow(pp_pre), nrow(pp_post)
))

# Pre-treatment split:
# - hold out first 50% (earliest events) for KDE background estimation
# - use second 50% for model estimation / carry-over in SEM
pp_pre_all <- pp_pre[order(pp_pre$t), ]
n_pre_total <- nrow(pp_pre_all)

# Gutenberg-Richter beta from all supported pre-treatment events. Before the
# intervention the whole Oklahoma domain is the control process, so future
# tessellation treatment labels must not determine this sample.
if (!is.finite(BETA_GR) || BETA_GR <= 0) {
  mag_pre_ok <- pp_pre_all$mag
  mag_pre_ok <- mag_pre_ok[is.finite(mag_pre_ok) & mag_pre_ok >= ETAS_M0]
  if (length(mag_pre_ok) < 50L) {
    stop("Too few supported Oklahoma pre-treatment events to estimate BETA_GR (n=",
         length(mag_pre_ok), ").")
  }
  BETA_GR <- 1 / mean(mag_pre_ok - ETAS_M0)
  BETA_GR_SOURCE <- sprintf(
    "pre_treatment_oklahoma_MLE(n=%d, m0=%.3f)", length(mag_pre_ok), ETAS_M0
  )
}
if (!is.finite(BETA_GR) || BETA_GR <= 0) {
  stop("BETA_GR must be finite and positive; got ", BETA_GR)
}
cat(sprintf(
  "  Magnitude cutoff ETAS_M0=%.3f | BETA_GR=%.4f (b=%.3f) [%s]\n",
  ETAS_M0, BETA_GR, BETA_GR / log(10), BETA_GR_SOURCE
))
# Keep default starts / inits inside the GR productivity domain alpha < beta.
.alpha_cap <- max(0.05, BETA_GR - 0.05)
if (length(VANILLA_STARTS)) {
  VANILLA_STARTS <- lapply(VANILLA_STARTS, function(st) {
    if (!is.null(st$alpha_m) && is.finite(st$alpha_m)) {
      st$alpha_m <- min(as.numeric(st$alpha_m), .alpha_cap)
    }
    st
  })
}

n_pre_holdout <- floor(n_pre_total * 0.5)
if (n_pre_total > 0 && n_pre_holdout < 1) n_pre_holdout <- 1
holdout_idx <- if (n_pre_holdout > 0) seq_len(n_pre_holdout) else integer(0)
keep_idx <- if (n_pre_total > n_pre_holdout) (n_pre_holdout + 1):n_pre_total else integer(0)
if (length(keep_idx) < 1) {
  stop("Pre-treatment split left zero estimation events; cannot continue.")
}
pp_pre_holdout <- pp_pre_all[holdout_idx, ]
pp_pre <- pp_pre_all[keep_idx, ]

windowT_post <- c(0, post_end_days)
windowT_fit <- c(min(pp_pre$t, na.rm = TRUE), post_end_days)
# Structural decays from first-half *AOI-exterior* pre data (robust to the
# later treated-region rate ramp). The AOI polygon is the zero-background
# region for this calibration, not the primary treated tiles.
aoi_exterior <- !spatstat.geom::inside.owin(pp_pre_all$x, pp_pre_all$y, aoi_treat_ss)
pp_pre_ctrl_all <- pp_pre_all[aoi_exterior, , drop = FALSE]
pp_pre_ctrl_all <- pp_pre_ctrl_all[order(pp_pre_ctrl_all$t), ]
n_pre_ctrl_total <- nrow(pp_pre_ctrl_all)
n_pre_ctrl_struct_init <- floor(n_pre_ctrl_total * 0.5)
if (n_pre_ctrl_total > 0 && n_pre_ctrl_struct_init < 1) n_pre_ctrl_struct_init <- 1
ctrl_struct_idx <- if (n_pre_ctrl_struct_init > 0) seq_len(n_pre_ctrl_struct_init) else integer(0)
pp_pre_ctrl_struct_init <- pp_pre_ctrl_all[ctrl_struct_idx, ]
n_pre_ctrl_struct_init <- nrow(pp_pre_ctrl_struct_init)

if (n_pre_ctrl_struct_init < 5) {
  stop("Insufficient first-half AOI-exterior pre-treatment events for structural parameter estimation.")
}

# The full pre-treatment ETAS init used to warm-start every downstream fit
# is estimated separately on *all* pre-treatment events with whole-domain
# control background (pre-treatment everywhere is control), matching the
# tstar diagnostic snapshot.
estimate_structural_init <- function() {
  starts <- VANILLA_STARTS
  best <- NULL
  best_val <- -Inf
  wT_holdout <- c(min(pp_pre_ctrl_struct_init$t, na.rm = TRUE), max(pp_pre_ctrl_struct_init$t, na.rm = TRUE))
  if (!all(is.finite(wT_holdout)) || diff(wT_holdout) <= 0) {
    return(STRUCT_DEFAULTS)
  }
  for (s in starts) {
    fit_try <- tryCatch(
      fit_etas(
        params_init = s,
        realiz = pp_pre_ctrl_struct_init,
        windowT = wT_holdout,
        windowS = win_km,
        m0 = ETAS_M0,
        beta_gr = BETA_GR,
        enforce_finite_trigger_moments = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
        p_lower_bound = ETAS_P_LOWER_BOUND,
        q_lower_bound = ETAS_Q_LOWER_BOUND,
        max_branching_ratio = ETAS_BRANCHING_MAX,
        t_trunc = SEM_T_TRUNC_DAYS,
        maxit = VANILLA_MAXIT,
        fixed_params = with_gamma_fixed(NULL),
        zero_background_region = aoi_treat_ss
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_try) && is.finite(fit_try$value) && fit_try$value > best_val) {
      best <- fit_try
      best_val <- fit_try$value
    }
  }
  if (is.null(best) || is.null(best$par)) return(STRUCT_DEFAULTS)
  out <- as.list(best$par[c("c", "p", "D", "gamma", "q")])
  out <- out[!vapply(out, is.null, logical(1))]
  needed <- c("c", "p", "D", "gamma", "q")
  for (nm in needed) {
    if (is.null(out[[nm]]) || !is.finite(out[[nm]])) out[[nm]] <- STRUCT_DEFAULTS[[nm]]
  }
  out
}
STRUCT_INIT <- estimate_structural_init()
STRUCT_INIT$gamma <- GAMMA_FIXED
FIXED_STRUCTURAL <- with_gamma_fixed(as.list(STRUCT_INIT[c("c", "p", "D", "q")]))

estimate_pre_full_etas_init <- function() {
  df <- pp_pre_all[, c("x", "y", "t", "mag"), drop = FALSE]
  finite_rows <- is.finite(df$x) & is.finite(df$y) & is.finite(df$t) & is.finite(df$mag)
  df <- df[finite_rows, , drop = FALSE]
  df <- df[order(df$t), , drop = FALSE]
  if (nrow(df) < 5L) {
    return(list(mu = 1.0, A = 0.2, alpha_m = 0.8,
                c = STRUCT_INIT$c, p = STRUCT_INIT$p,
                D = STRUCT_INIT$D, gamma = GAMMA_FIXED, q = STRUCT_INIT$q))
  }
  w_start <- suppressWarnings(min(df$t, na.rm = TRUE))
  # End at treatment time (same window as the tstar control snapshot).
  wT <- c(w_start, 0)
  if (!all(is.finite(wT)) || diff(wT) <= 0) {
    return(list(mu = 1.0, A = 0.2, alpha_m = 0.8,
                c = STRUCT_INIT$c, p = STRUCT_INIT$p,
                D = STRUCT_INIT$D, gamma = GAMMA_FIXED, q = STRUCT_INIT$q))
  }
  start_from_struct <- list(
    mu = 1.0, A = 0.5, alpha_m = 0.5,
    c = STRUCT_INIT$c, p = STRUCT_INIT$p,
    D = STRUCT_INIT$D, gamma = GAMMA_FIXED, q = STRUCT_INIT$q
  )
  starts <- c(list(start_from_struct), VANILLA_STARTS)
  best <- NULL
  best_val <- -Inf
  for (s in starts) {
    fit_try <- tryCatch(
      fit_etas(
        params_init = s,
        realiz = df,
        windowT = wT,
        windowS = win_km,
        m0 = ETAS_M0,
        beta_gr = BETA_GR,
        enforce_finite_trigger_moments = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
        p_lower_bound = ETAS_P_LOWER_BOUND,
        q_lower_bound = ETAS_Q_LOWER_BOUND,
        max_branching_ratio = ETAS_BRANCHING_MAX,
        t_trunc = SEM_T_TRUNC_DAYS,
        maxit = VANILLA_MAXIT,
        # Pre-treatment everywhere is control: do not carve out the eventual
        # treated counties from the background support.
        fixed_params = with_gamma_fixed(NULL)
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_try) && is.finite(fit_try$value) && fit_try$value > best_val) {
      best <- fit_try
      best_val <- fit_try$value
    }
  }
  pre_needed <- c("mu", "A", "alpha_m", "c", "p", "D", "gamma", "q")
  if (is.null(best) || is.null(best$par)) {
    return(list(mu = 1.0, A = 0.2, alpha_m = 0.8,
                c = STRUCT_INIT$c, p = STRUCT_INIT$p,
                D = STRUCT_INIT$D, gamma = GAMMA_FIXED, q = STRUCT_INIT$q)[pre_needed])
  }
  pre_full <- as.list(best$par)
  for (nm in pre_needed) {
    if (is.null(pre_full[[nm]]) || !is.finite(pre_full[[nm]])) {
      pre_full[[nm]] <- if (nm %in% names(STRUCT_DEFAULTS)) STRUCT_DEFAULTS[[nm]] else VANILLA_STARTS[[1]][[nm]]
    }
  }
  pre_full$gamma <- GAMMA_FIXED
  pre_full[pre_needed]
}
PRE_CTRL_BOOT_PARAMS <- estimate_pre_full_etas_init()

cat(sprintf("  Structural init (from first 50%% AOI-exterior pre): c=%.4f, p=%.4f, D=%.4f, gamma=%.4f, q=%.4f\n",
            STRUCT_INIT$c, STRUCT_INIT$p, STRUCT_INIT$D,
            STRUCT_INIT$gamma, STRUCT_INIT$q))
cat(sprintf("  Pre-treatment full ETAS init (all pre, whole-domain control): mu=%.4f, A=%.4f, alpha_m=%.4f, c=%.4f, p=%.4f, D=%.4f, gamma=%.4f, q=%.4f\n",
            PRE_CTRL_BOOT_PARAMS$mu, PRE_CTRL_BOOT_PARAMS$A, PRE_CTRL_BOOT_PARAMS$alpha_m,
            PRE_CTRL_BOOT_PARAMS$c, PRE_CTRL_BOOT_PARAMS$p, PRE_CTRL_BOOT_PARAMS$D,
            PRE_CTRL_BOOT_PARAMS$gamma, PRE_CTRL_BOOT_PARAMS$q))
cat(sprintf("  SEM t_trunc: %.4f days (source=%s)\n", SEM_T_TRUNC_DAYS, SEM_T_TRUNC_SOURCE))

# Data-informed triggering range from the KDE holdout sample (Q90 NN).
# Needed before the primary grid tessellation is built.
X_bg_scale <- spatstat.geom::ppp(pp_pre_holdout$x, pp_pre_holdout$y, window = win_km)
nn_training_early <- spatstat.geom::nndist(X_bg_scale)
trigger_range_km <- as.numeric(stats::quantile(nn_training_early, probs = 0.9, na.rm = TRUE))
if (!is.finite(trigger_range_km) || trigger_range_km <= 0) {
  trigger_range_km <- STRUCT_DEFAULTS$D
}
cat(sprintf("  Estimated triggering range (Q90 NN distance): %.2f km\n", trigger_range_km))

build_grid_partition <- function(diameter, win, aoi_owin, label, max_tiles = NULL) {
  if (is.na(diameter)) {
    target_tiles <- 100L
    aspect <- diff(win$xrange) / diff(win$yrange)
    nx <- max(2L, round(sqrt(target_tiles * aspect)))
    ny <- max(2L, round(target_tiles / nx))
  } else {
    nx <- max(2L, ceiling(diff(win$xrange) / diameter))
    ny <- max(2L, ceiling(diff(win$yrange) / diameter))
    if (is.null(max_tiles)) {
      max_tiles <- if (TEST_MODE || QUICK_CHECK) 2000L else 5000L
    }
    if ((nx * ny) > max_tiles) {
      shrink <- sqrt((nx * ny) / max_tiles)
      nx <- max(2L, ceiling(nx / shrink))
      ny <- max(2L, ceiling(ny / shrink))
      cat(sprintf("  [%s] requested grid too fine; coarsened to %d x %d (<= %d tiles)\n",
                  label, nx, ny, max_tiles))
    }
  }
  grid_part <- spatstat.geom::quadrats(win, nx = nx, ny = ny)

  n_tiles <- grid_part$n
  grid_procs <- character(n_tiles)
  tile_list <- spatstat.geom::tiles(grid_part)
  for (i in seq_len(n_tiles)) {
    tile_cent <- c(mean(tile_list[[i]]$xrange), mean(tile_list[[i]]$yrange))
    grid_procs[i] <- if (spatstat.geom::inside.owin(tile_cent[1], tile_cent[2], aoi_owin)) {
      "treated"
    } else {
      "control"
    }
  }
  names(grid_procs) <- spatstat.geom::tilenames(grid_part)

  g_treated_idx <- grid_procs == "treated"
  g_ctrl_ss <- tryCatch(
    as.owin(grid_part[!g_treated_idx]),
    error = function(e) county_control_ss
  )
  g_treat_ss <- tryCatch(
    as.owin(grid_part[g_treated_idx]),
    error = function(e) county_treated_ss
  )
  g_state_spaces <- list(control = g_ctrl_ss, treated = g_treat_ss)

  cat(sprintf("  [%s] %d tiles (%d treated, %d control)\n",
              label, n_tiles, sum(g_treated_idx), sum(!g_treated_idx)))

  list(partition = grid_part, processes = grid_procs,
       state_spaces = g_state_spaces, label = label,
       treated_idx = g_treated_idx)
}

if (identical(PRIMARY_PARTITION, "county")) {
  primary_scheme <- county_scheme
} else if (identical(PRIMARY_PARTITION, "grid_coarse")) {
  primary_scheme <- build_grid_partition(NA, win_km, aoi_owin, "grid_coarse", 2000L)
} else {
  max_tiles_primary <- if (TEST_MODE || QUICK_CHECK) 2000L else 5000L
  primary_scheme <- build_grid_partition(
    diameter = trigger_range_km,
    win = win_km,
    aoi_owin = aoi_owin,
    label = "grid_1.0R",
    max_tiles = max_tiles_primary
  )
}

install_primary_scheme <- function(scheme) {
  partition <<- scheme$partition
  partition_processes <<- scheme$processes
  treated_idx <<- scheme$treated_idx
  treated_names <<- names(scheme$processes)[scheme$treated_idx]
  if (length(treated_names) < 1L) {
    treated_names <<- as.character(which(scheme$treated_idx))
  }
  control_ss <<- scheme$state_spaces$control
  treated_ss <<- scheme$state_spaces$treated
  state_spaces <<- scheme$state_spaces
  PRIMARY_PARTITION_LABEL <<- scheme$label
}
install_primary_scheme(primary_scheme)

relabel_keep <- function(df, context) {
  df$location_process <- oklahoma_assign_partition_process(
    df$x, df$y, partition, partition_processes
  )
  df$W <- 1.0
  df$n <- nrow(df)
  df$background <- TRUE
  out <- df[!is.na(df$location_process), , drop = FALSE]
  oklahoma_assert_label_support(out, control_ss, treated_ss, context = context)
  out
}
pp_pre_all <- relabel_keep(pp_pre_all, "pp_pre_all")
pp_pre_holdout <- relabel_keep(pp_pre_holdout, "pp_pre_holdout")
pp_pre <- relabel_keep(pp_pre, "pp_pre")
pp_post <- relabel_keep(pp_post, "pp_post")

# Carry-over convention: all pre-treatment events are control-process
# events regardless of eventual treated/control location status.
pp_pre$process  <- "control"
pp_post$process <- pp_post$location_process
pp_pre$inferred_process  <- "control"
pp_post$inferred_process <- pp_post$location_process
pp_all <- rbind(pp_pre, pp_post)
pp_all <- pp_all[order(pp_all$t), ]

if (identical(PRIMARY_PARTITION, "county")) {
  cat(sprintf("  Primary tessellation county: %d treated (%s)\n",
              sum(treated_idx), paste(treated_names, collapse = ", ")))
} else {
  cat(sprintf("  Primary tessellation %s: %d tiles (%d treated, %d control)\n",
              PRIMARY_PARTITION_LABEL, partition$n,
              sum(treated_idx), sum(!treated_idx)))
}

apply_pre_init_etas <- function(start_par, preserve = c("A", "mu")) {
  # Keep caller-provided productivity/background starts for multistart diversity;
  # pin remaining structural params to the pre-control bootstrap estimates.
  out <- start_par
  preserve <- intersect(preserve, etas_names)
  for (nm in etas_names) {
    if (nm %in% preserve) {
      v <- suppressWarnings(as.numeric(start_par[[nm]]))
      if (is.finite(v)) {
        out[[nm]] <- v
        next
      }
    }
    out[[nm]] <- PRE_CTRL_BOOT_PARAMS[[nm]]
  }
  out[["gamma"]] <- GAMMA_FIXED
  out
}

apply_pre_init_biv <- function(par_vec) {
  out <- par_vec
  out[c("mu_0", "mu_1")] <- PRE_CTRL_BOOT_PARAMS$mu
  # Keep cross-excitation starts weak but strictly positive to avoid
  # control-fixed variants getting stuck at zero interaction.
  out[c("A_00", "A_11")] <- PRE_CTRL_BOOT_PARAMS$A
  out[c("alpha_m_00", "alpha_m_11")] <- PRE_CTRL_BOOT_PARAMS$alpha_m
  diag_A <- suppressWarnings(as.numeric(out[c("A_00", "A_11")]))
  diag_alpha <- suppressWarnings(as.numeric(out[c("alpha_m_00", "alpha_m_11")]))
  if (!all(is.finite(diag_A)) || length(diag_A) < 2L) {
    diag_A <- rep(0.1, 2L)
  }
  if (!all(is.finite(diag_alpha)) || length(diag_alpha) < 2L) {
    diag_alpha <- rep(0.5, 2L)
  }
  cross_A_floor <- max(1e-3, 0.05 * mean(diag_A))
  cross_alpha_floor <- max(1e-3, 0.25 * mean(diag_alpha))
  for (nm in c("A_01", "A_10")) {
    v <- suppressWarnings(as.numeric(out[[nm]]))
    if (!is.finite(v) || is.na(v) || v <= 0) out[[nm]] <- cross_A_floor
  }
  for (nm in c("alpha_m_01", "alpha_m_10")) {
    v <- suppressWarnings(as.numeric(out[[nm]]))
    if (!is.finite(v) || is.na(v) || v <= 0) out[[nm]] <- cross_alpha_floor
  }
  out[c("c", "p", "D", "gamma", "q")] <- unlist(PRE_CTRL_BOOT_PARAMS[c("c", "p", "D", "gamma", "q")])
  out[["gamma"]] <- GAMMA_FIXED
  out
}

cat(sprintf("  Events in OK: pre=%d, post=%d, total=%d\n",
            nrow(pp_pre), nrow(pp_post), nrow(pp_all)))
cat(sprintf("  Pre-treatment split: holdout=%d (first 50%%), estimation=%d (second 50%%)\n",
            nrow(pp_pre_holdout), nrow(pp_pre)))
cat(sprintf("  Post-treatment: %d control, %d treated\n",
            sum(pp_post$location_process == "control"),
            sum(pp_post$location_process == "treated")))

# ----------------------------------------------------------------------------
# Control-parameter diagnostics over expanding horizons (report table only)
# ----------------------------------------------------------------------------
# These fits are diagnostic snapshots for understanding control dynamics.
# They do NOT alter SEM initialization or any downstream fitting machinery.
as_snapshot_df <- function(df) {
  out <- df[, c("x", "y", "t", "mag"), drop = FALSE]
  finite_rows <- is.finite(out$x) & is.finite(out$y) & is.finite(out$t) & is.finite(out$mag)
  out <- out[finite_rows, , drop = FALSE]
  out[order(out$t), , drop = FALSE]
}
fit_control_snapshot <- function(df, label, window_end) {
  if (is.null(df) || nrow(df) < 10L) {
    return(list(
      label = label,
      n_events = if (is.null(df)) 0L else nrow(df),
      window_start = NA_real_,
      window_end = as.numeric(window_end),
      status = "insufficient_data",
      loglik = NA_real_,
      params = NULL
    ))
  }
  starts <- c(list(PRE_CTRL_BOOT_PARAMS), VANILLA_STARTS)
  w_start <- suppressWarnings(min(df$t, na.rm = TRUE))
  w_end <- suppressWarnings(as.numeric(window_end))
  wT <- c(w_start, w_end)
  if (!all(is.finite(wT)) || diff(wT) <= 0) {
    return(list(
      label = label,
      n_events = nrow(df),
      window_start = w_start,
      window_end = w_end,
      status = "invalid_window",
      loglik = NA_real_,
      params = NULL
    ))
  }
  best <- NULL
  best_val <- -Inf
  for (s in starts) {
    fit_try <- tryCatch(
      fit_etas(
        params_init = s,
        realiz = df,
        windowT = wT,
        windowS = win_km,
        m0 = ETAS_M0,
        beta_gr = BETA_GR,
        enforce_finite_trigger_moments = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
        p_lower_bound = ETAS_P_LOWER_BOUND,
        q_lower_bound = ETAS_Q_LOWER_BOUND,
        max_branching_ratio = ETAS_BRANCHING_MAX,
        t_trunc = SEM_T_TRUNC_DAYS,
        maxit = VANILLA_MAXIT,
        fixed_params = with_gamma_fixed(NULL)
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_try) && is.finite(fit_try$value) && fit_try$value > best_val) {
      best <- fit_try
      best_val <- fit_try$value
    }
  }
  list(
    label = label,
    n_events = nrow(df),
    window_start = w_start,
    window_end = w_end,
    status = if (is.null(best) || is.null(best$par)) "failed" else "ok",
    loglik = if (is.null(best)) NA_real_ else best$value,
    params = if (is.null(best)) NULL else as.list(best$par)
  )
}
pp_all_control_full <- rbind(
  pp_pre_all[, c("x", "y", "t", "mag"), drop = FALSE],
  pp_post[, c("x", "y", "t", "mag"), drop = FALSE]
)
pp_all_control_full <- pp_all_control_full[order(pp_all_control_full$t), , drop = FALSE]
skipped_control_snapshot <- function(label) {
  list(
    label = label,
    n_events = NA_integer_,
    window_start = NA_real_,
    window_end = NA_real_,
    status = "skipped",
    loglik = NA_real_,
    params = NULL
  )
}
if (isTRUE(SKIP_CONTROL_SNAPSHOTS)) {
  cat("  Skipping diagnostic control snapshots (not used for SEM init).\n")
  control_snapshot_fits <- list(
    pre50 = skipped_control_snapshot("50% pretreatment (all control, all regions control)"),
    tstar = skipped_control_snapshot("Treatment time (all control, all regions control)"),
    end = skipped_control_snapshot("End time (all points treated as control, all regions control)")
  )
} else {
  control_snapshot_fits <- list(
    pre50 = fit_control_snapshot(
      as_snapshot_df(pp_pre_holdout),
      "50% pretreatment (all control, all regions control)",
      if (nrow(pp_pre_holdout) > 0) max(pp_pre_holdout$t, na.rm = TRUE) else NA_real_
    ),
    tstar = fit_control_snapshot(
      as_snapshot_df(pp_pre_all),
      "Treatment time (all control, all regions control)",
      0
    ),
    end = fit_control_snapshot(
      as_snapshot_df(pp_all_control_full),
      "End time (all points treated as control, all regions control)",
      post_end_days
    )
  )
}
for (nm in names(control_snapshot_fits)) {
  cs <- control_snapshot_fits[[nm]]
  csp <- cs$params
  cat(sprintf("  [control-snapshot:%s] status=%s n=%d window=[%.2f, %.2f] mu=%s A=%s alpha_m=%s\n",
              nm, cs$status, as.integer(cs$n_events),
              as.numeric(cs$window_start), as.numeric(cs$window_end),
              if (!is.null(csp$mu)) as.character(signif(csp$mu, 4)) else "NA",
              if (!is.null(csp$A)) as.character(signif(csp$A, 4)) else "NA",
              if (!is.null(csp$alpha_m)) as.character(signif(csp$alpha_m, 4)) else "NA"))
}

if (length(T_TRUNC_SENS_DAYS) > 0L) {
  cat(sprintf("  Fit C/D t_trunc sensitivity grid (days): %s | SEM inner=%d\n",
              paste(signif(T_TRUNC_SENS_DAYS, 4), collapse = ", "),
              as.integer(SENS_SEM_INNER_ITER)))
}

# ============================================================================
# 3. Partition and point pattern plots
# ============================================================================
cat("\n--- Step 3: Plots ---\n")

tryCatch({
  counties_plot <- counties_sf_valid
  counties_plot$treatment <- county_processes

  p_partition <- ggplot() +
    geom_sf(data = counties_plot, aes(fill = treatment),
            color = "grey40", linewidth = 0.4) +
    scale_fill_manual(values = c(control = "#deebf7", treated = "#fc9272"),
                      name = "Assignment") +
    geom_point(data = pp_post,
               aes(x = x * 1000, y = y * 1000),
               size = 0.5, alpha = 0.5) +
    labs(title = "Oklahoma County Partition",
         subtitle = "OCC directive AOI_20150318 (treated = red)",
         x = "Easting (m)", y = "Northing (m)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  store_analysis_plot("partition_map_county", p_partition)
  cat("  Stored county partition plot in results$plots$partition_map_county\n")
}, error = function(e) cat("  Partition plot error:", e$message, "\n"))

tryCatch({
  p_pre <- ggplot(pp_pre, aes(x = x, y = y, alpha = t)) +
    geom_point(color = "#377eb8", size = 1.5, shape = 21, fill = "#377eb8") +
    scale_alpha_continuous(name = "Time (days)", range = c(0.15, 0.8)) +
    labs(title = "Pre-treatment Events", x = "X (km)", y = "Y (km)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  store_analysis_plot("pp_pre_treatment", p_pre)
  cat("  Stored pre-treatment point plot in results$plots$pp_pre_treatment\n")

  pp_post_plot <- pp_post
  pp_post_plot$Process <- factor(pp_post_plot$location_process,
                                 levels = c("control", "treated"))
  p_post <- ggplot(pp_post_plot, aes(x = x, y = y, fill = Process, alpha = t)) +
    geom_point(size = 1.8, shape = 21, stroke = 0.3) +
    scale_fill_manual(values = c(control = "#377eb8", treated = "#e41a1c")) +
    scale_alpha_continuous(name = "Time (days)", range = c(0.3, 1)) +
    labs(title = "Post-treatment Events (Location Labels)",
         x = "X (km)", y = "Y (km)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  store_analysis_plot("pp_post_treatment", p_post)
  cat("  Stored post-treatment point plot in results$plots$pp_post_treatment\n")
}, error = function(e) cat("  PP plot error:", e$message, "\n"))

# ============================================================================
# 4A. Fit A: Naive independent ETAS  [DISABLED — kept for reference]
# ============================================================================
# Univariate fits are no longer run. The bivariate model is the focus.

fit_best_indep <- function(realiz, zbr, starts, maxit) {
  best_fit <- NULL; best_val <- -Inf
  starts_use <- lapply(starts, apply_pre_init_etas)
  for (s in starts_use) {
    fit <- tryCatch(
      fit_etas(params_init = s, realiz = realiz, windowT = windowT_post,
               windowS = win_km, m0 = ETAS_M0, maxit = maxit,
               fixed_params = with_gamma_fixed(NULL), zero_background_region = zbr,
               beta_gr = BETA_GR,
               enforce_finite_trigger_moments = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
               p_lower_bound = ETAS_P_LOWER_BOUND,
               q_lower_bound = ETAS_Q_LOWER_BOUND,
               max_branching_ratio = ETAS_BRANCHING_MAX,
               t_trunc = SEM_T_TRUNC_DAYS,
               log_transform = TRUE,
               soft_branching_barrier = TRUE,
               polish_productivity = TRUE,
               interior_restart = TRUE),
      error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$value) && fit$value > best_val) {
      best_fit <- fit; best_val <- fit$value
    }
  }
  best_fit
}

# Quick independent fits used only to initialize bivariate parameters
cat("\n--- Quick independent fits (for bivariate initialization) ---\n")
naive_control <- pp_post[pp_post$location_process == "control", ]
naive_treated <- pp_post[pp_post$location_process == "treated", ]

if (!TEST_MODE && !QUICK_CHECK && N_CORES > 1) {
  cat("  Running independent control/treated initial fits in parallel...\n")
  init_jobs <- list(
    list(id = "control", realiz = naive_control, zbr = treated_ss),
    list(id = "treated", realiz = naive_treated, zbr = control_ss)
  )
  init_out <- run_parallel(init_jobs, function(j) {
    list(id = j$id, fit = fit_best_indep(j$realiz, j$zbr, VANILLA_STARTS, VANILLA_MAXIT))
  }, cores = min(2L, N_CORES), label = "indep-init")
  fitA_ctrl <- init_out[[which(vapply(init_out, function(z) z$id == "control", logical(1)))]]$fit
  fitA_treat <- init_out[[which(vapply(init_out, function(z) z$id == "treated", logical(1)))]]$fit
} else {
  cat("  Fitting control...\n")
  fitA_ctrl <- fit_best_indep(naive_control, treated_ss, VANILLA_STARTS, VANILLA_MAXIT)
  cat("  Fitting treated...\n")
  fitA_treat <- fit_best_indep(naive_treated, control_ss, VANILLA_STARTS, VANILLA_MAXIT)
}

A_ctrl <- if (!is.null(fitA_ctrl)) as.list(fitA_ctrl$par) else VANILLA_STARTS[[1]]
A_treat <- if (!is.null(fitA_treat)) as.list(fitA_treat$par) else VANILLA_STARTS[[1]]
A_ctrl <- apply_pre_init_etas(A_ctrl)
A_treat <- apply_pre_init_etas(A_treat)
names(A_ctrl) <- etas_names; names(A_treat) <- etas_names

cat("  Control:", paste(etas_names, round(unlist(A_ctrl), 4), sep = "=", collapse = ", "), "\n")
cat("  Treated:", paste(etas_names, round(unlist(A_treat), 4), sep = "=", collapse = ", "), "\n")

# ============================================================================
# 4A. Fit A: Naive bivariate ETAS
# ============================================================================
cat("\n--- Fit A: Naive bivariate ETAS ---\n")

biv_init <- apply_pre_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
biv_names <- names(biv_init)
valid_biv_params <- function(x) {
  PPDisentangle:::.etas_biv_params_ok(
    x, BETA_GR, ETAS_BRANCHING_MAX, biv_names
  )
}
valid_biv_fit <- function(fit) {
  PPDisentangle:::.etas_biv_fit_ok(
    fit, BETA_GR, ETAS_BRANCHING_MAX, biv_names
  )
}
valid_biv_sem <- function(sem) {
  PPDisentangle:::.etas_biv_sem_ok(
    sem, BETA_GR, ETAS_BRANCHING_MAX, biv_names
  )
}
fit_b <- function() {
  tryCatch({
    fit_etas_bivariate(
      params_init = biv_init, realiz = pp_all,
      windowT = windowT_fit, windowS = win_km, m0 = ETAS_M0,
      control_state_space = control_ss, treated_state_space = treated_ss,
      treated_background_zero_before = 0,
      # Pre-treatment everywhere is control (flat background: area ratio).
      control_background_everywhere_before = 0,
      beta_gr = BETA_GR,
      enforce_finite_trigger_moments = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
      p_lower_bound = ETAS_P_LOWER_BOUND,
      q_lower_bound = ETAS_Q_LOWER_BOUND,
      max_branching_radius = ETAS_BRANCHING_MAX,
      # Homogeneous A: all free except gamma=0 (same freedom as C/D KDE all-free).
      maxit = VANILLA_MAXIT, fixed_params = with_gamma_fixed(NULL), trace = 0,
      t_trunc = SEM_T_TRUNC_DAYS
    )
  }, error = function(e) { cat("  Bivariate fit error:", e$message, "\n"); NULL })
}

# ============================================================================
# 4C. Fit C: SEM independent ETAS  [DISABLED — kept for reference]
# ============================================================================
# Univariate SEM is no longer run. Set placeholders for downstream code.
semC <- NULL
C_ctrl <- A_ctrl; C_treat <- A_treat

# ============================================================================
# 4B. Fit B: SEM bivariate ETAS
# ============================================================================
cat("\n--- Fit B: SEM bivariate ETAS ---\n")

biv_init_D <- apply_pre_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
# Homogeneous B: all free except gamma=0.
biv_fixed <- with_gamma_fixed(NULL)

enforce_control_decay_start <- function(params_vec, control_params = A_ctrl) {
  if (is.null(params_vec)) return(params_vec)
  decay_names <- c("c", "p", "D", "gamma", "q")
  out <- params_vec
  if (is.null(names(out))) return(out)
  ctrl_vals <- suppressWarnings(as.numeric(unlist(control_params[decay_names])))
  names(ctrl_vals) <- decay_names
  for (nm in decay_names) {
    if (nm %in% names(out) && nm %in% names(ctrl_vals) && is.finite(ctrl_vals[[nm]])) {
      out[[nm]] <- ctrl_vals[[nm]]
    }
  }
  out
}

run_sem_fit <- function(pp_data_in,
                        partition_in,
                        partition_processes_in,
                        state_spaces_in,
                        init_params_in,
                        init_ctrl_params_in = A_ctrl,
                        init_treat_params_in = A_treat,
                        model_type_in = "etas_bivariate",
                        fixed_params_in = with_gamma_fixed(biv_fixed),
                        background_rate_var_in = NULL,
                        use_pre_history_for_biv_in = TRUE,
                        treated_background_zero_before_in = 0,
                        control_background_everywhere_before_in = 0,
                        control_background_pre_mass_ratio_in = NULL,
                        sem_t_trunc_in = SEM_T_TRUNC_DAYS,
                        sem_inner_iter_in = SEM_INNER_ITER,
                        sem_change_factor_in = SEM_CHANGE_FACTOR,
                        sem_change_factor_min_mult_in = SEM_CHANGE_FACTOR_MIN_MULT,
                        sem_change_factor_max_mult_in = SEM_CHANGE_FACTOR_MAX_MULT,
                        sem_optim_method_in = SEM_OPTIM_METHOD,
                        sem_selection_temperature_in = SEM_SELECTION_TEMPERATURE,
                        verbose_in = DF_VERBOSE,
                        label = "SEM",
                        sem_rng_label_in = NULL,
                        biv_n_threads_in = SEM_BIV_N_THREADS) {
  fixed_params_in <- with_gamma_fixed(fixed_params_in)
  t0 <- proc.time()[["elapsed"]]
  sem_rng_label <- if (!is.null(sem_rng_label_in)) {
    rl <- as.character(sem_rng_label_in)[1L]
    if (nzchar(rl)) rl else label
  } else {
    label
  }
  sem_verbose_effective <- isTRUE(verbose_in)
  sem_log_file <- NULL
  if (SEM_WORKER_LOGS) {
    sem_log_file <- file.path(
      SEM_WORKER_LOG_DIR,
      add_file_tag(sprintf("sem_%s_pid%d.log", sanitize_log_tag(label), Sys.getpid()))
    )
    if (SEM_WORKER_LOG_VERBOSE) sem_verbose_effective <- TRUE
  }
  cat(sprintf("  [%s] start (n=%d, pid=%d, mem=%s)\n",
              label, nrow(pp_data_in), Sys.getpid(), mem_snapshot()))
  if (!is.null(sem_log_file)) {
    cat(sprintf("  [%s] worker sem log: %s\n", label, sem_log_file))
  }
  sem_seed_base <- derive_run_seed(OK_GLOBAL_SEED, label = paste0("sem:", sem_rng_label))
  sem_seed_step <- 0L
  next_sem_seed <- function() {
    if (!is.finite(sem_seed_base) || is.na(sem_seed_base)) return(invisible(NULL))
    sem_seed_step <<- sem_seed_step + 1L
    set.seed(derive_run_seed(sem_seed_base, label = sem_rng_label, offset = sem_seed_step))
  }
  run_one_sem <- function(pp_data_sem, init_params_sem, fixed_params_sem, sem_label) {
    next_sem_seed()
    adaptive_SEM(
      pp_data = pp_data_sem, partition = partition_in,
      partition_processes = partition_processes_in,
      statespace = win_km, time_window = windowT_post, treatment_time = 0,
      hawkes_params_control = init_ctrl_params_in, hawkes_params_treated = init_treat_params_in,
      N_labellings = SEM_N_LABELLINGS, N_iter = SEM_N_ITER, verbose = sem_verbose_effective,
      model_type = model_type_in,
      adaptive_control = list(
        param_update_cadence = SEM_PARAM_UPDATE,
        proposal_update_cadence = 1,
        state_spaces = state_spaces_in,
        iter = sem_inner_iter_in, n_props = SEM_INNER_PROPS,
        change_factor = sem_change_factor_in,
        change_factor_min_mult = sem_change_factor_min_mult_in,
        change_factor_max_mult = sem_change_factor_max_mult_in,
        optim_method = sem_optim_method_in,
        selection_temperature = sem_selection_temperature_in,
        verbose = sem_verbose_effective,
        stagnation_trigger_every = SEM_STAGNATION_TRIGGER_EVERY,
        temporal_weight = SEM_TEMPORAL_WEIGHT,
        temporal_scale_days = SEM_TEMPORAL_SCALE_DAYS,
        update_starting_data = TRUE, include_starting_data = TRUE,
        include_starting_first_n = sem_inner_iter_in,
        max_relabel_step_frac = SEM_MAX_RELABEL_STEP_FRAC,
        force_param_update_flip_frac = SEM_FORCE_PARAM_UPDATE_FLIP_FRAC,
        update_control_params = TRUE, fixed_params = fixed_params_sem,
        require_monotone_complete_ll = SEM_MONOTONE_COMPLETE_LL,
        proposal_method = "simulation",
        single_flip_from_iter = SEM_SINGLE_FLIP_FROM_ITER,
        outer_maxit = SEM_OUTER_MAXIT, outer_maxit_biv = SEM_OUTER_MAXIT_BIV,
        biv_n_threads = max(1L, as.integer(biv_n_threads_in))
      ),
      m0 = ETAS_M0, beta_gr = BETA_GR,
      enforce_finite_trigger_moments = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
      p_lower_bound = ETAS_P_LOWER_BOUND,
      q_lower_bound = ETAS_Q_LOWER_BOUND,
      max_branching_ratio = ETAS_BRANCHING_MAX,
      max_branching_radius = ETAS_BRANCHING_MAX,
      hard_subcritical = TRUE,
      soft_branching_barrier = TRUE,
      log_transform = TRUE,
      init_branching_margin = 0.9,
      etas_bivariate_params = if (identical(model_type_in, "etas_bivariate")) init_params_sem else NULL,
      background_rate_var = background_rate_var_in,
      etas_use_filtration_history = TRUE,
      use_pre_history_for_biv = use_pre_history_for_biv_in,
      treated_background_zero_before = treated_background_zero_before_in,
      control_background_everywhere_before = control_background_everywhere_before_in,
      control_background_pre_mass_ratio = control_background_pre_mass_ratio_in,
      t_trunc = sem_t_trunc_in
    )
  }
  run_sem_body <- function() {
    tryCatch({
      pp_data_sem <- pp_data_in
      init_params_sem <- init_params_in
      if (identical(model_type_in, "etas_bivariate") && !is.null(init_params_sem) &&
          !isTRUE(SEM_START_FROM_C)) {
        init_params_sem <- enforce_control_decay_start(init_params_sem, init_ctrl_params_in)
      }
      if (SEM_WARMSTART_FIXED && identical(model_type_in, "etas_bivariate")) {
        fixed_all <- as.list(init_params_sem)
        cat(sprintf("  [%s] warm adaptive step with fixed pre-initialized parameters (mem=%s)...\n",
                    label, mem_snapshot()))
        warm_sem <- run_one_sem(pp_data_sem, init_params_sem, fixed_all, paste0(label, " warm"))
        if (!is.null(warm_sem) && !is.null(warm_sem$adaptive$adaptive_labelling)) {
          pp_data_sem <- warm_sem$adaptive$adaptive_labelling
        }
        if (!is.null(warm_sem) && !is.null(warm_sem$etas_bivariate_params)) {
          init_params_sem <- warm_sem$etas_bivariate_params
        }
      }
      run_one_sem(pp_data_sem, init_params_sem, fixed_params_in, label)
    }, error = function(e) {
      cat(sprintf("  [%s] error: %s\n", label, e$message))
      NULL
    })
  }
  out <- if (!is.null(sem_log_file)) {
    out_local <- NULL
    sem_log_con <- file(sem_log_file, open = "wt")
    sem_msg_con <- file(sem_log_file, open = "at")
    on.exit({
      try(sink(type = "message"), silent = TRUE)
      try(sink(), silent = TRUE)
      try(close(sem_msg_con), silent = TRUE)
      try(close(sem_log_con), silent = TRUE)
    }, add = TRUE)
    sink(sem_log_con, split = isTRUE(SEM_WORKER_LOG_SPLIT))
    sink(sem_msg_con, type = "message")
    cat(sprintf("[%s] sem-log-begin ts=%s pid=%d mem=%s verbose=%s\n",
                label, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), Sys.getpid(), mem_snapshot(),
                sem_verbose_effective))
    out_local <- run_sem_body()
    cat(sprintf("[%s] sem-log-end ts=%s pid=%d mem=%s status=%s\n",
                label, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), Sys.getpid(), mem_snapshot(),
                ifelse(is.null(out_local), "failed", "ok")))
    out_local
  } else {
    run_sem_body()
  }
  t1 <- proc.time()[["elapsed"]]
  cat(sprintf("  [%s] done in %.1fs (mem=%s)\n", label, t1 - t0, mem_snapshot()))
  out
}

fit_d <- function() {
  tryCatch({
    run_sem_fit(
      pp_data_in = pp_all,
      partition_in = partition,
      partition_processes_in = partition_processes,
      state_spaces_in = state_spaces,
      init_params_in = biv_init_D,
      background_rate_var_in = NULL,
      verbose_in = DF_VERBOSE,
      label = "Fit B"
    )
  }, error = function(e) { cat("  SEM-biv error:", e$message, "\n"); NULL })
}

fitB <- NULL
semD <- NULL

B_params <- NULL
B_loglik <- NA_real_
D_params <- NULL
D_ctrl <- NULL
D_treat <- NULL

# ============================================================================
# 4E. Non-parametric background rate from first 50% of control pre-treatment
# ============================================================================
cat("\n--- Step 4E: KDE background rate from held-out pre-treatment data ---\n")

# Pre-treatment is treated as control-process everywhere; for KDE background
# training, use the full held-out pre sample (50% by event-count split),
# not only events in control-labelled counties.
pp_pre_holdout_ctrl <- pp_pre_holdout[order(pp_pre_holdout$t), ]
n_pre_holdout_ctrl <- nrow(pp_pre_holdout_ctrl)
if (n_pre_holdout_ctrl < 2) {
  stop("Insufficient held-out pre-treatment events for KDE background estimation.")
}

cat(sprintf("  Held-out pre-treatment events for KDE: %d\n",
            n_pre_holdout_ctrl))

kde_training <- pp_pre_holdout_ctrl
X_bg <- ppp(x = kde_training$x, y = kde_training$y, window = win_km)
bw_diggle <- as.numeric(suppressWarnings(bw.diggle(X_bg)))
if (!is.finite(bw_diggle) || bw_diggle <= 0) {
  stop("bw.diggle failed; cannot form the KDE background.")
}
bw_scott_xy <- tryCatch(as.numeric(bw.scott(X_bg)), error = function(e) c(NA_real_, NA_real_))
bw_scott_iso <- tryCatch({
  if (exists("bw.scott.iso", mode = "function")) {
    as.numeric(bw.scott.iso(X_bg))
  } else {
    sqrt(prod(bw_scott_xy))
  }
}, error = function(e) sqrt(prod(bw_scott_xy)))
if (identical(KDE_BW_METHOD, "scott")) {
  if (length(bw_scott_xy) < 2L || any(!is.finite(bw_scott_xy)) || any(bw_scott_xy <= 0)) {
    stop("bw.scott failed; cannot form the Scott KDE background.")
  }
  bw_sigma <- bw_scott_xy
  kde_bw_label <- sprintf("Scott anisotropic (%.2f x %.2f km)", bw_sigma[[1]], bw_sigma[[2]])
} else if (identical(KDE_BW_METHOD, "scott-iso")) {
  if (!is.finite(bw_scott_iso) || bw_scott_iso <= 0) {
    stop("bw.scott.iso failed; cannot form the isotropic Scott KDE background.")
  }
  bw_sigma <- bw_scott_iso
  kde_bw_label <- sprintf("Scott isotropic (%.2f km)", bw_sigma)
} else if (identical(KDE_BW_METHOD, "diggle")) {
  bw_sigma <- bw_diggle
  kde_bw_label <- sprintf("Diggle CV (%.2f km)", bw_sigma)
} else {
  bw_sigma <- 2 * bw_diggle
  kde_bw_label <- sprintf("2 x Diggle CV (%.2f km)", bw_sigma)
}
bw_sigma_xy <- as.numeric(bw_sigma)
bw_sigma_scalar <- if (length(bw_sigma_xy) == 1L) {
  bw_sigma_xy[[1]]
} else {
  sqrt(prod(bw_sigma_xy))
}
lambda_im <- density(X_bg, sigma = bw_sigma, edge = TRUE, at = "pixels")
min_nz <- min(lambda_im$v[lambda_im$v > 0], na.rm = TRUE)
lambda_im$v[lambda_im$v <= 0] <- min_nz
cat(sprintf("  KDE bandwidth [%s]: %s (from %d held-out training events; 2*diggle=%.2f km)\n",
            KDE_BW_METHOD, kde_bw_label, n_pre_holdout_ctrl, 2 * bw_diggle))

# Full-domain / control-region background mass ratio. Pre-treatment the
# control background covers the whole domain (density-continuous), so each
# pre-cutoff day charges mu_0 * this ratio in the compensator.
kde_pre_mass_ratio <- function(covariate_im, ctrl_win) {
  total_mass <- spatstat.geom::integral.im(covariate_im)
  ctrl_mass <- spatstat.geom::integral.im(covariate_im[ctrl_win, drop = FALSE])
  if (!is.finite(total_mass) || !is.finite(ctrl_mass) || ctrl_mass <= 0) {
    stop("Cannot compute control-background pre-treatment mass ratio.")
  }
  max(1, total_mass / ctrl_mass)
}
CTRL_BG_PRE_MASS_RATIO <- kde_pre_mass_ratio(lambda_im, control_ss)
cat(sprintf("  Control-background pre-treatment KDE mass ratio: %.3f\n",
            CTRL_BG_PRE_MASS_RATIO))

# KDE background lookups normalized to spatial mean one on each observed
# region, but evaluable over the whole domain. Used by simulation-based
# estimands (ATE counterfactuals, bootstrap world sims) so simulated
# backgrounds follow the fitted KDE field instead of a uniform density.
make_kde_bg_lookup <- function(ref_win, covariate_im = lambda_im) {
  # Bind the image in this frame. PSOCK workers do not inherit .GlobalEnv, so
  # nested lookups of lambda_im fail with "object 'lambda_im' not found"
  # (job 8469262 t_trunc ATE).
  force(covariate_im)
  cov_in_window <- covariate_im[ref_win, drop = FALSE]
  total_mass <- spatstat.geom::integral.im(cov_in_window)
  target_area <- spatstat.geom::area(ref_win)
  if (!is.finite(total_mass) || total_mass <= 0 ||
      !is.finite(target_area) || target_area <= 0) {
    stop("Cannot normalize KDE background lookup on a state space.")
  }
  norm_factor <- target_area / total_mass
  im_local <- covariate_im
  force(norm_factor)
  force(im_local)
  function(x, y) {
    vals <- spatstat.geom::interp.im(im_local, x, y)
    vals[!is.finite(vals) | vals < 0] <- 0
    vals * norm_factor
  }
}
KDE_BG_LOOKUP <- list(
  control = make_kde_bg_lookup(control_ss),
  treated = make_kde_bg_lookup(treated_ss)
)

# Triggering range was computed in Step 2 from the same holdout sample
# (needed to build the primary grid tessellation).
if (!exists("trigger_range_km", inherits = FALSE) ||
    !is.finite(trigger_range_km) || trigger_range_km <= 0) {
  nn_training <- nndist(X_bg)
  trigger_range_km <- as.numeric(stats::quantile(nn_training, probs = 0.9, na.rm = TRUE))
  if (!is.finite(trigger_range_km) || trigger_range_km <= 0) {
    trigger_range_km <- STRUCT_DEFAULTS$D
  }
}
cat(sprintf("  Triggering range (Q90 NN distance): %.2f km\n", trigger_range_km))

kde_info <- list(
  bw_method = KDE_BW_METHOD,
  bw_label = kde_bw_label,
  bw_sigma = as.numeric(bw_sigma_scalar),
  bw_sigma_xy = as.numeric(bw_sigma_xy),
  bw_diggle = as.numeric(bw_diggle),
  bw_scott_xy = as.numeric(bw_scott_xy),
  bw_scott_iso = as.numeric(bw_scott_iso),
  n_training = n_pre_holdout_ctrl,
  n_pre_ctrl_holdout = n_pre_holdout_ctrl,
  n_pre_holdout = nrow(pp_pre_holdout),
  n_pre_used = nrow(pp_pre),
  trigger_range_km = trigger_range_km
)

normalize_bg_weights <- function(df_sub, win_sub, covariate_im, mark_name = "W") {
  if (nrow(df_sub) == 0) return(list(new_df = df_sub, norm = 0))
  cov_in_window <- covariate_im[win_sub, drop = FALSE]
  total_mass_raw <- integral.im(cov_in_window)
  target_area <- spatstat.geom::area(win_sub)
  norm_factor <- target_area / total_mass_raw
  vals_raw <- spatstat.geom::interp.im(covariate_im, df_sub$x, df_sub$y)
  vals_raw[is.na(vals_raw)] <- 0
  df_sub[[mark_name]] <- vals_raw * norm_factor
  min_val <- min(df_sub[[mark_name]][df_sub[[mark_name]] > 0], na.rm = TRUE)
  if (is.infinite(min_val) || is.na(min_val)) min_val <- 1e-9
  df_sub[[mark_name]][df_sub[[mark_name]] <= 0] <- min_val
  return(list(new_df = df_sub, norm = norm_factor))
}

build_background_weighted_data <- function(lambda_covariate, ctrl_ss, treat_ss) {
  bg_ctrl_res <- normalize_bg_weights(pp_post[pp_post$location_process == "control", ],
                                      ctrl_ss, lambda_covariate)
  bg_treat_res <- normalize_bg_weights(pp_post[pp_post$location_process == "treated", ],
                                       treat_ss, lambda_covariate)
  pp_post_bg_local <- rbind(bg_ctrl_res$new_df, bg_treat_res$new_df)
  pp_post_bg_local <- pp_post_bg_local[order(pp_post_bg_local$t), ]

  # Pre-treatment everywhere is control: every pre event carries the
  # control-normalized weight (valid anywhere on the KDE domain), so
  # pre-cutoff events in treated counties are usable control background.
  pp_pre_bg_res <- normalize_bg_weights(pp_pre, ctrl_ss, lambda_covariate)
  pp_all_bg_local <- rbind(pp_pre_bg_res$new_df, pp_post_bg_local)
  pp_all_bg_local <- pp_all_bg_local[order(pp_all_bg_local$t), ]
  list(pp_post_bg = pp_post_bg_local, pp_all_bg = pp_all_bg_local)
}

bg_ctrl_res <- normalize_bg_weights(pp_post[pp_post$location_process == "control", ],
                                     control_ss, lambda_im)
bg_treat_res <- normalize_bg_weights(pp_post[pp_post$location_process == "treated", ],
                                      treated_ss, lambda_im)
pp_post_bg <- rbind(bg_ctrl_res$new_df, bg_treat_res$new_df)
pp_post_bg <- pp_post_bg[order(pp_post_bg$t), ]

# Pre-treatment everywhere is control: every pre event carries the
# control-normalized KDE weight (the KDE is defined on the whole domain), so
# pre-cutoff events in treated counties are usable control background.
pp_pre_bg_res <- normalize_bg_weights(pp_pre, control_ss, lambda_im)
pp_all_bg <- rbind(pp_pre_bg_res$new_df, pp_post_bg)
pp_all_bg <- pp_all_bg[order(pp_all_bg$t), ]

cat(sprintf("  Background weights assigned: post=%d, all=%d\n",
            nrow(pp_post_bg), nrow(pp_all_bg)))

# Save a map of inhomogeneous background-rate estimate.
tryCatch({
  bg_df <- as.data.frame(lambda_im)
  names(bg_df)[3] <- "lambda"
  p_bg <- ggplot(bg_df, aes(x = x, y = y, fill = lambda)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_viridis_c(name = "KDE rate", option = "C") +
    coord_equal() +
    labs(title = "Inhomogeneous Background Rate (KDE)",
         subtitle = paste0("Held-out first 50% of pre-treatment events; ", kde_bw_label),
         x = "X (km)", y = "Y (km)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  store_analysis_plot("background_rate_kde", p_bg)
  cat("  Stored KDE background-rate plot in results$plots$background_rate_kde\n")
}, error = function(e) cat("  Background-rate plot error:", e$message, "\n"))

# Parameter profiles for primary-tessellation KDE bivariate fits.
# - all_free: free except gamma=0 (magnitude-independent spatial scale); letters C/D (PRIMARY)
# - productivity_free: all mu/A/alpha free; structural terms fixed from pre-treatment; E/F
# Internal objects still use fitE/semF / kde_variant_fits$E/$F for naive/SEM slots;
# public fits_named letters follow C–J below.
kde_variant_specs <- list(
  all_free = list(
    id = "all_free",
    label = "all free except gamma=0",
    fixed_params = with_gamma_fixed(NULL)
  ),
  productivity_free = list(
    id = "productivity_free",
    label = "mu/A/alpha free, structural fixed (gamma=0)",
    fixed_params = with_gamma_fixed(FIXED_STRUCTURAL)
  )
)
kde_primary_variant_id <- "all_free"
kde_variant_ids <- if (RUN_KDE_PROFILE_SWEEP) names(kde_variant_specs) else kde_primary_variant_id
kde_variant_fits <- list(E = list(), F = list())
# Sensitivity and bootstrap refits align with the all-free KDE pair (C/D).
SENSITIVITY_FIXED_PARAMS <- kde_variant_specs$all_free$fixed_params
# Keys E/F here are naive/SEM slots in kde_variant_fits; values are public letters.
KDE_FIT_LETTERS <- list(
  all_free = list(E = "C", F = "D"),
  productivity_free = list(E = "E", F = "F")
)
kde_fit_label <- function(fit_type, variant_id) {
  letter <- KDE_FIT_LETTERS[[variant_id]][[fit_type]]
  sprintf("Fit %s [%s]", letter, variant_id)
}
# KDE bandwidth SEM (all-free D) uses this RNG key so stochastic SEM matches main Fit D
# when sigma and data match (i.e. same bandwidth vs main primary KDE).
OK_BW_SEM_RNG_LABEL <- kde_fit_label("F", "all_free")
cat(sprintf(
  "  KDE primary-partition fit variants to run: %s\n",
  paste(kde_variant_ids, collapse = ", ")
))
cat(sprintf(
  "  Sensitivity fixed profile: %s\n",
  if (is.null(SENSITIVITY_FIXED_PARAMS) || length(SENSITIVITY_FIXED_PARAMS) < 1) "<none (all-free)>" else paste(names(SENSITIVITY_FIXED_PARAMS), collapse = ", ")
))

# Shared ATE machinery (Step 6 main ATEs and Step 6a fit-variability).
windowT_ate <- c(0, ATE_WINDOW_DAYS)
ate_crn_base_seed <- if (is.finite(OK_ATE_CRN_BASE) && !is.na(OK_ATE_CRN_BASE)) {
  as.integer(OK_ATE_CRN_BASE)
} else if (is.finite(OK_GLOBAL_SEED) && !is.na(OK_GLOBAL_SEED)) {
  as.integer(100000L + 1000L * OK_GLOBAL_SEED)
} else {
  100000L
}
extract_marginals <- function(biv_par) {
  if (is.null(biv_par)) return(NULL)
  c_val <- biv_par[["c"]]; if (!is.finite(c_val)) c_val <- STRUCT_DEFAULTS$c
  p_val <- biv_par[["p"]]; if (!is.finite(p_val)) p_val <- STRUCT_DEFAULTS$p
  D_val <- biv_par[["D"]]; if (!is.finite(D_val)) D_val <- STRUCT_DEFAULTS$D
  g_val <- GAMMA_FIXED
  q_val <- biv_par[["q"]]; if (!is.finite(q_val)) q_val <- STRUCT_DEFAULTS$q
  structural <- c(c = c_val, p = p_val, D = D_val, gamma = g_val, q = q_val)
  ctrl <- as.list(c(mu = biv_par[["mu_0"]], A = biv_par[["A_00"]],
    alpha_m = biv_par[["alpha_m_00"]], structural))
  treat <- as.list(c(mu = biv_par[["mu_1"]], A = biv_par[["A_11"]],
    alpha_m = biv_par[["alpha_m_11"]], structural))
  list(ctrl = ctrl, treat = treat)
}
# When the monotone gate is on, D must start from that job's C MLE. Otherwise
# the first M-step can leave C's basin while the gate only blocks downhill LL.
sem_init_from_c <- function(c_params, fallback_init,
                            fallback_ctrl = A_ctrl,
                            fallback_treat = A_treat) {
  if (!isTRUE(SEM_START_FROM_C) || is.null(c_params)) {
    return(list(init = fallback_init, ctrl = fallback_ctrl, treat = fallback_treat))
  }
  marg <- extract_marginals(c_params)
  list(
    init = c_params,
    ctrl = if (!is.null(marg)) marg$ctrl else fallback_ctrl,
    treat = if (!is.null(marg)) marg$treat else fallback_treat
  )
}
extract_marginals_indep <- function(indep_par) {
  if (is.null(indep_par) || is.null(indep_par$control) || is.null(indep_par$treated)) return(NULL)
  list(ctrl = indep_par$control, treat = indep_par$treated)
}
ate_estim_fast <- function(ctrl_pp, treat_pp, observed_data, label,
                           n_tiles_used = partition$n,
                           treated_idx_used = treated_idx,
                           filtration_history = NULL,
                           crn_base_seed = NA_integer_,
                           phase = "main_fit",
                           quiet = FALSE,
                           contrast = NULL,
                           covariate_lookup = NULL,
                           t_trunc = SEM_T_TRUNC_DAYS) {
  if (is.null(contrast) || !nzchar(as.character(contrast)[1])) {
    contrast <- OK_ATE_CONTRAST
  }
  contrast <- match.arg(as.character(contrast)[1], c("all_or_nothing", "observed"))
  phase_tag <- tolower(gsub("[^A-Za-z0-9]+", "_", as.character(phase)))
  if (!nzchar(phase_tag)) phase_tag <- "main_fit"
  phase_label <- switch(
    phase_tag,
    main_fit = "main-fit",
    sensitivity = "sensitivity",
    bootstrap = "bootstrap",
    phase_tag
  )
  ate_prefix <- sprintf("[ATE:%s/univariate/%s]", phase_label, contrast)
  if (!quiet) cat(sprintf("  %s Computing ATE for %s...\n", ate_prefix, label))
  tryCatch({
    safe_q <- function(x, probs) {
      x <- x[is.finite(x)]
      if (length(x) < 1L) return(rep(NA_real_, length(probs)))
      as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
    }
    summarize_process_stability <- function(par_obj, beta_gr = BETA_GR) {
      if (is.null(par_obj)) {
        return(list(mu = NA_real_, A = NA_real_, alpha_m = NA_real_,
                    beta_gr = beta_gr, beta_minus_alpha = NA_real_,
                    eta_exact = NA_real_, eta_proxy = NA_real_))
      }
      A_val <- suppressWarnings(as.numeric(par_obj$A))
      alpha_val <- suppressWarnings(as.numeric(par_obj$alpha_m))
      mu_val <- suppressWarnings(as.numeric(par_obj$mu))
      beta_minus_alpha <- beta_gr - alpha_val
      eta_exact <- if (is.finite(A_val) && is.finite(alpha_val) && is.finite(beta_gr) && beta_minus_alpha > 0) {
        A_val * beta_gr / beta_minus_alpha
      } else {
        Inf
      }
      eta_proxy <- if (is.finite(A_val) && is.finite(alpha_val) && is.finite(beta_gr) && beta_gr > 0) {
        A_val * exp(alpha_val / beta_gr)
      } else {
        NA_real_
      }
      list(
        mu = mu_val, A = A_val, alpha_m = alpha_val,
        beta_gr = beta_gr, beta_minus_alpha = beta_minus_alpha,
        eta_exact = eta_exact, eta_proxy = eta_proxy
      )
    }
    ate_event_warn_threshold <- suppressWarnings(as.numeric(Sys.getenv("OK_ATE_EXPLOSIVE_EVENT_THRESHOLD", "20000")))
    if (!is.finite(ate_event_warn_threshold) || is.na(ate_event_warn_threshold) || ate_event_warn_threshold <= 0) {
      ate_event_warn_threshold <- 20000
    }
    ate_flagged_frac_warn <- suppressWarnings(as.numeric(Sys.getenv("OK_ATE_EXPLOSIVE_FRACTION_WARN", "0.05")))
    if (!is.finite(ate_flagged_frac_warn) || is.na(ate_flagged_frac_warn) ||
        ate_flagged_frac_warn < 0 || ate_flagged_frac_warn > 1) {
      ate_flagged_frac_warn <- 0.05
    }
    ctrl_diag <- summarize_process_stability(ctrl_pp)
    treat_diag <- summarize_process_stability(treat_pp)
    ctrl_br <- ctrl_diag$eta_exact
    treat_br <- treat_diag$eta_exact
    if (!quiet) {
      cat(sprintf("    %s sim config: sims=%d cores=%d backend=%s eta(ctrl=%.3f, treat=%.3f)\n",
                  ate_prefix,
                  ATE_N_SIMS, max(1L, min(ATE_SIM_CORES, ATE_N_SIMS)),
                  PARALLEL_BACKEND, ctrl_br, treat_br))
      cat(sprintf("    %s fitted stability (ctrl): mu=%.4g A=%.4g alpha_m=%.4g beta-alpha=%.4g eta_exact=%.3f eta_proxy=%.3f\n",
                  ate_prefix, ctrl_diag$mu, ctrl_diag$A, ctrl_diag$alpha_m,
                  ctrl_diag$beta_minus_alpha, ctrl_diag$eta_exact, ctrl_diag$eta_proxy))
      cat(sprintf("    %s fitted stability (treat): mu=%.4g A=%.4g alpha_m=%.4g beta-alpha=%.4g eta_exact=%.3f eta_proxy=%.3f\n",
                  ate_prefix, treat_diag$mu, treat_diag$A, treat_diag$alpha_m,
                  treat_diag$beta_minus_alpha, treat_diag$eta_exact, treat_diag$eta_proxy))
    }
    if (any(!is.finite(c(ctrl_br, treat_br))) ||
        ctrl_br >= 1 || treat_br >= 1) {
      stop(sprintf(
        "%s has an explosive univariate ATE law (eta ctrl=%.3f, treat=%.3f).",
        label, ctrl_br, treat_br
      ))
    }
    boundary_tol <- 1e-4
    near_fit_boundary <- c(
      control = ctrl_br >= ETAS_BRANCHING_MAX - boundary_tol,
      treated = treat_br >= ETAS_BRANCHING_MAX - boundary_tol
    )
    if (any(near_fit_boundary)) {
      warning(sprintf(
        "%s is on the %.3f branching-margin boundary (eta ctrl=%.6f, treat=%.6f).",
        label, ETAS_BRANCHING_MAX, ctrl_br, treat_br
      ))
    }
    if (isTRUE(OK_ATE_CONDITIONAL_ON_PRE)) {
      if (is.null(filtration_history)) {
        pre_history <- pp_pre[, c("x", "y", "t", "mag"), drop = FALSE]
      } else {
        pre_history <- as.data.frame(filtration_history)[, c("x", "y", "t", "mag"), drop = FALSE]
      }
    } else {
      pre_history <- pp_pre[0, c("x", "y", "t", "mag"), drop = FALSE]
    }
    # Independent univariate treated fits have no pre-treatment treated
    # parents. Preserve that fitted law in every treated counterfactual.
    control_history <- pre_history
    treated_history <- pre_history[0, , drop = FALSE]
    flat_lookup <- function(x, y) rep(1, length(x))
    lookup_for <- function(process) {
      out <- if (is.list(covariate_lookup)) {
        covariate_lookup[[process]]
      } else {
        covariate_lookup
      }
      if (is.null(out)) flat_lookup else out
    }
    ctrl_lookup <- lookup_for("control")
    treat_lookup <- lookup_for("treated")
    ctrl_ref_area <- spatstat.geom::area(control_ss)
    treat_ref_area <- spatstat.geom::area(treated_ss)
    if (!is.finite(crn_base_seed) || is.na(crn_base_seed)) {
      if (is.finite(OK_ATE_CRN_BASE) && !is.na(OK_ATE_CRN_BASE)) {
        crn_base_seed <- as.integer(OK_ATE_CRN_BASE)
      } else if (is.finite(OK_GLOBAL_SEED) && !is.na(OK_GLOBAL_SEED)) {
        crn_base_seed <- as.integer(100000L + 1000L * OK_GLOBAL_SEED)
      } else {
        crn_base_seed <- as.integer(100000L)
      }
    }
    ate_label_slug <- gsub("[^A-Za-z0-9]+", "_", label)
    ate_label_slug <- gsub("^_+|_+$", "", ate_label_slug)
    if (!nzchar(ate_label_slug)) ate_label_slug <- "model"
    ate_parallel_label <- sprintf("ate-sim-%s", tolower(substr(ate_label_slug, 1L, 40L)))
    # Observed contrast under independent univariate margins:
    # left = control params on full window; right = ctrl on control support
    # + treat on treated support (no cross-excitation).
    sim_right_count <- function(seed_right) {
      if (identical(contrast, "all_or_nothing")) {
        if (!is.null(seed_right)) set.seed(seed_right)
        t_sim <- sim_etas(treat_pp, windowT_ate, windowS = win_km,
                          m0 = ETAS_M0, beta_gr = BETA_GR,
                          filtration = treated_history,
                          covariate_lookup = treat_lookup,
                          bg_ref_area = treat_ref_area,
                          t_trunc = t_trunc)
        return(length(t_sim$t))
      }
      if (!is.null(seed_right)) set.seed(seed_right)
      t_ctrl <- sim_etas(ctrl_pp, windowT_ate, windowS = control_ss,
                         m0 = ETAS_M0, beta_gr = BETA_GR,
                         filtration = control_history,
                         covariate_lookup = ctrl_lookup,
                         bg_ref_area = ctrl_ref_area,
                         t_trunc = t_trunc)
      if (!is.null(seed_right)) set.seed(as.integer(seed_right + 17L))
      t_treat <- sim_etas(treat_pp, windowT_ate, windowS = treated_ss,
                          m0 = ETAS_M0, beta_gr = BETA_GR,
                          filtration = treated_history,
                          covariate_lookup = treat_lookup,
                          bg_ref_area = treat_ref_area,
                          t_trunc = t_trunc)
      length(t_ctrl$t) + length(t_treat$t)
    }
    run_one_sim <- function(s) {
      s_int <- suppressWarnings(as.integer(s))
      if (!is.finite(s_int) || is.na(s_int)) s_int <- 1L
      if (isTRUE(OK_ATE_USE_CRN)) {
        seed_s <- as.integer(crn_base_seed + s_int)
        if (isTRUE(OK_ATE_CRN_PAIR)) {
          set.seed(seed_s)
          c_sim <- sim_etas(ctrl_pp, windowT_ate, windowS = win_km,
                            m0 = ETAS_M0, beta_gr = BETA_GR,
                            filtration = control_history,
                            covariate_lookup = ctrl_lookup,
                            bg_ref_area = ctrl_ref_area,
                            t_trunc = t_trunc)
          t_count <- sim_right_count(seed_s)
        } else {
          set.seed(seed_s)
          c_sim <- sim_etas(ctrl_pp, windowT_ate, windowS = win_km,
                            m0 = ETAS_M0, beta_gr = BETA_GR,
                            filtration = control_history,
                            covariate_lookup = ctrl_lookup,
                            bg_ref_area = ctrl_ref_area,
                            t_trunc = t_trunc)
          t_count <- sim_right_count(as.integer(seed_s + 1000000L))
        }
      } else {
        c_sim <- sim_etas(ctrl_pp, windowT_ate, windowS = win_km,
                          m0 = ETAS_M0, beta_gr = BETA_GR,
                          filtration = control_history,
                          covariate_lookup = ctrl_lookup,
                          bg_ref_area = ctrl_ref_area,
                          t_trunc = t_trunc)
        t_count <- sim_right_count(NULL)
      }
      c(c_count = length(c_sim$t), t_count = as.numeric(t_count))
    }
    sim_results <- if (N_CORES > 1 && ATE_N_SIMS > 1) {
      reusable_pool <- if (exists("ate_cl_reuse", inherits = TRUE)) {
        get("ate_cl_reuse", inherits = TRUE)
      } else {
        NULL
      }
      if (!is.null(reusable_pool)) {
        run_parallel_on_cluster(
          reusable_pool,
          as.list(seq_len(ATE_N_SIMS)),
          run_one_sim,
          label = ate_parallel_label
        )
      } else {
        run_parallel(
          as.list(seq_len(ATE_N_SIMS)), run_one_sim,
          cores = min(ATE_SIM_CORES, ATE_N_SIMS),
          label = ate_parallel_label
        )
      }
    } else {
      lapply(seq_len(ATE_N_SIMS), run_one_sim)
    }
    c_counts <- vapply(sim_results, function(z) as.numeric(z[["c_count"]]), numeric(1))
    t_counts <- vapply(sim_results, function(z) as.numeric(z[["t_count"]]), numeric(1))
    safe_p95 <- function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 1L) return(NA_real_)
      as.numeric(stats::quantile(x, 0.95, na.rm = TRUE, names = FALSE))
    }
    safe_max <- function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 1L) return(NA_real_)
      max(x, na.rm = TRUE)
    }
    if (!quiet) {
      cat(sprintf("    %s sim summary (%s): ctrl mean=%.1f p95=%.1f max=%.0f | treat mean=%.1f p95=%.1f max=%.0f\n",
                  ate_prefix,
                  label,
                  mean(c_counts, na.rm = TRUE),
                  safe_p95(c_counts),
                  safe_max(c_counts),
                  mean(t_counts, na.rm = TRUE),
                  safe_p95(t_counts),
                  safe_max(t_counts)))
      ctrl_q <- safe_q(c_counts, c(0, 0.5, 0.9, 0.95, 0.99, 1))
      treat_q <- safe_q(t_counts, c(0, 0.5, 0.9, 0.95, 0.99, 1))
      total_counts <- c_counts + t_counts
      total_q <- safe_q(total_counts, c(0, 0.5, 0.9, 0.95, 0.99, 1))
      cat(sprintf(
        "    %s sim counts ctrl[min,med,p90,p95,p99,max]=[%.0f, %.0f, %.0f, %.0f, %.0f, %.0f]\n",
        ate_prefix, ctrl_q[1], ctrl_q[2], ctrl_q[3], ctrl_q[4], ctrl_q[5], ctrl_q[6]
      ))
      cat(sprintf(
        "    %s sim counts treat[min,med,p90,p95,p99,max]=[%.0f, %.0f, %.0f, %.0f, %.0f, %.0f]\n",
        ate_prefix, treat_q[1], treat_q[2], treat_q[3], treat_q[4], treat_q[5], treat_q[6]
      ))
      cat(sprintf(
        "    %s sim counts total[min,med,p90,p95,p99,max]=[%.0f, %.0f, %.0f, %.0f, %.0f, %.0f]\n",
        ate_prefix, total_q[1], total_q[2], total_q[3], total_q[4], total_q[5], total_q[6]
      ))
    }
    flagged_mask <- (is.finite(c_counts) & c_counts > ate_event_warn_threshold) |
      (is.finite(t_counts) & t_counts > ate_event_warn_threshold)
    flagged_n <- sum(flagged_mask, na.rm = TRUE)
    if (flagged_n > 0L) {
      flagged_frac <- flagged_n / max(1L, length(c_counts))
      cat(sprintf("    [warning:ate-explosive] %s: %d/%d sims exceed %.0f events (frac=%.3f)\n",
                  label, flagged_n, length(c_counts), ate_event_warn_threshold, flagged_frac))
      if (flagged_frac >= ate_flagged_frac_warn) {
        cat(sprintf("    [warning:ate-explosive] %s: high fraction of extreme sims (>= %.3f); inspect params/branching.\n",
                    label, ate_flagged_frac_warn))
      }
    }
    if (!quiet) {
      top_n <- min(5L, length(c_counts))
      if (top_n > 0L) {
        ord <- order(c_counts + t_counts, decreasing = TRUE, na.last = TRUE)
        top_idx <- ord[seq_len(top_n)]
        top_desc <- paste(vapply(top_idx, function(ii) {
          sprintf("#%d:c=%d,t=%d,total=%d", ii, as.integer(c_counts[ii]), as.integer(t_counts[ii]), as.integer(c_counts[ii] + t_counts[ii]))
        }, character(1)), collapse = " | ")
        cat(sprintf("    %s top total-count sims: %s\n", ate_prefix, top_desc))
      }
    }
    total_saved <- c_counts - t_counts
    all_nothing_sim <- data.frame(
      c_total = c_counts,
      t_total = t_counts,
      total_saved = total_saved,
      total_effect = total_saved,
      c_mean = c_counts / n_tiles_used,
      t_mean = t_counts / n_tiles_used,
      saved_per_tile = total_saved / n_tiles_used,
      ATE = total_saved / n_tiles_used
    )
    n_ctrl_loc <- sum(observed_data$location_process == "control", na.rm = TRUE)
    n_treat_loc <- sum(observed_data$location_process == "treated", na.rm = TRUE)
    n_ctrl_tiles <- sum(!treated_idx_used)
    n_treat_tiles <- sum(treated_idx_used)
    saved_naive <- NA_real_
    saved_spillover <- NA_real_
    total_saved_naive <- NA_real_
    total_saved_spillover <- NA_real_
    if (n_ctrl_tiles > 0 && n_treat_tiles > 0) {
      saved_naive <- n_ctrl_loc / n_ctrl_tiles - n_treat_loc / n_treat_tiles
      total_saved_naive <- saved_naive * n_tiles_used
      saved_spillover <- if ("inferred_process" %in% names(observed_data)) {
        n_ctrl_loc / n_ctrl_tiles -
          sum(observed_data$inferred_process == "control" &
              observed_data$location_process == "control", na.rm = TRUE) / n_ctrl_tiles
      } else { 0 }
      total_saved_spillover <- saved_spillover * n_tiles_used
    }
    analytic <- ATE_analytic_etas(ctrl_pp, treat_pp,
                                  windowT = windowT_ate, n_tiles = n_tiles_used,
                                  beta_gr = BETA_GR, m0 = ETAS_M0)
    analytic_saved <- if (!is.null(analytic)) {
      c(
        eta_ctrl_minus_treat = analytic$eta_ctrl - analytic$eta_treat,
        total_ctrl_minus_treat = (analytic$eta_ctrl - analytic$eta_treat) * n_tiles_used
      )
    } else {
      c(eta_ctrl_minus_treat = NA_real_, total_ctrl_minus_treat = NA_real_)
    }
    list(all_nothing_sim = all_nothing_sim,
         saved_naive = saved_naive, saved_spillover = saved_spillover,
         total_saved_naive = total_saved_naive, total_saved_spillover = total_saved_spillover,
         ATE_naive = saved_naive, ATE_spillover = saved_spillover,
         total_naive = total_saved_naive, total_spillover = total_saved_spillover,
         n_tiles_used = n_tiles_used,
         treated_pp = treat_pp, control_pp = ctrl_pp,
         analytic = analytic,
         analytic_saved = analytic_saved,
         stability = list(
           eta_control = ctrl_br,
           eta_treated = treat_br,
           near_fit_boundary = near_fit_boundary,
           fit_margin = ETAS_BRANCHING_MAX
         ),
         ate_method = paste0("univariate_", contrast),
         contrast = contrast)
  }, error = function(e) { cat("    Error:", e$message, "\n"); NULL })
}
ate_cl_reuse <- NULL
ate_cl_reuse_guard <- new.env(parent = emptyenv())
ate_cl_reuse_guard$cluster <- NULL
reg.finalizer(
  ate_cl_reuse_guard,
  function(e) {
    if (!is.null(e$cluster)) {
      try(parallel::stopCluster(e$cluster), silent = TRUE)
      e$cluster <- NULL
    }
  },
  onexit = TRUE
)
stop_ate_psock_pool <- function() {
  if (!is.null(ate_cl_reuse)) {
    try(parallel::stopCluster(ate_cl_reuse), silent = TRUE)
    ate_cl_reuse <<- NULL
    ate_cl_reuse_guard$cluster <- NULL
  }
  invisible(NULL)
}
ensure_ate_psock_pool <- function() {
  if (!is.null(ate_cl_reuse)) return(invisible(NULL))
  ate_cl_reuse_cores <- max(1L, min(ATE_SIM_CORES, ATE_N_SIMS))
  if (PARALLEL_BACKEND == "psock" && N_CORES > 1 && ATE_N_SIMS > 1 && ate_cl_reuse_cores > 1L) {
    cat(sprintf("  [parallel:ate-sim] initializing reusable PSOCK pool: cores=%d\n", ate_cl_reuse_cores))
    ate_cl_reuse <<- parallel::makePSOCKcluster(ate_cl_reuse_cores, outfile = "")
    ate_cl_reuse_guard$cluster <- ate_cl_reuse
    if (exists("REPO_DIR", envir = .GlobalEnv, inherits = FALSE)) {
      parallel::clusterExport(ate_cl_reuse, varlist = c("REPO_DIR"), envir = .GlobalEnv)
    }
    parallel::clusterEvalQ(ate_cl_reuse, {
      suppressPackageStartupMessages({
        library(spatstat)
        library(sf)
        library(tigris)
        library(data.table)
        library(dplyr)
        library(ggplot2)
        library(parallel)
      })
      if (exists("REPO_DIR", inherits = TRUE) && requireNamespace("pkgload", quietly = TRUE)) {
        try(pkgload::load_all(REPO_DIR, quiet = TRUE, export_all = TRUE, helpers = FALSE, attach_testthat = FALSE),
            silent = TRUE)
      }
      NULL
    })
  }
  invisible(NULL)
}

# ============================================================================
# 4C. Fit C: Naive bivariate ETAS with KDE background
# ============================================================================
cat("\n--- Fit C: Naive bivariate ETAS with KDE background ---\n")

biv_init_E <- apply_pre_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
fit_e <- function(init_params = biv_init_E,
                  fixed_params = NULL,
                  fit_label = "Fit C") {
  fixed_params <- with_gamma_fixed(fixed_params)
  tryCatch({
    cat(sprintf("  [%s] fixed params: %s\n",
                fit_label,
                if (length(fixed_params) > 0) paste(names(fixed_params), collapse = ", ") else "<none>"))
    fit_etas_bivariate(
      params_init = init_params, realiz = pp_all_bg,
      windowT = windowT_fit, windowS = win_km, m0 = ETAS_M0,
      control_state_space = control_ss, treated_state_space = treated_ss,
      background_rate_var = "W",
      treated_background_zero_before = 0,
      control_background_everywhere_before = 0,
      control_background_pre_mass_ratio = CTRL_BG_PRE_MASS_RATIO,
      beta_gr = BETA_GR,
      enforce_finite_trigger_moments = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
      p_lower_bound = ETAS_P_LOWER_BOUND,
      q_lower_bound = ETAS_Q_LOWER_BOUND,
      max_branching_radius = ETAS_BRANCHING_MAX,
      maxit = VANILLA_MAXIT, fixed_params = fixed_params, trace = 0,
      t_trunc = SEM_T_TRUNC_DAYS
    )
  }, error = function(e) {
    cat(sprintf("  [%s] bivariate+KDE fit error: %s\n", fit_label, e$message))
    NULL
  })
}

# ============================================================================
# 4D. Fit D: SEM bivariate ETAS with KDE background
# ============================================================================
cat("\n--- Fit D: SEM bivariate ETAS with KDE background ---\n")

biv_init_F <- apply_pre_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
fit_f <- function(init_params = biv_init_F,
                  fixed_params = NULL,
                  fit_label = "Fit D",
                  init_ctrl_params = A_ctrl,
                  init_treat_params = A_treat) {
  fixed_params <- with_gamma_fixed(fixed_params)
  tryCatch({
    init_params_use <- if (isTRUE(SEM_START_FROM_C)) {
      init_params
    } else {
      enforce_control_decay_start(init_params, init_ctrl_params)
    }
    cat(sprintf(
      "  [%s] start=%s  c=%.4f p=%.4f D=%.4f gamma=%.4f q=%.4f\n",
      fit_label,
      if (isTRUE(SEM_START_FROM_C)) "Fit C MLE" else "control-decay",
      as.numeric(init_params_use[["c"]]),
      as.numeric(init_params_use[["p"]]),
      as.numeric(init_params_use[["D"]]),
      as.numeric(init_params_use[["gamma"]]),
      as.numeric(init_params_use[["q"]])
    ))
    run_sem_fit(
      pp_data_in = pp_all_bg,
      partition_in = partition,
      partition_processes_in = partition_processes,
      state_spaces_in = state_spaces,
      init_params_in = init_params_use,
      init_ctrl_params_in = init_ctrl_params,
      init_treat_params_in = init_treat_params,
      fixed_params_in = fixed_params,
      background_rate_var_in = "W",
      control_background_pre_mass_ratio_in = CTRL_BG_PRE_MASS_RATIO,
      verbose_in = DF_VERBOSE,
      label = fit_label
    )
  }, error = function(e) {
    cat(sprintf("  [%s] SEM-biv+KDE error: %s\n", fit_label, e$message))
    NULL
  })
}

primary_kde_spec <- kde_variant_specs[[kde_primary_variant_id]]

sem_pilot_tuning <- NULL
if (!isTRUE(FIT_VARIABILITY_ONLY) && RUN_SEM_PILOT) {
  cat("\n--- SEM pilot tuning (Model D only) ---\n")
  pilot_cf <- parse_numeric_vector_env("OK_SEM_PILOT_CHANGE_FACTORS", c(0.005, 0.01, 0.02), lower = 1e-10)
  pilot_min <- parse_numeric_vector_env("OK_SEM_PILOT_MIN_MULTS", c(0.1, 0.2, 0.4), lower = 1e-10)
  pilot_max <- parse_numeric_vector_env("OK_SEM_PILOT_MAX_MULTS", c(1.5, 2.0, 3.0), lower = 1e-10)
  pilot_temp <- parse_numeric_vector_env("OK_SEM_PILOT_TEMPS", c(0.04, 0.08, 0.12), lower = 1e-10)

  pilot_grid <- expand.grid(
    change_factor = pilot_cf,
    min_mult = pilot_min,
    max_mult = pilot_max,
    temperature = pilot_temp,
    stringsAsFactors = FALSE
  )
  pilot_grid <- pilot_grid[pilot_grid$max_mult >= pilot_grid$min_mult, , drop = FALSE]

  pilot_cores_use <- if (!is.na(SEM_PILOT_CORES)) SEM_PILOT_CORES else N_CORES
  pilot_max_jobs <- min(nrow(pilot_grid), as.integer(pilot_cores_use))
  pilot_cap <- min(as.integer(SEM_PILOT_MAX_COMBOS), as.integer(pilot_max_jobs))
  if (nrow(pilot_grid) > pilot_cap) {
    set.seed(derive_run_seed(OK_GLOBAL_SEED, label = "ok_sem_pilot_grid"))
    keep_idx <- sort(sample(seq_len(nrow(pilot_grid)), pilot_cap))
    pilot_grid <- pilot_grid[keep_idx, , drop = FALSE]
  }
  pilot_jobs <- split(pilot_grid, seq_len(nrow(pilot_grid)))
  cat(sprintf("  pilot combos: %d (cores=%d, inner_iter=%d)\n",
              length(pilot_jobs), pilot_cores_use, SEM_PILOT_INNER_ITER))
  if (length(pilot_jobs) < pilot_cores_use) {
    cat(sprintf("  pilot note: grid has fewer combos than cores (%d < %d).\n",
                length(pilot_jobs), pilot_cores_use))
  }

  run_one_pilot <- function(job_df) {
    cf <- as.numeric(job_df$change_factor[[1]])
    mn <- as.numeric(job_df$min_mult[[1]])
    mx <- as.numeric(job_df$max_mult[[1]])
    tt <- as.numeric(job_df$temperature[[1]])
    label <- sprintf("pilot_D_cf%.4f_min%.3f_max%.3f_T%.3f", cf, mn, mx, tt)
    fit <- tryCatch(
      run_sem_fit(
        pp_data_in = pp_all_bg,
        partition_in = partition,
        partition_processes_in = partition_processes,
        state_spaces_in = state_spaces,
        init_params_in = biv_init_F,
        fixed_params_in = primary_kde_spec$fixed_params,
        background_rate_var_in = "W",
        control_background_pre_mass_ratio_in = CTRL_BG_PRE_MASS_RATIO,
        sem_inner_iter_in = SEM_PILOT_INNER_ITER,
        sem_change_factor_in = cf,
        sem_change_factor_min_mult_in = mn,
        sem_change_factor_max_mult_in = mx,
        sem_optim_method_in = "sample_weighted",
        sem_selection_temperature_in = tt,
        verbose_in = FALSE,
        label = label
      ),
      error = function(e) NULL
    )
    if (is.null(fit) || is.null(fit$adaptive)) {
      return(data.frame(
        change_factor = cf, min_mult = mn, max_mult = mx, temperature = tt,
        ok = FALSE, final_metric = NA_real_,
        mean_avg_flips_tail = NA_real_, mean_acc_flips_tail = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    ad <- fit$adaptive
    metric_vec <- as.numeric(ad$metrics)
    metric_vec <- metric_vec[is.finite(metric_vec)]
    tail_n <- min(20L, length(ad$average_flips), length(ad$max_metric_flips))
    mean_avg_tail <- if (tail_n > 0L) mean(tail(as.numeric(ad$average_flips), tail_n), na.rm = TRUE) else NA_real_
    mean_acc_tail <- if (tail_n > 0L) mean(tail(as.numeric(ad$max_metric_flips), tail_n), na.rm = TRUE) else NA_real_
    data.frame(
      change_factor = cf, min_mult = mn, max_mult = mx, temperature = tt,
      ok = TRUE,
      final_metric = if (length(metric_vec) > 0L) tail(metric_vec, 1L) else NA_real_,
      mean_avg_flips_tail = mean_avg_tail,
      mean_acc_flips_tail = mean_acc_tail,
      stringsAsFactors = FALSE
    )
  }

  pilot_out <- if (length(pilot_jobs) > 0L) {
    run_parallel(
      pilot_jobs, run_one_pilot,
      cores = pilot_cores_use,
      label = "sem-pilot-fitD"
    )
  } else {
    list()
  }
  pilot_df <- if (length(pilot_out) > 0L) do.call(rbind, pilot_out) else NULL
  if (!is.null(pilot_df) && nrow(pilot_df) > 0) {
    ok_df <- pilot_df[pilot_df$ok & is.finite(pilot_df$final_metric), , drop = FALSE]
    if (nrow(ok_df) > 0) {
      ok_df <- ok_df %>%
        mutate(
          rank_metric = rank(-.data$final_metric, ties.method = "min"),
          rank_acc = rank(.data$mean_acc_flips_tail, ties.method = "min"),
          rank_avg = rank(.data$mean_avg_flips_tail, ties.method = "min"),
          pilot_score = .data$rank_metric + .data$rank_acc + .data$rank_avg
        ) %>%
        arrange(.data$pilot_score, desc(.data$final_metric), .data$mean_acc_flips_tail)
      best <- ok_df[1, , drop = FALSE]
      SEM_CHANGE_FACTOR <- as.numeric(best$change_factor[[1]])
      SEM_CHANGE_FACTOR_MIN_MULT <- as.numeric(best$min_mult[[1]])
      SEM_CHANGE_FACTOR_MAX_MULT <- as.numeric(best$max_mult[[1]])
      SEM_OPTIM_METHOD <- "sample_weighted"
      SEM_SELECTION_TEMPERATURE <- as.numeric(best$temperature[[1]])
      cat(sprintf(
        "  pilot best: change_factor=%.4f min_mult=%.3f max_mult=%.3f temp=%.3f final_metric=%.3f avg_tail=%.2f acc_tail=%.2f\n",
        SEM_CHANGE_FACTOR, SEM_CHANGE_FACTOR_MIN_MULT, SEM_CHANGE_FACTOR_MAX_MULT, SEM_SELECTION_TEMPERATURE,
        as.numeric(best$final_metric[[1]]), as.numeric(best$mean_avg_flips_tail[[1]]), as.numeric(best$mean_acc_flips_tail[[1]])
      ))
      sem_pilot_tuning <- list(
        enabled = TRUE,
        inner_iter = SEM_PILOT_INNER_ITER,
        cores = pilot_cores_use,
        grid_size = nrow(pilot_grid),
        results = ok_df,
        best = best
      )
    } else {
      cat("  pilot warning: no successful pilot fits; using existing SEM defaults.\n")
      sem_pilot_tuning <- list(enabled = TRUE, inner_iter = SEM_PILOT_INNER_ITER, cores = pilot_cores_use, grid_size = nrow(pilot_grid), results = pilot_df, best = NULL)
    }
  }
}

if (!isTRUE(FIT_VARIABILITY_ONLY) && !isTRUE(BOOTSTRAP_ONLY) && !isTRUE(T_TRUNC_SENS_ONLY)) {
cat("\n--- Step 4 unified dispatch: running all primary-tessellation fits in parallel ---\n")
fit_jobs_all <- list()
if (!isTRUE(CD_ONLY) && !isTRUE(UNIV_KDE_ONLY)) {
  fit_jobs_all <- c(
    fit_jobs_all,
    list(
      list(kind = "A_hom_naive", variant_id = NA_character_),
      list(kind = "B_hom_sem", variant_id = NA_character_)
    )
  )
} else {
  cat(sprintf("  %s: skipping homogeneous Fits A/B.\n",
              if (isTRUE(UNIV_KDE_ONLY)) "UNIV_KDE_ONLY" else "CD_ONLY"))
}
if (!isTRUE(AB_ONLY) && !isTRUE(UNIV_KDE_ONLY)) {
  for (vid in kde_variant_ids) {
    fit_jobs_all[[length(fit_jobs_all) + 1L]] <- list(kind = "C_kde_naive", variant_id = vid)
    fit_jobs_all[[length(fit_jobs_all) + 1L]] <- list(kind = "D_kde_sem", variant_id = vid)
  }
} else {
  cat(sprintf("  %s: skipping KDE Fits C/D.\n",
              if (isTRUE(UNIV_KDE_ONLY)) "UNIV_KDE_ONLY" else "AB_ONLY"))
}
c_start_params <- list()
run_all_fit_job <- function(job) {
  t0 <- proc.time()[["elapsed"]]
  out <- NULL
  fit_label <- NA_character_
  kind <- job$kind
  vid <- job$variant_id
  cat(sprintf("    [fit-job:%s/%s] start pid=%d mem=%s\n",
              kind, ifelse(is.na(vid), "base", vid), Sys.getpid(), mem_snapshot()))
  if (kind == "A_hom_naive") {
    fit_label <- "Fit A"
    out <- fit_b()
  } else if (kind == "B_hom_sem") {
    fit_label <- "Fit B"
    out <- fit_d()
  } else if (kind == "C_kde_naive") {
    spec <- kde_variant_specs[[vid]]
    fit_label <- kde_fit_label("E", vid)
    out <- fit_e(
      init_params = biv_init_E,
      fixed_params = spec$fixed_params,
      fit_label = fit_label
    )
  } else if (kind == "D_kde_sem") {
    spec <- kde_variant_specs[[vid]]
    fit_label <- kde_fit_label("F", vid)
    init_d <- biv_init_F
    init_ctrl_d <- A_ctrl
    init_treat_d <- A_treat
    if (isTRUE(SEM_START_FROM_C)) {
      c_par <- c_start_params[[as.character(vid)]]
      if (!is.null(c_par)) {
        init_d <- c_par
        marg <- extract_marginals(c_par)
        if (!is.null(marg)) {
          init_ctrl_d <- marg$ctrl
          init_treat_d <- marg$treat
        }
      }
    }
    out <- fit_f(
      init_params = init_d,
      fixed_params = spec$fixed_params,
      fit_label = fit_label,
      init_ctrl_params = init_ctrl_d,
      init_treat_params = init_treat_d
    )
  }
  complete <- if (kind %in% c("A_hom_naive", "C_kde_naive")) {
    valid_biv_fit(out)
  } else {
    valid_biv_sem(out)
  }
  if (!isTRUE(complete)) out <- NULL
  elapsed <- proc.time()[["elapsed"]] - t0
  cat(sprintf("    [fit-job:%s/%s] done in %.1fs status=%s mem=%s\n",
              kind, ifelse(is.na(vid), "base", vid), elapsed,
              ifelse(is.null(out), "failed", "ok"), mem_snapshot()))
  list(
    kind = kind,
    variant_id = vid,
    fit_label = fit_label,
    obj = out,
    elapsed = elapsed
  )
}
if (length(fit_jobs_all) < 1L) {
  fit_all_out <- list()
} else if (isTRUE(SEM_START_FROM_C)) {
  c_jobs <- Filter(function(j) identical(j$kind, "C_kde_naive"), fit_jobs_all)
  d_jobs <- Filter(function(j) identical(j$kind, "D_kde_sem"), fit_jobs_all)
  other_jobs <- Filter(
    function(j) !j$kind %in% c("C_kde_naive", "D_kde_sem"),
    fit_jobs_all
  )
  cat("  OK_SEM_START_FROM_C: Fit C first, then start D from that C MLE.\n")
  run_job_batch <- function(jobs, label) {
    if (length(jobs) < 1L) return(list())
    if (N_CORES > 1L && length(jobs) > 1L) {
      run_parallel(jobs, run_all_fit_job, cores = min(length(jobs), N_CORES), label = label)
    } else {
      lapply(jobs, run_all_fit_job)
    }
  }
  other_out <- run_job_batch(other_jobs, "fit-other-primary")
  c_out <- run_job_batch(c_jobs, "fit-C-before-D")
  for (row in c_out) {
    if (!is.null(row$obj) && valid_biv_fit(row$obj)) {
      c_start_params[[as.character(row$variant_id)]] <- row$obj$par
    }
  }
  d_out <- run_job_batch(d_jobs, "fit-D-from-C")
  fit_all_out <- c(as.list(other_out), as.list(c_out), as.list(d_out))
} else if (N_CORES > 1L && length(fit_jobs_all) > 1L) {
  fit_all_out <- run_parallel(
    fit_jobs_all, run_all_fit_job,
    cores = min(length(fit_jobs_all), N_CORES),
    label = "fit-all-primary-jobs"
  )
} else {
  fit_all_out <- lapply(fit_jobs_all, run_all_fit_job)
}
fit_all_out <- as.list(fit_all_out)
get_fit_out <- function(kind, variant_id = NA_character_) {
  idx <- which(vapply(
    fit_all_out,
    function(z) identical(z$kind, kind) && identical(as.character(z$variant_id), as.character(variant_id)),
    logical(1)
  ))
  if (length(idx) < 1L) return(NULL)
  fit_all_out[[idx[1]]]
}

fit_A_row <- get_fit_out("A_hom_naive")
fit_B_row <- get_fit_out("B_hom_sem")
fitB <- if (!is.null(fit_A_row)) fit_A_row$obj else NULL
semD <- if (!is.null(fit_B_row)) fit_B_row$obj else NULL
fit_B_elapsed <- if (!is.null(fit_A_row)) fit_A_row$elapsed else NA_real_
fit_D_elapsed <- if (!is.null(fit_B_row)) fit_B_row$elapsed else NA_real_
add_timing_row(
  stage = "fit_B_naive_bivariate",
  elapsed_sec = fit_B_elapsed,
  status = if (!is.null(fitB)) "ok" else "failed",
  detail = "elapsed from unified fit dispatch"
)
add_timing_row(
  stage = "fit_D_sem_bivariate",
  elapsed_sec = fit_D_elapsed,
  status = if (!is.null(semD)) "ok" else "failed",
  detail = "elapsed from unified fit dispatch"
)

if ((isTRUE(CD_ONLY) || isTRUE(UNIV_KDE_ONLY)) && is.null(fitB) && is.null(semD)) {
  B_params <- NULL
  B_loglik <- NA_real_
  D_params <- NULL
  D_ctrl <- NULL
  D_treat <- NULL
  cat(sprintf("  %s: Fit A/B not run (params left NULL).\n",
              if (isTRUE(UNIV_KDE_ONLY)) "UNIV_KDE_ONLY" else "CD_ONLY"))
} else {
  B_params <- if (!is.null(fitB) && !is.null(fitB$par)) fitB$par else NULL
  B_loglik <- if (!is.null(fitB)) fitB$value else NA_real_
  if (!is.null(B_params)) {
    cat("  Fit A params:", paste(biv_names, round(B_params, 4), sep = "=", collapse = ", "), "\n")
  } else {
    cat("  Fit A failed; params=NULL and downstream ATE will be skipped.\n")
  }
  if (!is.null(semD) && valid_biv_sem(semD)) {
    D_params <- semD$etas_bivariate_params
    D_ctrl <- semD$hawkes_params_control
    D_treat <- semD$hawkes_params_treated
    if (!is.null(D_params)) {
      cat("  Fit B params:", paste(biv_names, round(D_params, 4), sep = "=", collapse = ", "), "\n")
    } else {
      cat("  Fit B returned no parameters; params=NULL and downstream ATE will be skipped.\n")
    }
  } else {
    D_params <- NULL
    D_ctrl <- NULL
    D_treat <- NULL
    cat("  Fit B failed; params=NULL and downstream ATE will be skipped.\n")
  }
}

fitE <- NULL
semF <- NULL
E_params <- NULL
E_loglik <- NA_real_
F_params <- NULL
F_ctrl <- NULL
F_treat <- NULL
for (vid in kde_variant_ids) {
  spec <- kde_variant_specs[[vid]]
  row_e <- get_fit_out("C_kde_naive", vid)
  row_f <- get_fit_out("D_kde_sem", vid)
  fitE_var <- if (!is.null(row_e)) row_e$obj else NULL
  semF_var <- if (!is.null(row_f)) row_f$obj else NULL
  E_var_params <- if (valid_biv_fit(fitE_var)) fitE_var$par else NULL
  E_var_loglik <- if (valid_biv_fit(fitE_var)) fitE_var$value else NA_real_
  F_var_params <- if (valid_biv_sem(semF_var)) {
    semF_var$etas_bivariate_params
  } else {
    NULL
  }
  F_var_ctrl <- if (valid_biv_sem(semF_var)) semF_var$hawkes_params_control else NULL
  F_var_treat <- if (valid_biv_sem(semF_var)) semF_var$hawkes_params_treated else NULL

  kde_variant_fits$E[[vid]] <- list(
    id = spec$id,
    label = spec$label,
    fixed_params = spec$fixed_params,
    fit = fitE_var,
    params = E_var_params,
    objective = E_var_loglik,
    status = if (!is.null(E_var_params)) "ok" else "failed"
  )
  kde_variant_fits$F[[vid]] <- list(
    id = spec$id,
    label = spec$label,
    fixed_params = spec$fixed_params,
    fit = semF_var,
    params = F_var_params,
    hawkes_params_control = F_var_ctrl,
    hawkes_params_treated = F_var_treat,
    status = if (!is.null(F_var_params)) "ok" else "failed"
  )

  if (identical(vid, kde_primary_variant_id)) {
    fitE <- fitE_var
    semF <- semF_var
    E_params <- E_var_params
    E_loglik <- E_var_loglik
    F_params <- F_var_params
    F_ctrl <- F_var_ctrl
    F_treat <- F_var_treat
    add_timing_row(
      stage = "fit_E_naive_bivariate_kde",
      elapsed_sec = if (!is.null(row_e)) row_e$elapsed else NA_real_,
      status = if (!is.null(fitE_var)) "ok" else "failed",
      detail = "elapsed from unified fit dispatch"
    )
    add_timing_row(
      stage = "fit_F_sem_bivariate_kde",
      elapsed_sec = if (!is.null(row_f)) row_f$elapsed else NA_real_,
      status = if (!is.null(semF_var)) "ok" else "failed",
      detail = "elapsed from unified fit dispatch"
    )
    if (!is.null(E_var_params)) {
      cat("  Fit C params:", paste(biv_names, round(E_var_params, 4), sep = "=", collapse = ", "), "\n")
    } else if (!isTRUE(UNIV_KDE_ONLY)) {
      cat("  Fit C failed; params=NULL and downstream ATE will be skipped.\n")
    }
    if (!is.null(F_var_params)) {
      cat("  Fit D params:", paste(biv_names, round(F_var_params, 4), sep = "=", collapse = ", "), "\n")
    } else if (!isTRUE(UNIV_KDE_ONLY)) {
      cat("  Fit D failed; params=NULL and downstream ATE will be skipped.\n")
    }
  } else {
    add_timing_row(
      stage = sprintf("fit_E_kde_variant_%s", spec$id),
      elapsed_sec = if (!is.null(row_e)) row_e$elapsed else NA_real_,
      status = if (!is.null(fitE_var)) "ok" else "failed",
      detail = "elapsed from unified fit dispatch"
    )
    add_timing_row(
      stage = sprintf("fit_F_kde_variant_%s", spec$id),
      elapsed_sec = if (!is.null(row_f)) row_f$elapsed else NA_real_,
      status = if (!is.null(semF_var)) "ok" else "failed",
      detail = "elapsed from unified fit dispatch"
    )
    cat(sprintf("  %s status: %s\n", kde_fit_label("E", spec$id),
                if (!is.null(E_var_params)) "ok" else "failed (params=NULL)"))
    cat(sprintf("  %s status: %s\n", kde_fit_label("F", spec$id),
                if (!is.null(F_var_params)) "ok" else "failed (params=NULL)"))
  }
}

report_pq_boundary <- function(label, params, tolerance = 1e-3) {
  if (is.null(params) || !all(c("p", "q") %in% names(params))) return(invisible(NULL))
  p_val <- suppressWarnings(as.numeric(params[["p"]]))
  q_val <- suppressWarnings(as.numeric(params[["q"]]))
  if (!all(is.finite(c(p_val, q_val)))) return(invisible(NULL))
  cat(sprintf(
    "  %s kernel margins: p-1=%.6g, q-1=%.6g\n",
    label, p_val - 1, q_val - 1
  ))
  near <- c(p = p_val - 1, q = q_val - 1) <= tolerance
  if (any(near)) {
    warning(sprintf(
      "%s is on the normalized-kernel boundary for %s; do not treat that decay estimate as interior.",
      label, paste(names(near)[near], collapse = ", ")
    ))
  }
  invisible(NULL)
}
report_pq_boundary("Fit A", B_params)
report_pq_boundary("Fit B", D_params)
report_pq_boundary("Fit C", E_params)
report_pq_boundary("Fit D", F_params)

G_params <- if (!is.null(kde_variant_fits$E$productivity_free) &&
                !is.null(kde_variant_fits$E$productivity_free$params)) {
  kde_variant_fits$E$productivity_free$params
} else {
  NULL
}
H_params <- if (!is.null(kde_variant_fits$F$productivity_free) &&
                !is.null(kde_variant_fits$F$productivity_free$params)) {
  kde_variant_fits$F$productivity_free$params
} else {
  NULL
}

# ============================================================================
# 4G–4J. Univariate ETAS comparison (no cross-excitation)
#   Public letters G/H: independent margins + KDE background (internal fitI/semJ)
#   Public letters I/J: independent margins, homogeneous background (internal fitK/semL)
# ============================================================================
cat("\n--- Step 4G-4J: Univariate ETAS dispatch (public G/H univ+KDE, I/J univ homog) ---\n")
if (isTRUE(CD_ONLY) || isTRUE(AB_ONLY)) {
  cat(sprintf("  %s: skipping univariate Fits G–J (I/J/K/L jobs).\n",
              if (isTRUE(AB_ONLY)) "AB_ONLY" else "CD_ONLY"))
  fitI <- NULL
  semJ <- NULL
  fitK <- NULL
  semL <- NULL
  I_params <- NULL
  J_params <- NULL
  K_params <- NULL
  L_params <- NULL
  I_ctrl <- NULL
  I_treat <- NULL
  J_ctrl <- NULL
  J_treat <- NULL
  K_ctrl <- NULL
  K_treat <- NULL
  L_ctrl <- NULL
  L_treat <- NULL
  pp_post_sem_J <- pp_post_bg
  pp_post_sem_L <- pp_post
} else {

fit_indep_pair <- function(pp_data_in,
                           background_rate_var_in = NULL,
                           fixed_params_in = NULL,
                           ctrl_init_in = A_ctrl,
                           treat_init_in = A_treat,
                           fit_label = "Fit I") {
  fixed_params_in <- with_gamma_fixed(fixed_params_in)
  all_dat <- as.data.frame(pp_data_in)
  pre_history <- all_dat[
    all_dat$t < windowT_post[1],
    c("x", "y", "t", "mag"),
    drop = FALSE
  ]
  treated_history <- pre_history[0, , drop = FALSE]
  dat <- all_dat[all_dat$t >= windowT_post[1], , drop = FALSE]
  if (nrow(dat) < 1L) return(NULL)
  ctrl_dat <- dat[dat$location_process == "control", , drop = FALSE]
  treat_dat <- dat[dat$location_process == "treated", , drop = FALSE]
  if (nrow(ctrl_dat) < 1L || nrow(treat_dat) < 1L) return(NULL)
  # Multistart over interior productivity values; seed from caller init + VANILLA_STARTS.
  start_pool <- c(list(ctrl_init_in), VANILLA_STARTS)
  starts_use <- lapply(start_pool, apply_pre_init_etas)
  valid_margin_fit <- function(fit) {
    if (is.null(fit) || is.null(fit$par) || is.null(names(fit$par)) ||
        !all(etas_names %in% names(fit$par)) ||
        !is.finite(fit$value) ||
        is.null(fit$convergence) || as.integer(fit$convergence)[1L] != 0L) {
      return(FALSE)
    }
    par <- suppressWarnings(as.numeric(fit$par[etas_names]))
    names(par) <- etas_names
    gap <- BETA_GR - par[["alpha_m"]]
    eta <- if (is.finite(gap) && gap > 0) {
      par[["A"]] * BETA_GR / gap
    } else {
      Inf
    }
    all(is.finite(par)) && is.finite(eta) && eta < ETAS_BRANCHING_MAX
  }
  fit_one_margin <- function(realiz_m, zbr_m, history_m, starts_m, side_label) {
    best_fit <- NULL
    best_val <- -Inf
    for (s in starts_m) {
      fit <- tryCatch(
        fit_etas(
          params_init = s, realiz = realiz_m,
          windowT = windowT_post, windowS = win_km, m0 = ETAS_M0,
          maxit = VANILLA_MAXIT, fixed_params = fixed_params_in,
          zero_background_region = zbr_m,
          background_rate_var = background_rate_var_in,
          beta_gr = BETA_GR,
          enforce_finite_trigger_moments = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
          p_lower_bound = ETAS_P_LOWER_BOUND,
          q_lower_bound = ETAS_Q_LOWER_BOUND,
          max_branching_ratio = ETAS_BRANCHING_MAX,
          t_trunc = SEM_T_TRUNC_DAYS,
          history = history_m,
          log_transform = TRUE,
          soft_branching_barrier = TRUE,
          polish_productivity = TRUE,
          interior_restart = TRUE
        ),
        error = function(e) {
          cat(sprintf(
            "  [%s] %s independent fit error: %s\n",
            fit_label, side_label, e$message
          ))
          NULL
        }
      )
      if (valid_margin_fit(fit) && fit$value > best_val) {
        best_fit <- fit
        best_val <- fit$value
      }
    }
    best_fit
  }
  treat_start_pool <- c(list(treat_init_in), VANILLA_STARTS)
  treat_starts_use <- lapply(treat_start_pool, apply_pre_init_etas)
  fit_ctrl <- fit_one_margin(
    ctrl_dat, treated_ss, pre_history, starts_use, "control"
  )
  fit_treat <- fit_one_margin(
    treat_dat, control_ss, treated_history, treat_starts_use, "treated"
  )
  if (!valid_margin_fit(fit_ctrl) || !valid_margin_fit(fit_treat)) {
    return(NULL)
  }
  ctrl_par <- as.list(fit_ctrl$par)
  treat_par <- as.list(fit_treat$par)
  names(ctrl_par) <- etas_names
  names(treat_par) <- etas_names
  list(
    control = ctrl_par,
    treated = treat_par,
    fit_control = fit_ctrl,
    fit_treated = fit_treat,
    objective = sum(c(
      if (!is.null(fit_ctrl)) as.numeric(fit_ctrl$value) else NA_real_,
      if (!is.null(fit_treat)) as.numeric(fit_treat$value) else NA_real_
    ), na.rm = TRUE)
  )
}

fit_i <- function() {
  tryCatch(
    fit_indep_pair(
      pp_data_in = pp_all_bg,
      background_rate_var_in = "W",
      # Structural free except gamma=0; finite-moment constraints via likelihood.
      fixed_params_in = with_gamma_fixed(NULL),
      ctrl_init_in = A_ctrl,
      treat_init_in = A_treat,
      fit_label = "Fit I"
    ),
    error = function(e) {
      cat("  Fit I (naive univariate+KDE) error:", e$message, "\n")
      NULL
    }
  )
}

fit_j <- function() {
  tryCatch(
    run_sem_fit(
      pp_data_in = pp_all_bg,
      partition_in = partition,
      partition_processes_in = partition_processes,
      state_spaces_in = state_spaces,
      init_params_in = NULL,
      init_ctrl_params_in = A_ctrl,
      init_treat_params_in = A_treat,
      model_type_in = "etas",
      fixed_params_in = with_gamma_fixed(NULL),
      background_rate_var_in = "W",
      sem_inner_iter_in = SEM_INNER_ITER,
      verbose_in = DF_VERBOSE,
      label = "Fit J"
    ),
    error = function(e) {
      cat("  Fit J (SEM univariate+KDE) error:", e$message, "\n")
      NULL
    }
  )
}

fit_k <- function() {
  tryCatch(
    fit_indep_pair(
      pp_data_in = pp_all,
      background_rate_var_in = NULL,
      fixed_params_in = with_gamma_fixed(NULL),
      ctrl_init_in = A_ctrl,
      treat_init_in = A_treat,
      fit_label = "Fit K"
    ),
    error = function(e) {
      cat("  Fit K (naive univariate, no KDE) error:", e$message, "\n")
      NULL
    }
  )
}

fit_l <- function() {
  tryCatch(
    run_sem_fit(
      pp_data_in = pp_all,
      partition_in = partition,
      partition_processes_in = partition_processes,
      state_spaces_in = state_spaces,
      init_params_in = NULL,
      init_ctrl_params_in = A_ctrl,
      init_treat_params_in = A_treat,
      model_type_in = "etas",
      fixed_params_in = with_gamma_fixed(NULL),
      background_rate_var_in = NULL,
      sem_inner_iter_in = SEM_INNER_ITER,
      verbose_in = DF_VERBOSE,
      label = "Fit L"
    ),
    error = function(e) {
      cat("  Fit L (SEM univariate, no KDE) error:", e$message, "\n")
      NULL
    }
  )
}

univ_jobs <- if (isTRUE(UNIV_KDE_ONLY)) c("I", "J") else c("I", "J", "K", "L")
if (isTRUE(UNIV_KDE_ONLY)) {
  cat("  UNIV_KDE_ONLY: running G/H (univ+KDE); skipping homogeneous I/J.\n")
}
run_one_univ_job <- function(tag) {
  t0 <- proc.time()[["elapsed"]]
  out_obj <- NULL
  cat(sprintf("    [univ-job:%s] start pid=%d mem=%s\n", tag, Sys.getpid(), mem_snapshot()))
  if (tag == "I") out_obj <- fit_i()
  if (tag == "J") out_obj <- fit_j()
  if (tag == "K") out_obj <- fit_k()
  if (tag == "L") out_obj <- fit_l()
  complete <- if (tag %in% c("I", "K")) {
    !is.null(out_obj) && !is.null(out_obj$control) && !is.null(out_obj$treated)
  } else {
    !is.null(out_obj) && !is.null(out_obj$hawkes_params_control) &&
      !is.null(out_obj$hawkes_params_treated)
  }
  if (!isTRUE(complete)) out_obj <- NULL
  elapsed <- proc.time()[["elapsed"]] - t0
  cat(sprintf("    [univ-job:%s] done in %.1fs status=%s mem=%s\n",
              tag, elapsed, ifelse(is.null(out_obj), "failed", "ok"), mem_snapshot()))
  list(tag = tag, obj = out_obj, elapsed = elapsed)
}
if (N_CORES > 1L && length(univ_jobs) > 1L) {
  univ_out <- run_parallel(
    univ_jobs, run_one_univ_job,
    cores = min(length(univ_jobs), N_CORES),
    label = "univariate-jobs"
  )
} else {
  univ_out <- lapply(univ_jobs, run_one_univ_job)
}
univ_out <- as.list(univ_out)
get_univ_job <- function(tag) {
  idx <- which(vapply(univ_out, function(z) identical(z$tag, tag), logical(1)))
  if (length(idx) < 1L) return(NULL)
  univ_out[[idx[1]]]
}
row_i <- get_univ_job("I")
row_j <- get_univ_job("J")
row_k <- get_univ_job("K")
row_l <- get_univ_job("L")
fitI <- if (!is.null(row_i)) row_i$obj else NULL
semJ <- if (!is.null(row_j)) row_j$obj else NULL
fitK <- if (!is.null(row_k)) row_k$obj else NULL
semL <- if (!is.null(row_l)) row_l$obj else NULL
fit_I_elapsed <- if (!is.null(row_i)) row_i$elapsed else NA_real_
fit_J_elapsed <- if (!is.null(row_j)) row_j$elapsed else NA_real_
fit_K_elapsed <- if (!is.null(row_k)) row_k$elapsed else NA_real_
fit_L_elapsed <- if (!is.null(row_l)) row_l$elapsed else NA_real_
add_timing_row(
  stage = "fit_I_naive_univariate_kde",
  elapsed_sec = fit_I_elapsed,
  status = if (!is.null(fitI)) "ok" else "failed",
  detail = "elapsed from univariate dispatch"
)
add_timing_row(
  stage = "fit_J_sem_univariate_kde",
  elapsed_sec = fit_J_elapsed,
  status = if (!is.null(semJ)) "ok" else "failed",
  detail = "elapsed from univariate dispatch"
)
add_timing_row(
  stage = "fit_K_naive_univariate_homog",
  elapsed_sec = fit_K_elapsed,
  status = if (!is.null(fitK)) "ok" else "failed",
  detail = "elapsed from univariate dispatch"
)
add_timing_row(
  stage = "fit_L_sem_univariate_homog",
  elapsed_sec = fit_L_elapsed,
  status = if (!is.null(semL)) "ok" else "failed",
  detail = "elapsed from univariate dispatch"
)
I_ctrl <- if (!is.null(fitI)) fitI$control else NULL
I_treat <- if (!is.null(fitI)) fitI$treated else NULL
J_ctrl <- if (!is.null(semJ)) semJ$hawkes_params_control else NULL
J_treat <- if (!is.null(semJ)) semJ$hawkes_params_treated else NULL
K_ctrl <- if (!is.null(fitK)) fitK$control else NULL
K_treat <- if (!is.null(fitK)) fitK$treated else NULL
L_ctrl <- if (!is.null(semL)) semL$hawkes_params_control else NULL
L_treat <- if (!is.null(semL)) semL$hawkes_params_treated else NULL
I_params <- if (!is.null(I_ctrl) && !is.null(I_treat)) list(control = I_ctrl, treated = I_treat) else NULL
J_params <- if (!is.null(J_ctrl) && !is.null(J_treat)) list(control = J_ctrl, treated = J_treat) else NULL
K_params <- if (!is.null(K_ctrl) && !is.null(K_treat)) list(control = K_ctrl, treated = K_treat) else NULL
L_params <- if (!is.null(L_ctrl) && !is.null(L_treat)) list(control = L_ctrl, treated = L_treat) else NULL
pp_post_sem_J <- if (!is.null(semJ) && !is.null(semJ$adaptive$adaptive_labelling)) {
  semJ$adaptive$adaptive_labelling
} else {
  pp_post_bg
}
pp_post_sem_J <- pp_post_sem_J[pp_post_sem_J$t >= 0, , drop = FALSE]
pp_post_sem_L <- if (!is.null(semL) && !is.null(semL$adaptive$adaptive_labelling)) {
  semL$adaptive$adaptive_labelling
} else {
  pp_post
}
pp_post_sem_L <- pp_post_sem_L[pp_post_sem_L$t >= 0, , drop = FALSE]

} # end !CD_ONLY univariate block

# ============================================================================
# 4K. KDE bandwidth sensitivity (primary partition, inhomogeneous C/D = all_free)
# ============================================================================
cat("\n--- Step 4K: KDE bandwidth sensitivity (primary partition, inhom C/D all_free) ---\n")
if (length(KDE_BW_SENS_KM) > 0L) {
  kde_bandwidth_specs <- lapply(KDE_BW_SENS_KM, function(km) {
    list(
      label = sprintf("%gkm", km),
      multiplier = NA_real_,
      sigma_km = as.numeric(km)
    )
  })
  cat(sprintf("  Absolute-km bandwidth grid: %s (SEM inner=%d, t_trunc=%.4g days)\n",
              paste(vapply(kde_bandwidth_specs, `[[`, character(1), "label"), collapse = ", "),
              as.integer(SENS_SEM_INNER_ITER), as.numeric(SEM_T_TRUNC_DAYS)))
} else if (RUN_SENSITIVITY) {
  kde_bandwidth_specs <- list(
    list(label = "diggle", multiplier = 1),
    list(label = "digglex2", multiplier = 2),
    list(label = "digglex5", multiplier = 5),
    list(label = "digglex10", multiplier = 10)
  )
} else {
  kde_bandwidth_specs <- list()
  cat("  Sensitivity disabled: skipping KDE bandwidth and partition sensitivity fits.\n")
}

run_kde_bandwidth_fit <- function(spec) {
  bw_label <- spec$label
  bw_mult <- spec$multiplier
  sigma_local <- if (!is.null(spec$sigma_km) && is.finite(spec$sigma_km) && spec$sigma_km > 0) {
    as.numeric(spec$sigma_km)
  } else {
    bw_mult * bw_diggle
  }
  cat(sprintf("  [BW %s] start (sigma=%.3f km, SEM inner=%d)\n",
              bw_label, sigma_local, as.integer(SENS_SEM_INNER_ITER)))

  lambda_local <- density(X_bg, sigma = sigma_local, edge = TRUE, at = "pixels")
  min_nz_local <- min(lambda_local$v[lambda_local$v > 0], na.rm = TRUE)
  lambda_local$v[lambda_local$v <= 0] <- min_nz_local

  bg_data <- build_background_weighted_data(lambda_local, control_ss, treated_ss)
  pp_post_bg_local <- bg_data$pp_post_bg
  pp_all_bg_local <- bg_data$pp_all_bg
  mass_ratio_local <- kde_pre_mass_ratio(lambda_local, control_ss)

  biv_init_local <- apply_pre_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
  fitE_local <- tryCatch({
    fit_etas_bivariate(
      params_init = biv_init_local, realiz = pp_all_bg_local,
      windowT = windowT_fit, windowS = win_km, m0 = ETAS_M0,
      control_state_space = control_ss, treated_state_space = treated_ss,
      background_rate_var = "W",
      treated_background_zero_before = 0,
      control_background_everywhere_before = 0,
      control_background_pre_mass_ratio = mass_ratio_local,
      beta_gr = BETA_GR,
      enforce_finite_trigger_moments = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
      p_lower_bound = ETAS_P_LOWER_BOUND,
      q_lower_bound = ETAS_Q_LOWER_BOUND,
      max_branching_radius = ETAS_BRANCHING_MAX,
      maxit = VANILLA_MAXIT, fixed_params = SENSITIVITY_FIXED_PARAMS, trace = 0,
      t_trunc = SEM_T_TRUNC_DAYS
    )
  }, error = function(e) {
    cat(sprintf("  [BW %s] Fit C error: %s\n", bw_label, e$message))
    NULL
  })
  E_params_local <- if (valid_biv_fit(fitE_local)) fitE_local$par else NULL
  E_loglik_local <- if (valid_biv_fit(fitE_local)) fitE_local$value else NA_real_
  d_start_local <- sem_init_from_c(E_params_local, biv_init_local)

  semF_local <- tryCatch({
    run_sem_fit(
      pp_data_in = pp_all_bg_local,
      partition_in = partition,
      partition_processes_in = partition_processes,
      state_spaces_in = state_spaces,
      init_params_in = d_start_local$init,
      init_ctrl_params_in = d_start_local$ctrl,
      init_treat_params_in = d_start_local$treat,
      fixed_params_in = SENSITIVITY_FIXED_PARAMS,
      background_rate_var_in = "W",
      control_background_pre_mass_ratio_in = mass_ratio_local,
      sem_inner_iter_in = SENS_SEM_INNER_ITER,
      verbose_in = FALSE,
      label = sprintf("BW %s Fit D", bw_label),
      sem_rng_label_in = OK_BW_SEM_RNG_LABEL,
      biv_n_threads_in = SENS_BIV_N_THREADS
    )
  }, error = function(e) {
    cat(sprintf("  [BW %s] Fit D error: %s\n", bw_label, e$message))
    NULL
  })
  F_params_local <- if (valid_biv_sem(semF_local)) {
    semF_local$etas_bivariate_params
  } else {
    NULL
  }

  pp_post_sem_local <- if (!is.null(semF_local) && !is.null(semF_local$adaptive$adaptive_labelling)) {
    semF_local$adaptive$adaptive_labelling
  } else {
    pp_post_bg_local
  }
  pp_post_sem_local <- pp_post_sem_local[pp_post_sem_local$t >= 0, ]
  n_relabel_local <- if (!is.null(semF_local) && !is.null(semF_local$adaptive$adaptive_labelling)) {
    lp <- semF_local$adaptive$adaptive_labelling
    lp <- lp[lp$t >= 0, , drop = FALSE]
    sum(lp$location_process != lp$inferred_process, na.rm = TRUE)
  } else { 0L }

  cat(sprintf("  [BW %s] done: relabelled=%d, loglik(E)=%.3f\n",
              bw_label, as.integer(n_relabel_local), as.numeric(E_loglik_local)))
  list(
    label = bw_label,
    multiplier = bw_mult,
    sigma = as.numeric(sigma_local),
    E_params = E_params_local,
    F_params = F_params_local,
    E_loglik = E_loglik_local,
    pp_post_bg = pp_post_bg_local,
    pp_post_sem = pp_post_sem_local,
    n_relabel = as.integer(n_relabel_local)
  )
}

# Executed later in a joint sensitivity dispatch with partition sensitivity.
kde_bandwidth_fits <- NULL

# ============================================================================
# 5. Alternative partitioning schemes
# ============================================================================
cat("\n--- Step 5: Alternative partitioning schemes ---\n")

D_scale <- trigger_range_km
sens_ids <- oklahoma_sensitivity_partition_labels(PRIMARY_PARTITION)
if (QUICK_CHECK) {
  want_grid <- intersect(sens_ids, "grid_coarse")
  grid_multipliers <- if (length(want_grid)) 1.0 else numeric(0)
  grid_diameters <- if (length(want_grid)) NA_real_ else numeric(0)
  grid_max_tiles <- if (length(want_grid)) 2000L else integer(0)
  if (length(want_grid)) {
    cat("  QUICK_CHECK mode: using one coarse grid (~100 tiles) as sensitivity\n")
  } else {
    cat("  QUICK_CHECK mode: primary is the coarse grid; no extra grid sensitivities\n")
  }
} else {
  want_grid <- grep("^grid_", sens_ids, value = TRUE)
  grid_multipliers <- if (length(want_grid)) {
    as.numeric(sub("^grid_([0-9.]+)R$", "\\1", want_grid))
  } else {
    numeric(0)
  }
  grid_diameters <- if (length(grid_multipliers)) grid_multipliers * trigger_range_km else numeric(0)
  grid_max_lookup <- if (TEST_MODE) {
    c(`1` = 2000L, `2` = 1600L, `5` = 1200L)
  } else {
    c(`1` = 5000L, `2` = 3000L, `5` = 1500L)
  }
  grid_max_tiles <- unname(grid_max_lookup[as.character(as.integer(grid_multipliers))])
  grid_max_tiles[is.na(grid_max_tiles)] <- if (TEST_MODE) 1600L else 3000L
  if (length(grid_diameters)) {
    cat(sprintf("  Triggering range = %.2f km; sensitivity grid diameters: %s km\n",
                trigger_range_km,
                paste(round(grid_diameters, 2), collapse = ", ")))
  } else {
    cat("  No grid partition sensitivities scheduled (primary already covers them).\n")
  }
}

grid_partitions <- if (length(grid_diameters)) {
  lapply(seq_along(grid_diameters), function(i) {
    build_grid_partition(
      diameter = grid_diameters[i],
      win = win_km,
      aoi_owin = aoi_owin,
      label = if (QUICK_CHECK) "grid_coarse" else sprintf("grid_%.1fR", grid_multipliers[i]),
      max_tiles = grid_max_tiles[i]
    )
  })
} else {
  list()
}

cat("\n  Building AOI partition (AOI = treated, rest = control)...\n")
aoi_partition <- tess(tiles = list(control = aoi_ctrl_ss, treated = aoi_treat_ss),
                      window = win_km)
aoi_procs <- c(control = "control", treated = "treated")
aoi_state_spaces <- list(control = aoi_ctrl_ss, treated = aoi_treat_ss)
aoi_treated_idx <- c(FALSE, TRUE)

aoi_scheme <- list(partition = aoi_partition, processes = aoi_procs,
                   state_spaces = aoi_state_spaces, label = "aoi_region",
                   treated_idx = aoi_treated_idx)

candidate_partitions <- c(
  list(county = county_scheme),
  setNames(grid_partitions, vapply(grid_partitions, `[[`, character(1), "label")),
  list(aoi_region = aoi_scheme)
)
all_partitions <- candidate_partitions[intersect(sens_ids, names(candidate_partitions))]

if (!RUN_SENSITIVITY || !RUN_PARTITION_SENSITIVITY) {
  all_partitions <- list()
  if (RUN_SENSITIVITY && !RUN_PARTITION_SENSITIVITY) {
    cat("  Partition sensitivity disabled: skipping county/grid/AOI refits.\n")
  }
}
cat(sprintf("  Total partition schemes scheduled: %d\n", length(all_partitions)))

analysis_plots$partition_map_ids <- vapply(all_partitions, `[[`, character(1), "label")
cat("  Partition map files are disabled; render partition map panels from results payload.\n")

# ============================================================================
# 5b. Joint sensitivity dispatch: bandwidth + partition fits (single parallel layer)
# ============================================================================
cat("\n--- Step 5b: Joint sensitivity dispatch (bandwidth + partition) ---\n")

assign_to_partition <- function(df, part_info) {
  df$location_process <- oklahoma_assign_partition_process(
    df$x, df$y, part_info$partition, part_info$processes
  )
  df$process <- df$location_process
  df$inferred_process <- df$location_process
  df$W <- 1.0
  df$background <- TRUE
  out <- df[!is.na(df$location_process), ]
  oklahoma_assert_label_support(
    out,
    part_info$state_spaces$control,
    part_info$state_spaces$treated,
    context = paste0("partition:", part_info$label)
  )
  out
}

run_biv_for_partition <- function(part_info) {
  label <- part_info$label
  cat(sprintf("\n  === Partition: %s ===\n", label))

  p_ctrl_ss <- part_info$state_spaces$control
  p_treat_ss <- part_info$state_spaces$treated
  p_state_spaces <- part_info$state_spaces
  p_treated_idx <- part_info$treated_idx

  pp_post_p <- assign_to_partition(pp_post[, c("x", "y", "t", "mag")], part_info)
  pp_pre_p  <- assign_to_partition(pp_pre[, c("x", "y", "t", "mag")], part_info)
  pp_post_p_bg <- rbind(
    normalize_bg_weights(pp_post_p[pp_post_p$location_process == "control", ],
                         p_ctrl_ss, lambda_im)$new_df,
    normalize_bg_weights(pp_post_p[pp_post_p$location_process == "treated", ],
                         p_treat_ss, lambda_im)$new_df
  )
  pp_post_p_bg <- pp_post_p_bg[order(pp_post_p_bg$t), ]
  # Pre-treatment everywhere is control: control-normalized weights for all
  # pre events (partition-specific control region).
  pp_pre_p_bg <- normalize_bg_weights(pp_pre_p, p_ctrl_ss, lambda_im)$new_df
  pp_pre_p_bg <- pp_pre_p_bg[order(pp_pre_p_bg$t), ]
  mass_ratio_p <- kde_pre_mass_ratio(lambda_im, p_ctrl_ss)
  # Enforce carry-over convention in partition-specific SEM runs too.
  pp_pre_p_bg$process <- "control"
  pp_pre_p_bg$inferred_process <- "control"
  pp_all_p  <- rbind(pp_pre_p_bg, pp_post_p_bg)
  pp_all_p  <- pp_all_p[order(pp_all_p$t), ]
  windowT_p <- c(min(pp_pre_p_bg$t, na.rm = TRUE), post_end_days)

  cat(sprintf("    Post events: %d ctrl, %d treat\n",
              sum(pp_post_p$location_process == "control"),
              sum(pp_post_p$location_process == "treated")))

  # Naive bivariate + KDE background
  biv_init_p <- apply_pre_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
  fitE_p <- tryCatch({
    fit_etas_bivariate(
      params_init = biv_init_p, realiz = pp_all_p,
      windowT = windowT_p, windowS = win_km, m0 = ETAS_M0,
      control_state_space = p_ctrl_ss, treated_state_space = p_treat_ss,
      background_rate_var = "W",
      treated_background_zero_before = 0,
      control_background_everywhere_before = 0,
      control_background_pre_mass_ratio = mass_ratio_p,
      beta_gr = BETA_GR,
      enforce_finite_trigger_moments = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
      p_lower_bound = ETAS_P_LOWER_BOUND,
      q_lower_bound = ETAS_Q_LOWER_BOUND,
      max_branching_radius = ETAS_BRANCHING_MAX,
      maxit = VANILLA_MAXIT, fixed_params = SENSITIVITY_FIXED_PARAMS, trace = 0,
      t_trunc = SEM_T_TRUNC_DAYS
    )
  }, error = function(e) { cat(sprintf("    [%s] Inhom naive fit error: %s\n", label, e$message)); NULL })

  E_params_p <- if (valid_biv_fit(fitE_p)) fitE_p$par else NULL
  d_start_p <- sem_init_from_c(E_params_p, biv_init_p)

  # SEM bivariate + KDE background
  semF_p <- tryCatch({
    run_sem_fit(
      pp_data_in = pp_all_p,
      partition_in = part_info$partition,
      partition_processes_in = part_info$processes,
      state_spaces_in = p_state_spaces,
      init_params_in = d_start_p$init,
      init_ctrl_params_in = d_start_p$ctrl,
      init_treat_params_in = d_start_p$treat,
      fixed_params_in = SENSITIVITY_FIXED_PARAMS,
      background_rate_var_in = "W",
      control_background_pre_mass_ratio_in = mass_ratio_p,
      sem_inner_iter_in = SENS_SEM_INNER_ITER,
      verbose_in = FALSE,
      label = sprintf("Partition %s Fit D", label),
      biv_n_threads_in = SENS_BIV_N_THREADS
    )
  }, error = function(e) { cat(sprintf("    [%s] Inhom SEM error: %s\n", label, e$message)); NULL })

  F_params_p <- if (valid_biv_sem(semF_p)) {
    semF_p$etas_bivariate_params
  } else {
    NULL
  }
  pp_post_sem_p <- if (!is.null(semF_p) && !is.null(semF_p$adaptive$adaptive_labelling)) {
    semF_p$adaptive$adaptive_labelling
  } else {
    pp_post_p_bg
  }
  pp_post_sem_p <- pp_post_sem_p[pp_post_sem_p$t >= 0, ]
  n_relabel_p <- if (!is.null(semF_p) && !is.null(semF_p$adaptive$adaptive_labelling)) {
    lp <- semF_p$adaptive$adaptive_labelling
    lp <- lp[lp$t >= 0, , drop = FALSE]
    sum(lp$location_process != lp$inferred_process, na.rm = TRUE)
  } else { 0L }

  cat(sprintf("    [%s] C status: %s\n", label,
              if (!is.null(E_params_p)) "ok" else "failed (params=NULL)"))
  cat(sprintf("    [%s] D status: %s\n", label,
              if (!is.null(F_params_p)) "ok" else "failed (params=NULL)"))

  list(label = label,
       fitE = if (TRIM_SENS_OBJECTS) NULL else fitE_p,
       E_params = E_params_p,
       semF = if (TRIM_SENS_OBJECTS) NULL else semF_p,
       F_params = F_params_p,
       pp_post_sem = pp_post_sem_p,
       n_relabel = as.integer(n_relabel_p),
       n_tiles = part_info$partition$n,
       n_treated = sum(p_treated_idx),
       n_post_ctrl = sum(pp_post_p$location_process == "control"),
       n_post_treat = sum(pp_post_p$location_process == "treated"))
}

sensitivity_jobs <- c(
  lapply(kde_bandwidth_specs, function(spec) list(type = "bandwidth", id = spec$label, payload = spec)),
  lapply(all_partitions, function(part_info) list(type = "partition", id = part_info$label, payload = part_info))
)
kde_bandwidth_labels <- vapply(kde_bandwidth_specs, `[[`, character(1), "label")
sens_n_parallel <- max(1L, min(SENS_CORES, max(1L, length(sensitivity_jobs))))
SENS_BIV_N_THREADS <- max(1L, as.integer(floor(N_CORES / sens_n_parallel)))
cat(sprintf("  Sensitivity jobs: %d bandwidth + %d partition = %d total\n",
            length(kde_bandwidth_specs), length(all_partitions), length(sensitivity_jobs)))
cat(sprintf("  Sensitivity parallel: %d workers, %d SEM threads/job (N_CORES=%d, SENS_CORES=%d)\n",
            sens_n_parallel, SENS_BIV_N_THREADS, N_CORES, SENS_CORES))

run_sensitivity_job <- function(job) {
  t0 <- proc.time()[["elapsed"]]
  cat(sprintf("    [sensitivity:%s/%s] start pid=%d mem=%s\n",
              as.character(job$type), as.character(job$id), Sys.getpid(), mem_snapshot()))
  if (identical(job$type, "bandwidth")) {
    out <- run_kde_bandwidth_fit(job$payload)
  } else if (identical(job$type, "partition")) {
    out <- run_biv_for_partition(job$payload)
  } else {
    out <- NULL
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  cat(sprintf("    [sensitivity:%s/%s] done in %.1fs status=%s mem=%s\n",
              as.character(job$type), as.character(job$id), elapsed,
              ifelse(is.null(out), "failed", "ok"), mem_snapshot()))
  list(type = job$type, id = job$id, out = out)
}

run_sensitivity_dispatch <- function() {
  if (!TEST_MODE && !QUICK_CHECK && N_CORES > 1 && length(sensitivity_jobs) > 1) {
    t_sensitivity_dispatch <- proc.time()[["elapsed"]]
    sens_out <- run_parallel(
      sensitivity_jobs, run_sensitivity_job,
      cores = min(SENS_CORES, length(sensitivity_jobs)),
      label = "sensitivity"
    )
    sensitivity_dispatch_elapsed <- proc.time()[["elapsed"]] - t_sensitivity_dispatch
  } else {
    t_sensitivity_dispatch <- proc.time()[["elapsed"]]
    if ((TEST_MODE || QUICK_CHECK) && length(sensitivity_jobs) > 1) {
      cat("  TEST/QUICK mode: running sensitivity jobs sequentially for stability.\n")
    }
    sens_out <- lapply(sensitivity_jobs, run_sensitivity_job)
    sensitivity_dispatch_elapsed <- proc.time()[["elapsed"]] - t_sensitivity_dispatch
  }
  add_timing_row(
    stage = "sensitivity_dispatch_total",
    elapsed_sec = sensitivity_dispatch_elapsed,
    status = "ok",
    detail = sprintf("jobs=%d", length(sensitivity_jobs))
  )

  kde_bandwidth_fits_local <- lapply(kde_bandwidth_labels, function(lbl) {
    idx <- which(vapply(sens_out, function(z) identical(z$type, "bandwidth") && identical(z$id, lbl), logical(1)))
    if (length(idx) == 0) return(NULL)
    sens_out[[idx[1]]]$out
  })
  names(kde_bandwidth_fits_local) <- kde_bandwidth_labels

  partition_results_local <- lapply(sapply(all_partitions, `[[`, "label"), function(lbl) {
    idx <- which(vapply(sens_out, function(z) identical(z$type, "partition") && identical(z$id, lbl), logical(1)))
    if (length(idx) == 0) return(NULL)
    sens_out[[idx[1]]]$out
  })
  names(partition_results_local) <- sapply(all_partitions, `[[`, "label")
  rm(sens_out)
  invisible(gc(verbose = FALSE))

  list(
    kde_bandwidth_fits = kde_bandwidth_fits_local,
    partition_results = partition_results_local
  )
}
# Sensitivity execution is intentionally deferred until after main-fit SEM plots
# and ATE payloads are generated.
kde_bandwidth_fits <- NULL
partition_results <- NULL

# ============================================================================
# 5c. SEM diagnostic plots
# ============================================================================
cat("\n--- Step 5c: SEM plots ---\n")

for (nm in list(list(res = semD, label = "biv",   title = "SEM Bivariate"),
                list(res = semF, label = "biv_kde", title = "SEM Bivariate+KDE"))) {
  tryCatch({
    if (!is.null(nm$res)) {
      sem_lab <- nm$res$adaptive$adaptive_labelling
      sem_post <- if (!is.null(sem_lab)) sem_lab[sem_lab$t >= 0, ] else NULL
      if (!is.null(sem_post) && nrow(sem_post) > 0) {
        sem_post$Process <- factor(sem_post$inferred_process,
                                   levels = c("control", "treated"))
        p_s <- ggplot(sem_post, aes(x = x, y = y, fill = Process, alpha = t)) +
          geom_point(size = 1.8, shape = 21, stroke = 0.3) +
          scale_fill_manual(values = c(control = "#377eb8", treated = "#e41a1c")) +
          scale_alpha_continuous(name = "Time (days)", range = c(0.3, 1)) +
          labs(title = paste0("Post-treatment (", nm$title, " labels)"),
               x = "X (km)", y = "Y (km)") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        store_analysis_plot(sprintf("pp_post_sem_%s", nm$label), p_s)
        cat(sprintf("  Stored results$plots$pp_post_sem_%s\n", nm$label))
      }
      p_f <- plot_flips(nm$res)
      store_analysis_plot(sprintf("sem_flips_%s", nm$label), p_f)
      cat(sprintf("  Stored results$plots$sem_flips_%s\n", nm$label))
    }
  }, error = function(e) cat(sprintf("  %s plot error: %s\n", nm$title, e$message)))
}

# ============================================================================
# 6. ATE estimation (post-treatment horizon)
# ============================================================================
cat(sprintf("\n--- Step 6: ATE estimation (%d-day simulation horizon) ---\n", ATE_WINDOW_DAYS))
cat("  Stage meaning: computing MAIN-FIT ATE payloads (A-J) now; sensitivity ATEs start at Step 7a.\n")

pp_post_sem_D <- if (!is.null(semD)) semD$adaptive$adaptive_labelling else NULL
pp_post_sem_D <- if (!is.null(pp_post_sem_D)) pp_post_sem_D[pp_post_sem_D$t >= 0, ] else pp_post
pp_post_sem_F <- if (!is.null(semF)) semF$adaptive$adaptive_labelling else NULL
pp_post_sem_F <- if (!is.null(pp_post_sem_F)) pp_post_sem_F[pp_post_sem_F$t >= 0, ] else pp_post
if (!is.null(pp_post_sem_D) && !is.null(pp_post_sem_F) && nrow(pp_post_sem_D) == nrow(pp_post_sem_F)) {
  same_labels_DF <- mean(pp_post_sem_D$inferred_process == pp_post_sem_F$inferred_process, na.rm = TRUE)
  d_flip <- sum(pp_post_sem_D$location_process != pp_post_sem_D$inferred_process, na.rm = TRUE)
  f_flip <- sum(pp_post_sem_F$location_process != pp_post_sem_F$inferred_process, na.rm = TRUE)
  cat(sprintf("  Relabelling overlap B vs D: %.1f%% same inferred labels (B flips=%d, D flips=%d)\n",
              100 * same_labels_DF, d_flip, f_flip))
}

cat(sprintf("  ATE window: [0, %d] days\n", ATE_WINDOW_DAYS))

B_marginals <- extract_marginals(B_params)
D_marginals <- extract_marginals(D_params)
E_marginals <- extract_marginals(E_params)
F_marginals <- extract_marginals(F_params)
G_marginals <- extract_marginals(G_params)
H_marginals <- extract_marginals(H_params)
I_marginals <- extract_marginals_indep(I_params)
J_marginals <- extract_marginals_indep(J_params)
K_marginals <- extract_marginals_indep(K_params)
L_marginals <- extract_marginals_indep(L_params)

ensure_ate_psock_pool()
# Compute both AoN and observed contrasts; attach under $by_contrast.
# Primary payload (backward-compatible fields) follows OK_ATE_CONTRAST.
ate_with_both_contrasts <- function(compute_one) {
  ate_aon <- tryCatch(
    compute_one("all_or_nothing"),
    error = function(e) {
      cat("    ATE all_or_nothing error:", e$message, "\n")
      NULL
    }
  )
  ate_obs <- tryCatch(
    compute_one("observed"),
    error = function(e) {
      cat("    ATE observed error:", e$message, "\n")
      NULL
    }
  )
  primary <- if (identical(OK_ATE_CONTRAST, "observed")) ate_obs else ate_aon
  if (!is.null(primary)) {
    primary$by_contrast <- list(
      all_or_nothing = ate_aon,
      observed = ate_obs
    )
  }
  primary
}
ate_biv_or_marginal <- function(biv_params, marg, observed_data, label,
                                bg_lookup = NULL,
                                state_spaces_in = state_spaces,
                                n_tiles_in = partition$n,
                                filtration_history_in = if (isTRUE(OK_ATE_CONDITIONAL_ON_PRE)) pp_pre else NULL) {
  if (isTRUE(SKIP_ATE)) {
    cat(sprintf("    SKIP_ATE: not evaluating ATE for %s.\n", label))
    return(NULL)
  }
  if (is.null(biv_params) || is.null(marg)) {
    cat(sprintf("    Skipping ATE for %s: fitted parameters are NULL.\n", label))
    return(NULL)
  }
  ate_with_both_contrasts(function(contrast) {
    if (isTRUE(OK_ATE_BIVARIATE) && !is.null(biv_params)) {
      ate_estim_bivariate(
        biv_params = biv_params,
        windowT = windowT_ate,
        windowS = win_km,
        state_spaces_obs = state_spaces_in,
        label = label,
        n_sims = ATE_N_SIMS,
        n_cores = ATE_SIM_CORES,
        m0 = ETAS_M0,
        beta_gr = BETA_GR,
        filtration_history = filtration_history_in,
        t_trunc = SEM_T_TRUNC_DAYS,
        n_tiles = n_tiles_in,
        crn_base_seed = ate_crn_base_seed,
        use_crn = OK_ATE_USE_CRN,
        crn_pair = OK_ATE_CRN_PAIR,
        parallel_cluster = ate_cl_reuse,
        quiet = FALSE,
        contrast = contrast,
        covariate_lookup = bg_lookup
      )
    } else {
      ate_estim_fast(
        marg$ctrl, marg$treat, observed_data, label,
        contrast = contrast,
        covariate_lookup = bg_lookup,
        t_trunc = SEM_T_TRUNC_DAYS
      )
    }
  })
}
ate_univ_both <- function(marg, observed_data, label, bg_lookup = NULL) {
  if (isTRUE(SKIP_ATE) || is.null(marg)) return(NULL)
  ate_with_both_contrasts(function(contrast) {
    ate_estim_fast(
      marg$ctrl, marg$treat, observed_data, label,
      contrast = contrast,
      covariate_lookup = bg_lookup,
      t_trunc = SEM_T_TRUNC_DAYS
    )
  })
}
if (isTRUE(CD_ONLY) || isTRUE(UNIV_KDE_ONLY)) {
  cat(sprintf("  %s: skipping ATE for homogeneous Fits A/B.\n",
              if (isTRUE(UNIV_KDE_ONLY)) "UNIV_KDE_ONLY" else "CD_ONLY"))
  ate_B <- NULL
  ate_D <- NULL
} else {
  t_ate_B <- proc.time()[["elapsed"]]
  ate_B <- ate_biv_or_marginal(B_params, B_marginals, pp_post,
                               "Fit B (naive bivariate)")
  ate_B_elapsed <- proc.time()[["elapsed"]] - t_ate_B
  add_timing_row("ate_B", ate_B_elapsed, if (!is.null(ate_B)) "ok" else "failed")
  t_ate_D <- proc.time()[["elapsed"]]
  ate_D <- ate_biv_or_marginal(D_params, D_marginals, pp_post_sem_D,
                               "Fit D (SEM bivariate)")
  ate_D_elapsed <- proc.time()[["elapsed"]] - t_ate_D
  add_timing_row("ate_D", ate_D_elapsed, if (!is.null(ate_D)) "ok" else "failed")
}
t_ate_E <- proc.time()[["elapsed"]]
ate_E <- if (isTRUE(UNIV_KDE_ONLY)) {
  cat("  UNIV_KDE_ONLY: skipping ATE for bivariate KDE Fits C/D.\n")
  NULL
} else {
  ate_biv_or_marginal(E_params, E_marginals, pp_post_bg,
                      "Fit C (naive biv+KDE, all-free)",
                      bg_lookup = KDE_BG_LOOKUP)
}
ate_E_elapsed <- proc.time()[["elapsed"]] - t_ate_E
add_timing_row("ate_E", ate_E_elapsed, if (!is.null(ate_E)) "ok" else if (isTRUE(UNIV_KDE_ONLY)) "skipped" else "failed")
t_ate_F <- proc.time()[["elapsed"]]
ate_F <- if (isTRUE(UNIV_KDE_ONLY)) {
  NULL
} else {
  ate_biv_or_marginal(F_params, F_marginals, pp_post_sem_F,
                      "Fit D (SEM biv+KDE, all-free)",
                      bg_lookup = KDE_BG_LOOKUP)
}
ate_F_elapsed <- proc.time()[["elapsed"]] - t_ate_F
add_timing_row("ate_F", ate_F_elapsed, if (!is.null(ate_F)) "ok" else if (isTRUE(UNIV_KDE_ONLY)) "skipped" else "failed")

# Keep kde_variant_fits$E/F$all_free$ate aligned with the canonical bivariate
# estimands. fits_named$C/D are built from these lists, so a later marginal
# overwrite would otherwise clobber ate_E/ate_F in report/paper outputs.
if (!is.null(kde_variant_fits$E$all_free)) {
  kde_variant_fits$E$all_free$ate <- ate_E
  kde_variant_fits$E$all_free$marginals <- E_marginals
}
if (!is.null(kde_variant_fits$F$all_free)) {
  kde_variant_fits$F$all_free$ate <- ate_F
  kde_variant_fits$F$all_free$marginals <- F_marginals
}

# Compute ATE/saved estimands for remaining KDE fix-profile pairs (E/F = productivity_free).
for (vid in names(kde_variant_specs)) {
  if (identical(vid, "all_free")) next
  if (!is.null(kde_variant_fits$E[[vid]])) {
    p_e <- kde_variant_fits$E[[vid]]$params
    m_e <- extract_marginals(p_e)
    letter_e <- KDE_FIT_LETTERS[[vid]]$E
    kde_variant_fits$E[[vid]]$marginals <- m_e
    kde_variant_fits$E[[vid]]$ate <- ate_biv_or_marginal(
      p_e, m_e, pp_post_bg,
      label = sprintf("Fit %s (naive biv+KDE, %s)", letter_e, vid),
      bg_lookup = KDE_BG_LOOKUP
    )
  }
  if (!is.null(kde_variant_fits$F[[vid]])) {
    p_f <- kde_variant_fits$F[[vid]]$params
    m_f <- extract_marginals(p_f)
    sem_obj <- kde_variant_fits$F[[vid]]$fit
    post_sem_local <- if (!is.null(sem_obj) && !is.null(sem_obj$adaptive$adaptive_labelling)) {
      sem_obj$adaptive$adaptive_labelling
    } else {
      pp_post_bg
    }
    post_sem_local <- post_sem_local[post_sem_local$t >= 0, , drop = FALSE]
    letter_f <- KDE_FIT_LETTERS[[vid]]$F
    kde_variant_fits$F[[vid]]$marginals <- m_f
    kde_variant_fits$F[[vid]]$ate <- ate_biv_or_marginal(
      p_f, m_f, post_sem_local,
      label = sprintf("Fit %s (SEM biv+KDE, %s)", letter_f, vid),
      bg_lookup = KDE_BG_LOOKUP
    )
  }
}
ate_I <- if (!is.null(I_marginals)) {
  t_ate_I <- proc.time()[["elapsed"]]
  ate_univ_both(
    I_marginals, pp_post_bg, "Fit G (Naive univariate ETAS + KDE)",
    bg_lookup = KDE_BG_LOOKUP
  )
} else NULL
if (exists("t_ate_I", inherits = FALSE)) {
  ate_I_elapsed <- proc.time()[["elapsed"]] - t_ate_I
  add_timing_row("ate_I", ate_I_elapsed, if (!is.null(ate_I)) "ok" else "failed")
}
ate_J <- if (!is.null(J_marginals)) {
  t_ate_J <- proc.time()[["elapsed"]]
  ate_univ_both(
    J_marginals, pp_post_sem_J, "Fit H (SEM univariate ETAS + KDE)",
    bg_lookup = KDE_BG_LOOKUP
  )
} else NULL
if (exists("t_ate_J", inherits = FALSE)) {
  ate_J_elapsed <- proc.time()[["elapsed"]] - t_ate_J
  add_timing_row("ate_J", ate_J_elapsed, if (!is.null(ate_J)) "ok" else "failed")
}
ate_K <- if (!is.null(K_marginals)) {
  t_ate_K <- proc.time()[["elapsed"]]
  ate_univ_both(K_marginals, pp_post, "Fit I (Naive univariate ETAS, no KDE)")
} else NULL
if (exists("t_ate_K", inherits = FALSE)) {
  ate_K_elapsed <- proc.time()[["elapsed"]] - t_ate_K
  add_timing_row("ate_K", ate_K_elapsed, if (!is.null(ate_K)) "ok" else "failed")
}
ate_L <- if (!is.null(L_marginals)) {
  t_ate_L <- proc.time()[["elapsed"]]
  ate_univ_both(L_marginals, pp_post_sem_L, "Fit J (SEM univariate ETAS, no KDE)")
} else NULL
if (exists("t_ate_L", inherits = FALSE)) {
  ate_L_elapsed <- proc.time()[["elapsed"]] - t_ate_L
  add_timing_row("ate_L", ate_L_elapsed, if (!is.null(ate_L)) "ok" else "failed")
}

} else {
  if (isTRUE(BOOTSTRAP_ONLY)) {
    cat("\n--- BOOTSTRAP_ONLY: skipped Steps 4-6 (county fits / main ATEs); hydrating C/D from patch ---\n")
    boot_patch <- readRDS(normalizePath(BOOTSTRAP_PATCH_FILE, winslash = "/", mustWork = TRUE))
    # Prefer new C/D (all_free); fall back to legacy E/F-as-primary or fitE/fitF aliases.
    E_params <- if (!is.null(boot_patch$fits_named$C$params)) {
      boot_patch$fits_named$C$params
    } else if (!is.null(boot_patch$fits_named$E$params)) {
      boot_patch$fits_named$E$params
    } else {
      boot_patch$fitE$params
    }
    F_params <- if (!is.null(boot_patch$fits_named$D$params)) {
      boot_patch$fits_named$D$params
    } else if (!is.null(boot_patch$fits_named$F$params)) {
      boot_patch$fits_named$F$params
    } else {
      boot_patch$fitF$params
    }
    if (!valid_biv_params(E_params) || !valid_biv_params(F_params)) {
      stop("BOOTSTRAP_ONLY: failed to hydrate C/D (all_free) params from patch.")
    }
    E_marginals <- extract_marginals(E_params)
    F_marginals <- extract_marginals(F_params)
    ate_E <- boot_patch$fits_named$C$ate
    if (is.null(ate_E)) ate_E <- boot_patch$fits_named$E$ate
    if (is.null(ate_E)) ate_E <- boot_patch$fitE$ate
    ate_F <- boot_patch$fits_named$D$ate
    if (is.null(ate_F)) ate_F <- boot_patch$fits_named$F$ate
    if (is.null(ate_F)) ate_F <- boot_patch$fitF$ate
    trunc_from_patch <- suppressWarnings(as.numeric(boot_patch$config$SEM_T_TRUNC_DAYS))
    if (length(trunc_from_patch) == 1L && is.finite(trunc_from_patch) && trunc_from_patch > 0) {
      SEM_T_TRUNC_DAYS <- trunc_from_patch
      SEM_T_TRUNC_SOURCE <- "patch"
      cat(sprintf("  Using SEM_T_TRUNC_DAYS=%.4f from patch config\n", SEM_T_TRUNC_DAYS))
    }
    cat(sprintf(
      "  Hydrated C/D (all_free) params for bivariate bootstrap (C ate_method=%s, D ate_method=%s)\n",
      if (!is.null(ate_E$ate_method)) ate_E$ate_method else "NA",
      if (!is.null(ate_F$ate_method)) ate_F$ate_method else "NA"
    ))
    rm(boot_patch)
  } else if (isTRUE(T_TRUNC_SENS_ONLY)) {
    cat("\n--- T_TRUNC_SENS_ONLY: skipped Steps 4-6 (county fits / main ATEs); will run Fit C/D t_trunc sensitivity only ---\n")
    trunc_patch <- readRDS(normalizePath(T_TRUNC_SENS_PATCH_FILE, winslash = "/", mustWork = TRUE))
    # Adopt patch CRN base if present; t_trunc is already pinned (auto disabled).
    crn_from_patch <- suppressWarnings(as.integer(trunc_patch$config$ATE_CRN_BASE))
    if (length(crn_from_patch) == 1L && is.finite(crn_from_patch) && !is.na(crn_from_patch)) {
      OK_ATE_CRN_BASE <- crn_from_patch
      ate_crn_base_seed <- as.integer(OK_ATE_CRN_BASE)
      cat(sprintf("  Using ATE_CRN_BASE=%d from patch config\n", OK_ATE_CRN_BASE))
    }
    rm(trunc_patch)
  } else {
    cat("\n--- FIT_VARIABILITY_ONLY: skipped Steps 4-6 (county fits, univ, partition sensitivities, main ATEs) ---\n")
  }
}

# ============================================================================
# 6a. Fit variability (repeat county all-free C/D fits; internal E/F objects)
# ============================================================================
fit_variability <- NULL
fit_variability_elapsed <- NA_real_
if (RUN_FIT_VARIABILITY && FIT_VARIABILITY_REPS > 0L) {
  t_fitvar <- proc.time()[["elapsed"]]
  cat(sprintf("\n--- Step 6a: Fit variability (C/D all_free repeats; reps=%d, cores=%d) ---\n",
              FIT_VARIABILITY_REPS, FIT_VARIABILITY_CORES))
  ensure_ate_psock_pool()

  fitvar_ate_stats <- function(ate_obj, n_tiles = partition$n) {
    out <- list(
      raw_total_saved = NA_real_,
      mc_total_saved_mean = NA_real_,
      mc_total_saved_sd = NA_real_,
      eta_ctrl = NA_real_,
      eta_treat = NA_real_
    )
    if (is.null(ate_obj)) return(out)
    if (!is.null(ate_obj$analytic)) {
      out$eta_ctrl <- suppressWarnings(as.numeric(ate_obj$analytic$eta_ctrl))
      out$eta_treat <- suppressWarnings(as.numeric(ate_obj$analytic$eta_treat))
      if (is.finite(out$eta_ctrl) && is.finite(out$eta_treat)) {
        out$raw_total_saved <- as.numeric(n_tiles) * (out$eta_ctrl - out$eta_treat)
      }
    }
    sim_total <- suppressWarnings(as.numeric(ate_obj$all_nothing_sim$total_saved))
    sim_total <- sim_total[is.finite(sim_total)]
    if (length(sim_total) > 0L) {
      out$mc_total_saved_mean <- mean(sim_total, na.rm = TRUE)
      out$mc_total_saved_sd <- if (length(sim_total) > 1L) stats::sd(sim_total, na.rm = TRUE) else 0
    }
    out
  }

  run_fitvar_rep <- function(rep_id) {
    t0_rep <- proc.time()[["elapsed"]]
    cat(sprintf("    [fitvar:%d] start pid=%d mem=%s\n",
                as.integer(rep_id), Sys.getpid(), mem_snapshot()))
    out <- list(rep = as.integer(rep_id))

    fitE_var <- tryCatch(
      fit_etas_bivariate(
        params_init = biv_init_E, realiz = pp_all_bg,
        windowT = windowT_fit, windowS = win_km, m0 = ETAS_M0,
        control_state_space = control_ss, treated_state_space = treated_ss,
        background_rate_var = "W",
        treated_background_zero_before = 0,
        control_background_everywhere_before = 0,
        control_background_pre_mass_ratio = CTRL_BG_PRE_MASS_RATIO,
        beta_gr = BETA_GR,
        enforce_finite_trigger_moments = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
        p_lower_bound = ETAS_P_LOWER_BOUND,
        q_lower_bound = ETAS_Q_LOWER_BOUND,
        max_branching_radius = ETAS_BRANCHING_MAX,
        maxit = VANILLA_MAXIT, fixed_params = SENSITIVITY_FIXED_PARAMS, trace = 0,
        t_trunc = SEM_T_TRUNC_DAYS
      ),
      error = function(e) NULL
    )
    E_params_var <- if (valid_biv_fit(fitE_var)) fitE_var$par else NULL
    E_marg_var <- extract_marginals(E_params_var)
    ate_E_var <- if (is.null(E_marg_var)) NULL else tryCatch(
      ate_estim_fast(
        E_marg_var$ctrl, E_marg_var$treat, pp_post_bg,
        label = sprintf("FitVar E #%d", rep_id),
        phase = "fit_variability",
        n_tiles_used = partition$n,
        treated_idx_used = treated_idx,
        quiet = TRUE,
        covariate_lookup = KDE_BG_LOOKUP,
        t_trunc = SEM_T_TRUNC_DAYS
      ),
      error = function(e) NULL
    )

    semF_var <- tryCatch(
      run_sem_fit(
        pp_data_in = pp_all_bg,
        partition_in = partition,
        partition_processes_in = partition_processes,
        state_spaces_in = state_spaces,
        init_params_in = biv_init_F,
        fixed_params_in = SENSITIVITY_FIXED_PARAMS,
        background_rate_var_in = "W",
        control_background_pre_mass_ratio_in = CTRL_BG_PRE_MASS_RATIO,
        sem_inner_iter_in = SEM_INNER_ITER,
        verbose_in = FALSE,
        label = sprintf("FitVar F #%d", rep_id)
      ),
      error = function(e) NULL
    )
    F_params_var <- if (valid_biv_sem(semF_var)) {
      semF_var$etas_bivariate_params
    } else {
      NULL
    }
    pp_post_sem_var <- if (!is.null(semF_var) && !is.null(semF_var$adaptive$adaptive_labelling)) {
      tmp <- semF_var$adaptive$adaptive_labelling
      tmp[tmp$t >= 0, , drop = FALSE]
    } else {
      pp_post_bg
    }
    n_relabel_var <- if (!is.null(semF_var) && !is.null(semF_var$adaptive$adaptive_labelling)) {
      lp <- semF_var$adaptive$adaptive_labelling
      lp <- lp[lp$t >= 0, , drop = FALSE]
      sum(lp$location_process != lp$inferred_process, na.rm = TRUE)
    } else { 0L }
    F_marg_var <- extract_marginals(F_params_var)
    ate_F_var <- if (is.null(F_marg_var)) NULL else tryCatch(
      ate_estim_fast(
        F_marg_var$ctrl, F_marg_var$treat, pp_post_sem_var,
        label = sprintf("FitVar F #%d", rep_id),
        phase = "fit_variability",
        n_tiles_used = partition$n,
        treated_idx_used = treated_idx,
        quiet = TRUE,
        covariate_lookup = KDE_BG_LOOKUP,
        t_trunc = SEM_T_TRUNC_DAYS
      ),
      error = function(e) NULL
    )

    elapsed_rep <- proc.time()[["elapsed"]] - t0_rep
    out$E <- list(
      ok = !is.null(E_params_var),
      params = E_params_var,
      ate = ate_E_var,
      ate_stats = fitvar_ate_stats(ate_E_var),
      loglik = if (!is.null(fitE_var) && !is.null(fitE_var$value)) as.numeric(fitE_var$value) else NA_real_
    )
    out$F <- list(
      ok = !is.null(F_params_var),
      params = F_params_var,
      ate = ate_F_var,
      ate_stats = fitvar_ate_stats(ate_F_var),
      n_relabel = as.integer(n_relabel_var)
    )
    out$elapsed_sec <- elapsed_rep
    cat(sprintf("    [fitvar:%d] done in %.1fs E_ok=%s F_ok=%s relabel=%d mem=%s\n",
                as.integer(rep_id), elapsed_rep, out$E$ok, out$F$ok,
                as.integer(out$F$n_relabel), mem_snapshot()))
    out
  }

  fitvar_jobs <- as.list(seq_len(FIT_VARIABILITY_REPS))
  fitvar_out <- if (FIT_VARIABILITY_CORES > 1L && length(fitvar_jobs) > 1L) {
    run_parallel(
      fitvar_jobs, run_fitvar_rep,
      cores = min(FIT_VARIABILITY_CORES, length(fitvar_jobs)),
      label = "fit-variability"
    )
  } else {
    lapply(fitvar_jobs, run_fitvar_rep)
  }

  fitvar_rows <- list()
  fitvar_param <- function(params, name) {
    if (is.null(params) || is.null(params[[name]])) return(NA_real_)
    value <- suppressWarnings(as.numeric(params[[name]]))
    if (length(value) != 1L || !is.finite(value)) NA_real_ else value
  }
  for (z in fitvar_out) {
    for (model_nm in c("E", "F")) {
      zi <- z[[model_nm]]
      if (is.null(zi)) next
      fitvar_rows[[length(fitvar_rows) + 1L]] <- data.frame(
        rep = as.integer(z$rep),
        model = as.character(model_nm),
        success = isTRUE(zi$ok),
        elapsed_sec = as.numeric(z$elapsed_sec),
        loglik = if (!is.null(zi$loglik)) as.numeric(zi$loglik) else NA_real_,
        n_relabel = if (identical(model_nm, "F")) {
          if (!is.null(zi$n_relabel)) as.integer(zi$n_relabel) else NA_integer_
        } else {
          NA_integer_
        },
        raw_total_saved = as.numeric(zi$ate_stats$raw_total_saved),
        mc_total_saved_mean = as.numeric(zi$ate_stats$mc_total_saved_mean),
        mc_total_saved_sd = as.numeric(zi$ate_stats$mc_total_saved_sd),
        eta_ctrl = as.numeric(zi$ate_stats$eta_ctrl),
        eta_treat = as.numeric(zi$ate_stats$eta_treat),
        mu_0 = fitvar_param(zi$params, "mu_0"),
        mu_1 = fitvar_param(zi$params, "mu_1"),
        A_00 = fitvar_param(zi$params, "A_00"),
        A_11 = fitvar_param(zi$params, "A_11"),
        alpha_m_00 = fitvar_param(zi$params, "alpha_m_00"),
        alpha_m_11 = fitvar_param(zi$params, "alpha_m_11"),
        c = fitvar_param(zi$params, "c"),
        p = fitvar_param(zi$params, "p"),
        D = fitvar_param(zi$params, "D"),
        gamma = fitvar_param(zi$params, "gamma"),
        q = fitvar_param(zi$params, "q"),
        stringsAsFactors = FALSE
      )
    }
  }
  fitvar_df <- if (length(fitvar_rows) > 0L) do.call(rbind, fitvar_rows) else data.frame()
  fit_variability <- list(
    config = list(
      reps = FIT_VARIABILITY_REPS,
      cores = FIT_VARIABILITY_CORES,
      sem_inner_iter = SEM_INNER_ITER,
      ate_n_sims = ATE_N_SIMS,
      seed = OK_GLOBAL_SEED
    ),
    replicate_summary = fitvar_df,
    replicates = fitvar_out
  )
  fit_variability_elapsed <- proc.time()[["elapsed"]] - t_fitvar
  add_timing_row(
    "fit_variability_total",
    fit_variability_elapsed,
    "ok",
    sprintf("reps=%d cores=%d", FIT_VARIABILITY_REPS, FIT_VARIABILITY_CORES)
  )
  if (nrow(fitvar_df) > 0L) {
    fitvar_ok <- fitvar_df[fitvar_df$success, , drop = FALSE]
    if (nrow(fitvar_ok) > 0L) {
      for (m in c("E", "F")) {
        mm <- fitvar_ok[fitvar_ok$model == m, , drop = FALSE]
        if (nrow(mm) < 1L) next
        cat(sprintf(
          "  FitVar %s: reps=%d raw_total_saved=%.1f +/- %.1f, mc_total_saved=%.1f +/- %.1f\n",
          m, nrow(mm),
          mean(mm$raw_total_saved, na.rm = TRUE), stats::sd(mm$raw_total_saved, na.rm = TRUE),
          mean(mm$mc_total_saved_mean, na.rm = TRUE), stats::sd(mm$mc_total_saved_mean, na.rm = TRUE)
        ))
      }
    }
  }
  if (nzchar(FIT_VARIABILITY_PATCH_FILE)) {
    patch_file <- normalizePath(FIT_VARIABILITY_PATCH_FILE, winslash = "/", mustWork = FALSE)
    patch_ok <- FALSE
    if (file.exists(patch_file)) {
      patched <- tryCatch(readRDS(patch_file), error = function(e) NULL)
      if (!is.null(patched)) {
        patched$fit_variability <- fit_variability
        if (is.null(patched$config)) patched$config <- list()
        patched$config$RUN_FIT_VARIABILITY <- TRUE
        patched$config$FIT_VARIABILITY_REPS <- FIT_VARIABILITY_REPS
        patched$config$FIT_VARIABILITY_CORES <- FIT_VARIABILITY_CORES
        patched$config$FIT_VARIABILITY_ONLY <- isTRUE(FIT_VARIABILITY_ONLY)
        saveRDS(patched, patch_file)
        patch_ok <- TRUE
        cat(sprintf("Patched fit variability into existing results: %s\n", patch_file))
      } else {
        cat(sprintf("Fit variability patch skipped: failed to read %s\n", patch_file))
      }
    } else {
      cat(sprintf("Fit variability patch skipped: file does not exist: %s\n", patch_file))
    }
    if (isTRUE(FIT_VARIABILITY_ONLY) && !isTRUE(patch_ok)) {
      stop("FIT_VARIABILITY_ONLY: failed to merge fit_variability into ", patch_file)
    }
  }
} else {
  add_timing_row("fit_variability_total", NA_real_, "skipped", "disabled")
}

if (isTRUE(FIT_VARIABILITY_ONLY)) {
  cat("\n=== FIT_VARIABILITY_ONLY: finished; exiting before Step 6b checkpoint / sensitivity / bootstrap ===\n")
  quit(save = "no", status = 0)
}

if (isTRUE(BOOTSTRAP_ONLY)) {
  ensure_ate_psock_pool()
  cat("\n--- BOOTSTRAP_ONLY: skipping Step 6b/7 sensitivity checkpoints; jumping to Step 8 bootstrap ---\n")
}
t_trunc_sensitivity <- NULL

# Save checkpoint with all main-fit ATE payloads, but without sensitivity payloads.
if (!isTRUE(BOOTSTRAP_ONLY)) {
if (!isTRUE(T_TRUNC_SENS_ONLY)) {
cat("\n--- Step 6b checkpoint: saving pre-sensitivity results ---\n")
pre_sens_saved_at <- as.character(Sys.time())
# Public A–J lettering. Legacy top-level fitE/fitF alias the primary all_free
# pair (fits_named$C/$D) for older checkpoint / bootstrap-only code paths.
fits_named_pre_sensitivity <- list(
  A = list(
    letter = "A", label = sprintf("Naive bivariate (%s)", PRIMARY_PARTITION_LABEL),
    method = "county_bivariate", algorithm = "naive",
    params = B_params, fit = fitB, ate = ate_B
  ),
  B = list(
    letter = "B", label = sprintf("SEM bivariate (%s)", PRIMARY_PARTITION_LABEL),
    method = "county_bivariate", algorithm = "sem",
    params = D_params, fit = semD, ctrl = D_ctrl, treat = D_treat, ate = ate_D
  ),
  C = c(list(letter = "C", label = "Naive biv+KDE (all_free)",
             method = "kde_all_free", algorithm = "naive"),
        kde_variant_fits$E$all_free),
  D = c(list(letter = "D", label = "SEM biv+KDE (all_free)",
             method = "kde_all_free", algorithm = "sem"),
        kde_variant_fits$F$all_free),
  E = c(list(letter = "E", label = "Naive biv+KDE (productivity_free)",
             method = "kde_productivity_free", algorithm = "naive"),
        kde_variant_fits$E$productivity_free),
  F = c(list(letter = "F", label = "SEM biv+KDE (productivity_free)",
             method = "kde_productivity_free", algorithm = "sem"),
        kde_variant_fits$F$productivity_free),
  G = list(
    letter = "G", label = sprintf("Naive univariate ETAS + KDE (%s)", PRIMARY_PARTITION_LABEL),
    method = "county_univariate_kde", algorithm = "naive",
    params = I_params, fit = fitI, ate = ate_I
  ),
  H = list(
    letter = "H", label = sprintf("SEM univariate ETAS + KDE (%s)", PRIMARY_PARTITION_LABEL),
    method = "county_univariate_kde", algorithm = "sem",
    params = J_params, fit = semJ, ctrl = J_ctrl, treat = J_treat, ate = ate_J
  ),
  I = list(
    letter = "I", label = sprintf("Naive univariate ETAS (%s, no KDE)", PRIMARY_PARTITION_LABEL),
    method = "county_univariate_homog", algorithm = "naive",
    params = K_params, fit = fitK, ate = ate_K
  ),
  J = list(
    letter = "J", label = sprintf("SEM univariate ETAS (%s, no KDE)", PRIMARY_PARTITION_LABEL),
    method = "county_univariate_homog", algorithm = "sem",
    params = L_params, fit = semL, ctrl = L_ctrl, treat = L_treat, ate = ate_L
  )
)
results_pre_sensitivity <- list(
  fitB = list(params = B_params, loglik = B_loglik, fit = fitB, ate = ate_B),
  fitD = list(params = D_params, ctrl = D_ctrl, treat = D_treat, sem = semD, ate = ate_D),
  # Legacy aliases: fitE/fitF hold primary all_free content (= fits_named$C/$D).
  fitE = list(params = E_params, loglik = E_loglik, fit = fitE, ate = ate_E),
  fitF = list(params = F_params, ctrl = F_ctrl, treat = F_treat, sem = semF, ate = ate_F),
  fitI = list(params = I_params, fit = fitI, ate = ate_I),
  fitJ = list(params = J_params, fit = semJ, ctrl = J_ctrl, treat = J_treat, ate = ate_J),
  fitK = list(params = K_params, fit = fitK, ate = ate_K),
  fitL = list(params = L_params, fit = semL, ctrl = L_ctrl, treat = L_treat, ate = ate_L),
  fits_named = fits_named_pre_sensitivity,
  kde_variants = kde_variant_fits,
  fitG = NULL,
  fitH = NULL,
  bootstrap_ate = NULL,
  fit_variability = fit_variability,
  partition_results = NULL,
  ate_partitions = NULL,
  kde_bandwidth_sensitivity = NULL,
  t_trunc_sensitivity = NULL,
  pp_data = list(pp_pre = pp_pre, pp_pre_holdout = pp_pre_holdout, pp_post = pp_post),
  kde_info = kde_info,
  control_snapshot_fits = control_snapshot_fits,
  sem_pilot_tuning = sem_pilot_tuning,
  plots = analysis_plots,
  counties = list(
    names = counties_sf_valid$NAME,
    treated_names = county_treated_names,
    n_counties = county_partition$n,
    n_treated = sum(county_treated_idx)
  ),
  metadata = list(stage = "pre_sensitivity", saved_at = pre_sens_saved_at),
  checkpoint = list(
    stage = "pre_sensitivity",
    saved_at = pre_sens_saved_at
  ),
  fit_name_map = list(
    A = "fits_named$A",
    B = "fits_named$B",
    C = "fits_named$C",
    D = "fits_named$D",
    E = "fits_named$E",
    F = "fits_named$F",
    G = "fits_named$G",
    H = "fits_named$H",
    I = "fits_named$I",
    J = "fits_named$J"
  ),
  config = list(
    ETAS_M0 = ETAS_M0, BETA_GR = BETA_GR,
    BETA_GR_SOURCE = BETA_GR_SOURCE,
    GAMMA_FIXED = GAMMA_FIXED,
    FIXED_STRUCTURAL = FIXED_STRUCTURAL,
    SEM_N_ITER = SEM_N_ITER, SEM_INNER_ITER = SEM_INNER_ITER,
    SENS_SEM_INNER_ITER = SENS_SEM_INNER_ITER,
    SEM_INNER_PROPS = SEM_INNER_PROPS,
    SEM_N_LABELLINGS = SEM_N_LABELLINGS,
    SEM_OUTER_MAXIT = SEM_OUTER_MAXIT,
    SEM_OUTER_MAXIT_BIV = SEM_OUTER_MAXIT_BIV,
    SEM_WARMSTART_FIXED = SEM_WARMSTART_FIXED,
    SEM_CHANGE_FACTOR = SEM_CHANGE_FACTOR,
    SEM_CHANGE_FACTOR_MIN_MULT = SEM_CHANGE_FACTOR_MIN_MULT,
    SEM_CHANGE_FACTOR_MAX_MULT = SEM_CHANGE_FACTOR_MAX_MULT,
    SEM_OPTIM_METHOD = SEM_OPTIM_METHOD,
    SEM_MONOTONE_COMPLETE_LL = SEM_MONOTONE_COMPLETE_LL,
    SEM_START_FROM_C = SEM_START_FROM_C,
    SEM_BIV_N_THREADS = SEM_BIV_N_THREADS,
    SEM_SINGLE_FLIP_FROM_ITER = SEM_SINGLE_FLIP_FROM_ITER,
    SEM_SELECTION_TEMPERATURE = SEM_SELECTION_TEMPERATURE,
    RUN_SEM_PILOT = RUN_SEM_PILOT,
    SEM_PILOT_INNER_ITER = SEM_PILOT_INNER_ITER,
    SEM_PILOT_CORES = SEM_PILOT_CORES,
    SEM_PILOT_MAX_COMBOS = SEM_PILOT_MAX_COMBOS,
    SEM_TEMPORAL_WEIGHT = SEM_TEMPORAL_WEIGHT,
    SEM_TEMPORAL_SCALE_DAYS = SEM_TEMPORAL_SCALE_DAYS,
    SEM_T_TRUNC_DAYS = SEM_T_TRUNC_DAYS,
    SEM_T_TRUNC_SOURCE = SEM_T_TRUNC_SOURCE,
    SEM_T_TRUNC_REL = SEM_T_TRUNC_REL,
    T_TRUNC_SENS_DAYS = T_TRUNC_SENS_DAYS,
    RUN_T_TRUNC_SENSITIVITY = RUN_T_TRUNC_SENSITIVITY,
    RUN_SENSITIVITY = RUN_SENSITIVITY,
    RUN_FIT_VARIABILITY = RUN_FIT_VARIABILITY,
    FIT_VARIABILITY_REPS = FIT_VARIABILITY_REPS,
    FIT_VARIABILITY_CORES = FIT_VARIABILITY_CORES,
    KDE_VARIANT_MODE = KDE_VARIANT_MODE,
    KDE_BW_METHOD = KDE_BW_METHOD,
    PRIMARY_PARTITION = PRIMARY_PARTITION,
    KDE_BW_SENS_KM = KDE_BW_SENS_KM,
    RUN_PARTITION_SENSITIVITY = RUN_PARTITION_SENSITIVITY,
    CD_ONLY = CD_ONLY,
    UNIV_KDE_ONLY = UNIV_KDE_ONLY,
    RUN_KDE_PROFILE_SWEEP = RUN_KDE_PROFILE_SWEEP,
    MEMORY_SAFE = MEMORY_SAFE,
    TRIM_SENS_OBJECTS = TRIM_SENS_OBJECTS,
    SENS_CORES = SENS_CORES,
    SENS_BIV_N_THREADS = SENS_BIV_N_THREADS,
    ATE_SIM_CORES = ATE_SIM_CORES,
    ATE_USE_CRN = OK_ATE_USE_CRN,
    ATE_CRN_PAIR = OK_ATE_CRN_PAIR,
    ATE_CONDITIONAL_ON_PRE = OK_ATE_CONDITIONAL_ON_PRE,
    ATE_CRN_BASE = OK_ATE_CRN_BASE,
    SEM_MAX_RELABEL_STEP_FRAC = SEM_MAX_RELABEL_STEP_FRAC,
    SEM_FORCE_PARAM_UPDATE_FLIP_FRAC = SEM_FORCE_PARAM_UPDATE_FLIP_FRAC,
    ATE_N_SIMS = ATE_N_SIMS, ATE_WINDOW_DAYS = ATE_WINDOW_DAYS,
    ATE_BIVARIATE = OK_ATE_BIVARIATE, ATE_CONTRAST = OK_ATE_CONTRAST,
    RUN_BOOTSTRAP_ATE = RUN_BOOTSTRAP_ATE,
    BOOT_N_REPS = BOOT_N_REPS,
    BOOT_TARGETS = BOOT_TARGETS,
    SENS_TARGETS = SENS_TARGETS,
    BOOT_REFIT_SCOPE = BOOT_REFIT_SCOPE,
    BOOT_SEM_INNER_ITER = BOOT_SEM_INNER_ITER,
    BOOT_OUTER_CORES = BOOT_OUTER_CORES,
    BOOT_SEED = BOOT_SEED,
    BOOT_IDENTICAL_RANDOMNESS = OK_BOOT_IDENTICAL_RANDOMNESS,
    BOOT_GUARD_DEGENERATE = OK_BOOT_GUARD_DEGENERATE,
    ETAS_BRANCHING_MAX = ETAS_BRANCHING_MAX,
    ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
    ETAS_P_LOWER_BOUND = ETAS_P_LOWER_BOUND,
    ETAS_Q_LOWER_BOUND = ETAS_Q_LOWER_BOUND,
    BOOT_BRANCHING_MAX = BOOT_BRANCHING_MAX,
    BOOT_MAX_PRE_EVENTS = BOOT_MAX_PRE_EVENTS,
    BOOT_MAX_POST_EVENTS_PER_PROC = BOOT_MAX_POST_EVENTS_PER_PROC,
    BOOT_MAX_TOTAL_EVENTS = BOOT_MAX_TOTAL_EVENTS,
    TEST_MODE = TEST_MODE,
    windowT_post = windowT_post,
    n_pre = nrow(pp_pre), n_pre_holdout = nrow(pp_pre_holdout), n_pre_total = nrow(pp_pre_all),
    n_post = nrow(pp_post),
    n_counties = county_partition$n,
    n_treated_counties = sum(county_treated_idx),
    n_tiles = partition$n,
    n_treated = sum(treated_idx),
    grid_diameters = NULL,
    grid_multipliers = NULL,
    kde_bandwidth_multipliers = NULL,
    spatial_scale_D = trigger_range_km
  )
)
pre_sensitivity_out_file <- file.path(OUT_DIR, add_file_tag("oklahoma_results_pre_sensitivity.rds"))
saveRDS(results_pre_sensitivity, pre_sensitivity_out_file)
cat(sprintf("Pre-sensitivity checkpoint saved to: %s\n", pre_sensitivity_out_file))
rm(results_pre_sensitivity)
invisible(gc(verbose = FALSE))

# Execute sensitivity runs after main-fit SEM plots + ATE payloads so
# the final two major stages are sensitivities, then bootstraps.
cat("\n--- Step 7: Sensitivity dispatch (bandwidth + partition) ---\n")
sens_dispatch <- run_sensitivity_dispatch()
kde_bandwidth_fits <- sens_dispatch$kde_bandwidth_fits
partition_results <- sens_dispatch$partition_results
rm(sens_dispatch)
invisible(gc(verbose = FALSE))

# Sensitivity ATE payloads are computed before pre-bootstrap checkpoint so
# pre-bootstrap results can fully render partition/bandwidth sections.
cat("\n--- Step 7a: Sensitivity ATE payloads (for checkpoint + final report) ---\n")
cat("  Stage meaning: MAIN-FIT ATE stage is complete; now computing sensitivity ATEs only.\n")

# ATE sensitivity by KDE bandwidth (primary partition, inhomogeneous C/D all_free)
kde_bandwidth_sensitivity <- lapply(kde_bandwidth_fits, function(kf) {
  if (is.null(kf) || (is.null(kf$E_params) && is.null(kf$F_params))) return(NULL)
  ate_saved_stats <- function(ate_obj) {
    sim_saved <- if (!is.null(ate_obj) && !is.null(ate_obj$all_nothing_sim$total_saved)) {
      suppressWarnings(as.numeric(ate_obj$all_nothing_sim$total_saved))
    } else {
      numeric(0)
    }
    sim_saved <- sim_saved[is.finite(sim_saved)]
    list(
      mean = if (length(sim_saved) > 0L) mean(sim_saved) else NA_real_,
      sd = if (length(sim_saved) > 1L) stats::sd(sim_saved) else NA_real_,
      n = length(sim_saved)
    )
  }
  ate_E <- if (!is.null(kf$E_params)) {
    tryCatch(
      ate_biv_or_marginal(
        kf$E_params, extract_marginals(kf$E_params), kf$pp_post_bg,
        sprintf("BW %s Fit C naive biv+KDE", kf$label),
        bg_lookup = KDE_BG_LOOKUP
      ),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  ate_F <- if (!is.null(kf$F_params)) {
    tryCatch(
      ate_biv_or_marginal(
        kf$F_params, extract_marginals(kf$F_params), kf$pp_post_sem,
        sprintf("BW %s Fit D SEM biv+KDE", kf$label),
        bg_lookup = KDE_BG_LOOKUP
      ),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  stats_E <- ate_saved_stats(ate_E)
  stats_F <- ate_saved_stats(ate_F)
  list(
    label = kf$label,
    multiplier = kf$multiplier,
    sigma = kf$sigma,
    E_params = kf$E_params,
    F_params = kf$F_params,
    n_relabel = kf$n_relabel,
    ate_E = ate_E,
    ate_F = ate_F,
    C_ate_mean_saved = stats_E$mean,
    C_ate_sd_saved = stats_E$sd,
    C_ate_n_sims = stats_E$n,
    # Legacy names remain Fit D so older reports still read D.
    ate_mean_saved = stats_F$mean,
    ate_sd_saved = stats_F$sd,
    ate_n_sims = stats_F$n
  )
})
kde_bandwidth_sensitivity <- Filter(Negate(is.null), kde_bandwidth_sensitivity)

# ATE for alternative partitions (C/D all-free, bivariate AoN)
ate_partitions <- lapply(partition_results, function(pr) {
  if (is.null(pr) || identical(pr$label, PRIMARY_PARTITION_LABEL)) return(NULL)
  if (is.null(pr$E_params) && is.null(pr$F_params)) return(NULL)
  part_info <- all_partitions[[pr$label]]
  if (is.null(part_info)) return(NULL)
  n_tiles_p <- if (!is.null(pr$n_tiles)) pr$n_tiles else part_info$partition$n
  bg_lookup_p <- list(
    control = make_kde_bg_lookup(part_info$state_spaces$control),
    treated = make_kde_bg_lookup(part_info$state_spaces$treated)
  )
  pp_pre_part <- if (isTRUE(OK_ATE_CONDITIONAL_ON_PRE)) {
    assign_to_partition(pp_pre[, c("x", "y", "t", "mag")], part_info)
  } else {
    NULL
  }
  ate_E <- if (!is.null(pr$E_params)) {
    tryCatch(
      ate_biv_or_marginal(
        pr$E_params, extract_marginals(pr$E_params), NULL,
        sprintf("%s naive biv+KDE", pr$label),
        bg_lookup = bg_lookup_p,
        state_spaces_in = part_info$state_spaces,
        n_tiles_in = n_tiles_p,
        filtration_history_in = pp_pre_part
      ),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  ate_F <- if (!is.null(pr$F_params)) {
    tryCatch(
      ate_biv_or_marginal(
        pr$F_params, extract_marginals(pr$F_params), NULL,
        sprintf("%s SEM biv+KDE", pr$label),
        bg_lookup = bg_lookup_p,
        state_spaces_in = part_info$state_spaces,
        n_tiles_in = n_tiles_p,
        filtration_history_in = pp_pre_part
      ),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  list(ate_E = ate_E, ate_F = ate_F)
})
ate_partitions <- Filter(Negate(is.null), ate_partitions)
} else {
  cat("\n--- T_TRUNC_SENS_ONLY: skipping Step 6b checkpoint and bandwidth/partition sensitivity ---\n")
  kde_bandwidth_sensitivity <- NULL
  partition_results <- NULL
  ate_partitions <- NULL
  ensure_ate_psock_pool()
}

# ============================================================================
# 7a2. Fit C/D temporal-truncation sensitivity (ATE vs t_trunc)
#      Optional: multi-seed smoke at fixed t_trunc (no ATE) when OK_SMOKE_SEM_D_SEEDS > 0
# ============================================================================
smoke_sem_d <- NULL
t_trunc_sensitivity <- NULL
if (SMOKE_SEM_D_SEEDS > 0L) {
  cat(sprintf("\n--- Step 7a2-smoke: Fit D SEM multi-seed at t_trunc=%.4g (%d seeds, no ATE) ---\n",
              SMOKE_SEM_D_TRUNC, SMOKE_SEM_D_SEEDS))
  biv_init_trunc <- apply_pre_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
  cap_tol <- 1e-4
  run_smoke_d_seed <- function(seed_i) {
    seed_i <- as.integer(seed_i)
    t0 <- proc.time()[["elapsed"]]
    rng_label <- sprintf("smoke_d_trunc%.4g_seed%d", SMOKE_SEM_D_TRUNC, seed_i)
    cat(sprintf("    [smoke_d:%s] start pid=%d mem=%s\n",
                rng_label, Sys.getpid(), mem_snapshot()))
    sem_local <- tryCatch({
      run_sem_fit(
        pp_data_in = pp_all_bg,
        partition_in = partition,
        partition_processes_in = partition_processes,
        state_spaces_in = state_spaces,
        init_params_in = biv_init_trunc,
        fixed_params_in = SENSITIVITY_FIXED_PARAMS,
        background_rate_var_in = "W",
        control_background_pre_mass_ratio_in = CTRL_BG_PRE_MASS_RATIO,
        sem_t_trunc_in = SMOKE_SEM_D_TRUNC,
        sem_inner_iter_in = SEM_INNER_ITER,
        verbose_in = FALSE,
        sem_rng_label_in = rng_label,
        label = sprintf("Fit D smoke trunc %.4g seed %d", SMOKE_SEM_D_TRUNC, seed_i)
      )
    }, error = function(e) {
      cat(sprintf("    [smoke_d:%s] SEM error: %s\n", rng_label, e$message))
      NULL
    })
    params_local <- if (valid_biv_sem(sem_local)) sem_local$etas_bivariate_params else NULL
    rho_local <- if (!is.null(params_local)) {
      tryCatch(
        PPDisentangle:::.etas_biv_spectral_radius(params_local, BETA_GR),
        error = function(e) NA_real_
      )
    } else {
      NA_real_
    }
    is_explosive <- !is.null(params_local) && (
      !is.finite(rho_local) ||
        rho_local >= 1 ||
        rho_local >= (ETAS_BRANCHING_MAX - cap_tol)
    )
    status_local <- if (is.null(params_local)) {
      "failed"
    } else if (isTRUE(is_explosive)) {
      "explosive"
    } else {
      "ok"
    }
    elapsed <- proc.time()[["elapsed"]] - t0
    cat(sprintf("    [smoke_d:%s] done in %.1fs status=%s rho=%s mem=%s\n",
                rng_label, elapsed, status_local,
                format(rho_local, digits = 6), mem_snapshot()))
    list(
      seed = seed_i,
      rng_label = rng_label,
      t_trunc_days = as.numeric(SMOKE_SEM_D_TRUNC),
      status = status_local,
      elapsed_sec = as.numeric(elapsed),
      rho = as.numeric(rho_local),
      branching_metric = as.numeric(rho_local),
      params = params_local,
      n_relabel = if (!is.null(sem_local) && !is.null(sem_local$adaptive$adaptive_labelling)) {
        lp <- sem_local$adaptive$adaptive_labelling
        lp <- lp[lp$t >= 0, , drop = FALSE]
        sum(lp$location_process != lp$inferred_process, na.rm = TRUE)
      } else {
        NA_integer_
      }
    )
  }
  smoke_jobs <- as.list(seq_len(SMOKE_SEM_D_SEEDS))
  smoke_sem_d_runs <- if (N_CORES > 1L && length(smoke_jobs) > 1L) {
    run_parallel(
      smoke_jobs, run_smoke_d_seed,
      cores = min(length(smoke_jobs), max(1L, N_CORES)),
      label = "smoke_sem_d"
    )
  } else {
    lapply(smoke_jobs, run_smoke_d_seed)
  }
  names(smoke_sem_d_runs) <- sprintf("seed_%d", seq_len(SMOKE_SEM_D_SEEDS))
  n_ok <- sum(vapply(smoke_sem_d_runs, function(x) identical(x$status, "ok"), logical(1)))
  n_exp <- sum(vapply(smoke_sem_d_runs, function(x) identical(x$status, "explosive"), logical(1)))
  n_fail <- sum(vapply(smoke_sem_d_runs, function(x) identical(x$status, "failed"), logical(1)))
  rhos <- vapply(smoke_sem_d_runs, function(x) as.numeric(x$rho), numeric(1))
  cat(sprintf(
    "  Fit D smoke complete: ok=%d explosive=%d failed=%d / %d | rho range=[%s, %s]\n",
    as.integer(n_ok), as.integer(n_exp), as.integer(n_fail), length(smoke_sem_d_runs),
    format(min(rhos, na.rm = TRUE), digits = 6),
    format(max(rhos, na.rm = TRUE), digits = 6)
  ))
  smoke_sem_d <- list(
    t_trunc_days = SMOKE_SEM_D_TRUNC,
    n_seeds = SMOKE_SEM_D_SEEDS,
    sem_inner_iter = SEM_INNER_ITER,
    branching_max = ETAS_BRANCHING_MAX,
    n_ok = as.integer(n_ok),
    n_explosive = as.integer(n_exp),
    n_failed = as.integer(n_fail),
    rho = rhos,
    runs = smoke_sem_d_runs
  )
  patch_file <- normalizePath(T_TRUNC_SENS_PATCH_FILE, winslash = "/", mustWork = TRUE)
  out_file <- file.path(OUT_DIR, add_file_tag("oklahoma_smoke_sem_d.rds"))
  saveRDS(list(
    smoke_sem_d = smoke_sem_d,
    metadata = list(
      stage = "smoke_sem_d",
      saved_at = as.character(Sys.time()),
      patched_from = patch_file,
      sem_inner_iter = SEM_INNER_ITER,
      t_trunc_days = SMOKE_SEM_D_TRUNC,
      n_seeds = SMOKE_SEM_D_SEEDS
    ),
    config = list(
      SMOKE_SEM_D_SEEDS = SMOKE_SEM_D_SEEDS,
      SMOKE_SEM_D_TRUNC = SMOKE_SEM_D_TRUNC,
      SEM_INNER_ITER = SEM_INNER_ITER,
      SEM_OUTER_MAXIT_BIV = SEM_OUTER_MAXIT_BIV,
      ETAS_BRANCHING_MAX = ETAS_BRANCHING_MAX,
      ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
      ETAS_P_LOWER_BOUND = ETAS_P_LOWER_BOUND,
      ETAS_Q_LOWER_BOUND = ETAS_Q_LOWER_BOUND
    )
  ), out_file)
  cat(sprintf(
    "\n=== SMOKE_SEM_D: wrote %s; exiting ===\n",
    out_file
  ))
  quit(save = "no", status = 0)
}

cat("\n--- Step 7a2: Fit C/D t_trunc sensitivity ---\n")
if (length(T_TRUNC_SENS_DAYS) < 1L) {
  cat("  Skipping Fit C/D t_trunc sensitivity (empty grid / disabled).\n")
} else {
  cat(sprintf("  Refitting Fit C (naive) and Fit D (SEM) + AoN ATE at t_trunc days: %s\n",
              paste(signif(T_TRUNC_SENS_DAYS, 4), collapse = ", ")))
  biv_init_trunc <- apply_pre_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
  run_t_trunc_sens_job <- function(t_trunc_val) {
    t0 <- proc.time()[["elapsed"]]
    label <- sprintf("t_trunc=%.4g", t_trunc_val)
    cat(sprintf("    [t_trunc_sens:%s] start pid=%d mem=%s\n",
                label, Sys.getpid(), mem_snapshot()))
    cap_tol <- 1e-4
    rho_of <- function(params) {
      if (is.null(params)) return(NA_real_)
      tryCatch(
        PPDisentangle:::.etas_biv_spectral_radius(params, BETA_GR),
        error = function(e) NA_real_
      )
    }
    is_explosive_rho <- function(params, rho_val) {
      !is.null(params) && (
        !is.finite(rho_val) ||
          rho_val >= 1 ||
          rho_val >= (ETAS_BRANCHING_MAX - cap_tol)
      )
    }
    ate_at_trunc <- function(params, fit_label) {
      ate_estim_bivariate(
        biv_params = params,
        windowT = windowT_ate,
        windowS = win_km,
        state_spaces_obs = state_spaces,
        label = fit_label,
        n_sims = ATE_N_SIMS,
        n_cores = 1L,
        m0 = ETAS_M0,
        beta_gr = BETA_GR,
        filtration_history = if (isTRUE(OK_ATE_CONDITIONAL_ON_PRE)) pp_pre else NULL,
        t_trunc = t_trunc_val,
        n_tiles = partition$n,
        crn_base_seed = ate_crn_base_seed,
        use_crn = OK_ATE_USE_CRN,
        crn_pair = OK_ATE_CRN_PAIR,
        quiet = TRUE,
        contrast = "all_or_nothing",
        covariate_lookup = KDE_BG_LOOKUP
      )
    }
    saved_stats <- function(ate_obj) {
      sim_saved <- if (!is.null(ate_obj) && !is.null(ate_obj$all_nothing_sim$total_saved)) {
        suppressWarnings(as.numeric(ate_obj$all_nothing_sim$total_saved))
      } else {
        numeric(0)
      }
      sim_saved <- sim_saved[is.finite(sim_saved)]
      list(
        mean = if (length(sim_saved) > 0L) mean(sim_saved) else NA_real_,
        sd = if (length(sim_saved) > 1L) stats::sd(sim_saved) else NA_real_,
        n = length(sim_saved)
      )
    }
    status_of <- function(params, is_explosive, n_sims) {
      if (isTRUE(is_explosive)) {
        "explosive"
      } else if (!is.null(params) && n_sims > 0L) {
        "ok"
      } else {
        "failed"
      }
    }

    fitC_local <- tryCatch({
      fit_etas_bivariate(
        params_init = biv_init_trunc, realiz = pp_all_bg,
        windowT = windowT_fit, windowS = win_km, m0 = ETAS_M0,
        control_state_space = control_ss, treated_state_space = treated_ss,
        background_rate_var = "W",
        treated_background_zero_before = 0,
        control_background_everywhere_before = 0,
        control_background_pre_mass_ratio = CTRL_BG_PRE_MASS_RATIO,
        beta_gr = BETA_GR,
        enforce_finite_trigger_moments = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
        p_lower_bound = ETAS_P_LOWER_BOUND,
        q_lower_bound = ETAS_Q_LOWER_BOUND,
        max_branching_radius = ETAS_BRANCHING_MAX,
        maxit = VANILLA_MAXIT, fixed_params = SENSITIVITY_FIXED_PARAMS, trace = 0,
        t_trunc = t_trunc_val
      )
    }, error = function(e) {
      cat(sprintf("    [t_trunc_sens:%s] Fit C error: %s\n", label, e$message))
      NULL
    })
    c_params_local <- if (valid_biv_fit(fitC_local)) fitC_local$par else NULL
    c_rho_local <- rho_of(c_params_local)
    c_explosive <- is_explosive_rho(c_params_local, c_rho_local)
    ate_c_local <- NULL
    if (!is.null(c_params_local) && !isTRUE(c_explosive)) {
      ate_c_local <- tryCatch({
        ate_at_trunc(c_params_local, sprintf("Fit C t_trunc %.4g", t_trunc_val))
      }, error = function(e) {
        cat(sprintf("    [t_trunc_sens:%s] Fit C ATE error: %s\n", label, e$message))
        NULL
      })
    } else if (isTRUE(c_explosive)) {
      cat(sprintf(
        "    [t_trunc_sens:%s] discarding explosive/near-critical Fit C (rho=%s, margin=%.3f)\n",
        label, format(c_rho_local, digits = 6), ETAS_BRANCHING_MAX
      ))
    }
    c_saved <- saved_stats(ate_c_local)
    c_status_local <- status_of(c_params_local, c_explosive, c_saved$n)

    sem_local <- tryCatch({
      run_sem_fit(
        pp_data_in = pp_all_bg,
        partition_in = partition,
        partition_processes_in = partition_processes,
        state_spaces_in = state_spaces,
        init_params_in = biv_init_trunc,
        fixed_params_in = SENSITIVITY_FIXED_PARAMS,
        background_rate_var_in = "W",
        control_background_pre_mass_ratio_in = CTRL_BG_PRE_MASS_RATIO,
        sem_t_trunc_in = t_trunc_val,
        sem_inner_iter_in = SENS_SEM_INNER_ITER,
        verbose_in = FALSE,
        # Hold SEM RNG fixed across the trunc grid; vary only t_trunc.
        sem_rng_label_in = OK_BW_SEM_RNG_LABEL,
        label = sprintf("Fit D t_trunc %.4g", t_trunc_val)
      )
    }, error = function(e) {
      cat(sprintf("    [t_trunc_sens:%s] SEM error: %s\n", label, e$message))
      NULL
    })
    params_local <- if (valid_biv_sem(sem_local)) sem_local$etas_bivariate_params else NULL
    rho_local <- rho_of(params_local)
    # Match bootstrap/INSTRUCTIONS explosive filter (rho >= 1), and also discard
    # hard-cap landings (projected just below ETAS_BRANCHING_MAX) as unstable.
    is_explosive <- is_explosive_rho(params_local, rho_local)
    ate_local <- NULL
    if (!is.null(params_local) && !isTRUE(is_explosive)) {
      ate_local <- tryCatch({
        ate_at_trunc(params_local, sprintf("Fit D t_trunc %.4g", t_trunc_val))
      }, error = function(e) {
        cat(sprintf("    [t_trunc_sens:%s] Fit D ATE error: %s\n", label, e$message))
        NULL
      })
    } else if (isTRUE(is_explosive)) {
      cat(sprintf(
        "    [t_trunc_sens:%s] discarding explosive/near-critical Fit D (rho=%s, margin=%.3f)\n",
        label, format(rho_local, digits = 6), ETAS_BRANCHING_MAX
      ))
    }
    d_saved <- saved_stats(ate_local)
    status_local <- status_of(params_local, is_explosive, d_saved$n)
    elapsed <- proc.time()[["elapsed"]] - t0
    cat(sprintf(
      "    [t_trunc_sens:%s] done in %.1fs C=%s rhoC=%s D=%s rhoD=%s mem=%s\n",
      label, elapsed, c_status_local, format(c_rho_local, digits = 6),
      status_local, format(rho_local, digits = 6), mem_snapshot()
    ))
    list(
      t_trunc_days = as.numeric(t_trunc_val),
      status = status_local,
      C_status = c_status_local,
      D_status = status_local,
      elapsed_sec = as.numeric(elapsed),
      params = if (identical(status_local, "ok")) params_local else NULL,
      C_params = if (identical(c_status_local, "ok")) c_params_local else NULL,
      rho = as.numeric(rho_local),
      C_rho = as.numeric(c_rho_local),
      branching_metric = as.numeric(rho_local),
      n_relabel = if (!is.null(sem_local) && !is.null(sem_local$adaptive$adaptive_labelling)) {
        lp <- sem_local$adaptive$adaptive_labelling
        lp <- lp[lp$t >= 0, , drop = FALSE]
        sum(lp$location_process != lp$inferred_process, na.rm = TRUE)
      } else {
        NA_integer_
      },
      C_ate_mean_saved = c_saved$mean,
      C_ate_sd_saved = c_saved$sd,
      C_ate_n_sims = c_saved$n,
      ate_mean_saved = d_saved$mean,
      ate_sd_saved = d_saved$sd,
      ate_n_sims = d_saved$n,
      ate = if (isTRUE(TRIM_SENS_OBJECTS)) NULL else ate_local,
      C_ate = if (isTRUE(TRIM_SENS_OBJECTS)) NULL else ate_c_local,
      sem = if (isTRUE(TRIM_SENS_OBJECTS)) NULL else sem_local
    )
  }
  trunc_jobs <- as.list(T_TRUNC_SENS_DAYS)
  t_trunc_sensitivity <- if (N_CORES > 1L && length(trunc_jobs) > 1L) {
    run_parallel(
      trunc_jobs, run_t_trunc_sens_job,
      cores = min(length(trunc_jobs), max(1L, N_CORES)),
      label = "t_trunc_sens"
    )
  } else {
    lapply(trunc_jobs, run_t_trunc_sens_job)
  }
  names(t_trunc_sensitivity) <- sprintf("t_trunc_%.4g", T_TRUNC_SENS_DAYS)
  sens_ok <- sum(vapply(t_trunc_sensitivity, function(x) identical(x$status, "ok"), logical(1)))
  sens_explosive <- sum(vapply(t_trunc_sensitivity, function(x) identical(x$status, "explosive"), logical(1)))
  sens_c_ok <- sum(vapply(t_trunc_sensitivity, function(x) identical(x$C_status, "ok"), logical(1)))
  cat(sprintf(
    "  Fit C/D t_trunc sensitivity complete: C %d/%d ok | D %d/%d ok, %d D explosive/near-critical discarded\n",
    as.integer(sens_c_ok), length(t_trunc_sensitivity),
    as.integer(sens_ok), length(t_trunc_sensitivity), as.integer(sens_explosive)
  ))
}

if (isTRUE(T_TRUNC_SENS_ONLY)) {
  patch_file <- normalizePath(T_TRUNC_SENS_PATCH_FILE, winslash = "/", mustWork = TRUE)
  patched <- readRDS(patch_file)
  patched$t_trunc_sensitivity <- t_trunc_sensitivity
  if (is.null(patched$config)) patched$config <- list()
  patched$config$T_TRUNC_SENS_DAYS <- T_TRUNC_SENS_DAYS
  patched$config$RUN_T_TRUNC_SENSITIVITY <- TRUE
  patched$config$T_TRUNC_SENS_ONLY <- TRUE
  patched$config$SEM_INNER_ITER <- SEM_INNER_ITER
  patched$config$SEM_OUTER_MAXIT_BIV <- SEM_OUTER_MAXIT_BIV
  patched$metadata <- list(
    stage = "t_trunc_sensitivity_only",
    saved_at = as.character(Sys.time()),
    patched_from = patch_file,
    sem_inner_iter = SEM_INNER_ITER
  )
  out_file <- file.path(OUT_DIR, add_file_tag("oklahoma_results_t_trunc_sens.rds"))
  saveRDS(patched, out_file)
  # Also refresh the staging patch file when it lives under OUT_DIR.
  if (dirname(patch_file) == normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE)) {
    saveRDS(patched, patch_file)
  }
  cat(sprintf(
    "\n=== T_TRUNC_SENS_ONLY: wrote %s (patched from %s); exiting before bootstrap ===\n",
    out_file, patch_file
  ))
  quit(save = "no", status = 0)
}

# Save a checkpoint before bootstrap so long runs retain core fit outputs
# even if bootstrap gets interrupted or OOM-killed.
cat("\n--- Step 7b checkpoint: saving pre-bootstrap results ---\n")
# Include report-facing fields (pp_data, counties, kde_info, full config) so Quarto can
# render from this checkpoint alone — same keys as final oklahoma_results.rds.
pre_boot_saved_at <- as.character(Sys.time())

# Public A–J lettering (same as pre-sensitivity). Legacy fitE/fitF alias C/D.
fits_named <- list(
  A = list(
    letter = "A", label = sprintf("Naive bivariate (%s)", PRIMARY_PARTITION_LABEL),
    method = "county_bivariate", algorithm = "naive",
    params = B_params, fit = fitB, ate = ate_B
  ),
  B = list(
    letter = "B", label = sprintf("SEM bivariate (%s)", PRIMARY_PARTITION_LABEL),
    method = "county_bivariate", algorithm = "sem",
    params = D_params, fit = semD, ctrl = D_ctrl, treat = D_treat, ate = ate_D
  ),
  C = c(list(letter = "C", label = "Naive biv+KDE (all_free)",
             method = "kde_all_free", algorithm = "naive"),
        kde_variant_fits$E$all_free),
  D = c(list(letter = "D", label = "SEM biv+KDE (all_free)",
             method = "kde_all_free", algorithm = "sem"),
        kde_variant_fits$F$all_free),
  E = c(list(letter = "E", label = "Naive biv+KDE (productivity_free)",
             method = "kde_productivity_free", algorithm = "naive"),
        kde_variant_fits$E$productivity_free),
  F = c(list(letter = "F", label = "SEM biv+KDE (productivity_free)",
             method = "kde_productivity_free", algorithm = "sem"),
        kde_variant_fits$F$productivity_free),
  G = list(
    letter = "G", label = sprintf("Naive univariate ETAS + KDE (%s)", PRIMARY_PARTITION_LABEL),
    method = "county_univariate_kde", algorithm = "naive",
    params = I_params, fit = fitI, ate = ate_I
  ),
  H = list(
    letter = "H", label = sprintf("SEM univariate ETAS + KDE (%s)", PRIMARY_PARTITION_LABEL),
    method = "county_univariate_kde", algorithm = "sem",
    params = J_params, fit = semJ, ctrl = J_ctrl, treat = J_treat, ate = ate_J
  ),
  I = list(
    letter = "I", label = sprintf("Naive univariate ETAS (%s, no KDE)", PRIMARY_PARTITION_LABEL),
    method = "county_univariate_homog", algorithm = "naive",
    params = K_params, fit = fitK, ate = ate_K
  ),
  J = list(
    letter = "J", label = sprintf("SEM univariate ETAS (%s, no KDE)", PRIMARY_PARTITION_LABEL),
    method = "county_univariate_homog", algorithm = "sem",
    params = L_params, fit = semL, ctrl = L_ctrl, treat = L_treat, ate = ate_L
  )
)

results_pre_bootstrap <- list(
  fitB = list(params = B_params, loglik = B_loglik, fit = fitB, ate = ate_B),
  fitD = list(params = D_params, ctrl = D_ctrl, treat = D_treat, sem = semD, ate = ate_D),
  # Legacy aliases: fitE/fitF hold primary all_free content (= fits_named$C/$D).
  fitE = list(params = E_params, loglik = E_loglik, fit = fitE, ate = ate_E),
  fitF = list(params = F_params, ctrl = F_ctrl, treat = F_treat, sem = semF, ate = ate_F),
  fitI = list(params = I_params, fit = fitI, ate = ate_I),
  fitJ = list(params = J_params, fit = semJ, ctrl = J_ctrl, treat = J_treat, ate = ate_J),
  fitK = list(params = K_params, fit = fitK, ate = ate_K),
  fitL = list(params = L_params, fit = semL, ctrl = L_ctrl, treat = L_treat, ate = ate_L),
  fits_named = fits_named,
  kde_variants = kde_variant_fits,
  fitG = NULL,
  fitH = NULL,
  bootstrap_ate = NULL,
  fit_variability = fit_variability,
  partition_results = partition_results,
  ate_partitions = ate_partitions,
  kde_bandwidth_sensitivity = kde_bandwidth_sensitivity,
  t_trunc_sensitivity = t_trunc_sensitivity,
  pp_data = list(pp_pre = pp_pre, pp_pre_holdout = pp_pre_holdout, pp_post = pp_post),
  kde_info = kde_info,
  control_snapshot_fits = control_snapshot_fits,
  sem_pilot_tuning = sem_pilot_tuning,
  plots = analysis_plots,
  counties = list(
    names = counties_sf_valid$NAME,
    treated_names = county_treated_names,
    n_counties = county_partition$n,
    n_treated = sum(county_treated_idx)
  ),
  metadata = list(stage = "pre_bootstrap", saved_at = pre_boot_saved_at),
  checkpoint = list(
    stage = "pre_bootstrap",
    saved_at = pre_boot_saved_at,
    boot_targets_requested = BOOT_TARGETS,
    boot_targets_run = intersect(BOOT_TARGETS, c("C", "D")),
    boot_reps = BOOT_N_REPS
  ),
  fit_name_map = list(
    A = "fits_named$A",
    B = "fits_named$B",
    C = "fits_named$C",
    D = "fits_named$D",
    E = "fits_named$E",
    F = "fits_named$F",
    G = "fits_named$G",
    H = "fits_named$H",
    I = "fits_named$I",
    J = "fits_named$J"
  ),
  config = list(
    ETAS_M0 = ETAS_M0, BETA_GR = BETA_GR,
    BETA_GR_SOURCE = BETA_GR_SOURCE,
    GAMMA_FIXED = GAMMA_FIXED,
    FIXED_STRUCTURAL = FIXED_STRUCTURAL,
    SEM_N_ITER = SEM_N_ITER, SEM_INNER_ITER = SEM_INNER_ITER,
    SENS_SEM_INNER_ITER = SENS_SEM_INNER_ITER,
    SEM_INNER_PROPS = SEM_INNER_PROPS,
    SEM_N_LABELLINGS = SEM_N_LABELLINGS,
    SEM_OUTER_MAXIT = SEM_OUTER_MAXIT,
    SEM_OUTER_MAXIT_BIV = SEM_OUTER_MAXIT_BIV,
    SEM_WARMSTART_FIXED = SEM_WARMSTART_FIXED,
    SEM_CHANGE_FACTOR = SEM_CHANGE_FACTOR,
    SEM_CHANGE_FACTOR_MIN_MULT = SEM_CHANGE_FACTOR_MIN_MULT,
    SEM_CHANGE_FACTOR_MAX_MULT = SEM_CHANGE_FACTOR_MAX_MULT,
    SEM_OPTIM_METHOD = SEM_OPTIM_METHOD,
    SEM_MONOTONE_COMPLETE_LL = SEM_MONOTONE_COMPLETE_LL,
    SEM_START_FROM_C = SEM_START_FROM_C,
    SEM_BIV_N_THREADS = SEM_BIV_N_THREADS,
    SEM_SINGLE_FLIP_FROM_ITER = SEM_SINGLE_FLIP_FROM_ITER,
    SEM_SELECTION_TEMPERATURE = SEM_SELECTION_TEMPERATURE,
    RUN_SEM_PILOT = RUN_SEM_PILOT,
    SEM_PILOT_INNER_ITER = SEM_PILOT_INNER_ITER,
    SEM_PILOT_CORES = SEM_PILOT_CORES,
    SEM_PILOT_MAX_COMBOS = SEM_PILOT_MAX_COMBOS,
    SEM_TEMPORAL_WEIGHT = SEM_TEMPORAL_WEIGHT,
    SEM_TEMPORAL_SCALE_DAYS = SEM_TEMPORAL_SCALE_DAYS,
    SEM_T_TRUNC_DAYS = SEM_T_TRUNC_DAYS,
    SEM_T_TRUNC_SOURCE = SEM_T_TRUNC_SOURCE,
    SEM_T_TRUNC_REL = SEM_T_TRUNC_REL,
    T_TRUNC_SENS_DAYS = T_TRUNC_SENS_DAYS,
    RUN_T_TRUNC_SENSITIVITY = RUN_T_TRUNC_SENSITIVITY,
    RUN_SENSITIVITY = RUN_SENSITIVITY,
    RUN_FIT_VARIABILITY = RUN_FIT_VARIABILITY,
    FIT_VARIABILITY_REPS = FIT_VARIABILITY_REPS,
    FIT_VARIABILITY_CORES = FIT_VARIABILITY_CORES,
    KDE_VARIANT_MODE = KDE_VARIANT_MODE,
    KDE_BW_METHOD = KDE_BW_METHOD,
    PRIMARY_PARTITION = PRIMARY_PARTITION,
    KDE_BW_SENS_KM = KDE_BW_SENS_KM,
    RUN_PARTITION_SENSITIVITY = RUN_PARTITION_SENSITIVITY,
    CD_ONLY = CD_ONLY,
    UNIV_KDE_ONLY = UNIV_KDE_ONLY,
    RUN_KDE_PROFILE_SWEEP = RUN_KDE_PROFILE_SWEEP,
    MEMORY_SAFE = MEMORY_SAFE,
    TRIM_SENS_OBJECTS = TRIM_SENS_OBJECTS,
    SENS_CORES = SENS_CORES,
    SENS_BIV_N_THREADS = SENS_BIV_N_THREADS,
    ATE_SIM_CORES = ATE_SIM_CORES,
    ATE_USE_CRN = OK_ATE_USE_CRN,
    ATE_CRN_PAIR = OK_ATE_CRN_PAIR,
    ATE_CONDITIONAL_ON_PRE = OK_ATE_CONDITIONAL_ON_PRE,
    ATE_CRN_BASE = OK_ATE_CRN_BASE,
    SEM_MAX_RELABEL_STEP_FRAC = SEM_MAX_RELABEL_STEP_FRAC,
    SEM_FORCE_PARAM_UPDATE_FLIP_FRAC = SEM_FORCE_PARAM_UPDATE_FLIP_FRAC,
    ATE_N_SIMS = ATE_N_SIMS, ATE_WINDOW_DAYS = ATE_WINDOW_DAYS,
    ATE_BIVARIATE = OK_ATE_BIVARIATE, ATE_CONTRAST = OK_ATE_CONTRAST,
    RUN_BOOTSTRAP_ATE = RUN_BOOTSTRAP_ATE,
    BOOT_N_REPS = BOOT_N_REPS,
    BOOT_TARGETS = BOOT_TARGETS,
    SENS_TARGETS = SENS_TARGETS,
    BOOT_REFIT_SCOPE = BOOT_REFIT_SCOPE,
    BOOT_SEM_INNER_ITER = BOOT_SEM_INNER_ITER,
    BOOT_OUTER_CORES = BOOT_OUTER_CORES,
    BOOT_SEED = BOOT_SEED,
    BOOT_IDENTICAL_RANDOMNESS = OK_BOOT_IDENTICAL_RANDOMNESS,
    BOOT_GUARD_DEGENERATE = OK_BOOT_GUARD_DEGENERATE,
    ETAS_BRANCHING_MAX = ETAS_BRANCHING_MAX,
    ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
    ETAS_P_LOWER_BOUND = ETAS_P_LOWER_BOUND,
    ETAS_Q_LOWER_BOUND = ETAS_Q_LOWER_BOUND,
    BOOT_BRANCHING_MAX = BOOT_BRANCHING_MAX,
    BOOT_MAX_PRE_EVENTS = BOOT_MAX_PRE_EVENTS,
    BOOT_MAX_POST_EVENTS_PER_PROC = BOOT_MAX_POST_EVENTS_PER_PROC,
    BOOT_MAX_TOTAL_EVENTS = BOOT_MAX_TOTAL_EVENTS,
    TEST_MODE = TEST_MODE,
    windowT_post = windowT_post,
    n_pre = nrow(pp_pre), n_pre_holdout = nrow(pp_pre_holdout), n_pre_total = nrow(pp_pre_all),
    n_post = nrow(pp_post),
    n_counties = county_partition$n,
    n_treated_counties = sum(county_treated_idx),
    n_tiles = partition$n,
    n_treated = sum(treated_idx),
    grid_diameters = grid_diameters,
    grid_multipliers = grid_multipliers,
    kde_bandwidth_multipliers = vapply(kde_bandwidth_specs, `[[`, numeric(1), "multiplier"),
    spatial_scale_D = D_scale
  )
)
pre_bootstrap_out_file <- file.path(OUT_DIR, add_file_tag("oklahoma_results_pre_bootstrap.rds"))
saveRDS(results_pre_bootstrap, pre_bootstrap_out_file)
cat(sprintf("Pre-bootstrap checkpoint saved to: %s\n", pre_bootstrap_out_file))
rm(results_pre_bootstrap)
invisible(gc(verbose = FALSE))
} # !BOOTSTRAP_ONLY: skip sensitivity / pre-bootstrap checkpoints

# Parametric bootstrap ATEs for the all-free KDE pair (C/D).
bootstrap_ate <- NULL
boot_targets_run <- intersect(BOOT_TARGETS, c("C", "D"))
if ("C" %in% boot_targets_run && is.null(E_params)) {
  cat("  Dropping bootstrap target C because fitted parameters are NULL.\n")
  boot_targets_run <- setdiff(boot_targets_run, "C")
}
if ("D" %in% boot_targets_run && is.null(F_params)) {
  cat("  Dropping bootstrap target D because fitted parameters are NULL.\n")
  boot_targets_run <- setdiff(boot_targets_run, "D")
}
bootstrap_elapsed <- NA_real_
if (RUN_BOOTSTRAP_ATE && BOOT_N_REPS > 0L && length(boot_targets_run) > 0L) {
  t_bootstrap <- proc.time()[["elapsed"]]
  cat(sprintf("\n--- Step 8: Parametric bootstrap ATEs (targets=%s, reps=%d, scope=%s, outer cores=%d, boot SEM inner=%d) ---\n",
              paste(boot_targets_run, collapse = ","), BOOT_N_REPS, BOOT_REFIT_SCOPE, BOOT_OUTER_CORES, BOOT_SEM_INNER_ITER))

  # R's default connection limit is 128. The reusable ATE PSOCK pool already
  # holds ATE_SIM_CORES sockets; opening BOOT_OUTER_CORES more workers on top
  # of that hits "all connections are in use" (job 8469262: 64+64).
  stop_ate_psock_pool()
  cat("  Closed reusable ATE PSOCK pool before bootstrap outer workers.\n")

  if (BOOT_REFIT_SCOPE == "full") {
    cat("  Full scope selected: for current targets (C/D), this runs per-replicate refits before ATE estimation.\n")
  }

  boot_data_cols <- c(
    "x", "y", "t", "mag", "location_process", "process",
    "inferred_process", "W", "n", "background"
  )
  canonicalize_boot_df <- function(df) {
    df <- as.data.frame(df)
    n <- nrow(df)
    defaults <- list(
      x = numeric(n), y = numeric(n), t = numeric(n), mag = numeric(n),
      location_process = rep(NA_character_, n),
      process = rep(NA_character_, n),
      inferred_process = rep(NA_character_, n),
      W = rep(1, n), n = rep.int(n, n), background = rep(TRUE, n)
    )
    for (nm in boot_data_cols) {
      if (is.null(df[[nm]])) df[[nm]] <- defaults[[nm]]
    }
    for (nm in c("location_process", "process", "inferred_process")) {
      df[[nm]] <- as.character(df[[nm]])
    }
    for (nm in c("x", "y", "t", "mag", "W")) {
      df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
    }
    df$n <- as.integer(df$n)
    df$background <- as.logical(df$background)
    df[, boot_data_cols, drop = FALSE]
  }
  as_pp_df <- function(sim_obj, location_process_value, inferred_process_value = NULL) {
    if (is.null(sim_obj) || is.null(sim_obj$t) || length(sim_obj$t) < 1) {
      return(canonicalize_boot_df(data.frame()))
    }
    out <- data.frame(
      x = sim_obj$x, y = sim_obj$y, t = sim_obj$t, mag = sim_obj$mag,
      stringsAsFactors = FALSE
    )
    out$location_process <- location_process_value
    out$process <- location_process_value
    out$inferred_process <- if (is.null(inferred_process_value)) location_process_value else inferred_process_value
    out$W <- if (!is.null(sim_obj$W) && length(sim_obj$W) == nrow(out)) sim_obj$W else 1
    out$background <- if (!is.null(sim_obj$background) &&
                          length(sim_obj$background) == nrow(out)) {
      sim_obj$background
    } else {
      FALSE
    }
    canonicalize_boot_df(out)
  }
  as_pp_bivariate_df <- function(sim_obj) {
    if (is.null(sim_obj) || is.null(sim_obj$t) || length(sim_obj$t) < 1L) {
      return(as_pp_df(NULL, "control", "control"))
    }
    out <- data.frame(
      x = sim_obj$x, y = sim_obj$y, t = sim_obj$t, mag = sim_obj$mag,
      stringsAsFactors = FALSE
    )
    in_treated <- spatstat.geom::inside.owin(
      out$x, out$y, w = treated_ss
    )
    out$location_process <- ifelse(in_treated, "treated", "control")
    out$process <- ifelse(as.integer(sim_obj$process_id) == 1L, "treated", "control")
    # Naive labels are geographic; SEM starts from the same observed labels.
    out$inferred_process <- out$location_process
    out$background <- if (!is.null(sim_obj$background) &&
                          length(sim_obj$background) == nrow(out)) {
      sim_obj$background
    } else {
      FALSE
    }
    canonicalize_boot_df(out)
  }

  post_window_boot <- windowT_post
  # Bootstrap targets C/D correspond to the all-free KDE pair (internal E_*/F_*).
  c_ctrl_seed <- E_marginals$ctrl
  c_treat_seed <- E_marginals$treat
  d_ctrl_seed <- F_marginals$ctrl
  d_treat_seed <- F_marginals$treat

  # Same normalization as the global ATE lookups (mean one on each observed
  # region, evaluable everywhere).
  boot_covariate_lookup <- KDE_BG_LOOKUP

  get_num <- function(obj, nm, default = NA_real_) {
    if (is.null(obj) || is.null(obj[[nm]])) return(default)
    v <- suppressWarnings(as.numeric(obj[[nm]]))
    if (length(v) < 1L || !is.finite(v[[1]])) return(default)
    v[[1]]
  }
  univ_branching_ratio <- function(par_obj, beta_gr = BETA_GR) {
    A <- get_num(par_obj, "A", NA_real_)
    alpha_m <- get_num(par_obj, "alpha_m", NA_real_)
    if (!is.finite(A) || A < 0) return(Inf)
    gap <- beta_gr - alpha_m
    if (!is.finite(gap) || gap <= 1e-8) return(Inf)
    A * beta_gr / gap
  }
  biv_branching_radius <- function(par_obj, beta_gr = BETA_GR) {
    eta_channel <- function(A_name, alpha_name) {
      A <- get_num(par_obj, A_name, NA_real_)
      alpha <- get_num(par_obj, alpha_name, NA_real_)
      gap <- beta_gr - alpha
      if (!is.finite(A) || A < 0 || !is.finite(gap) || gap <= 1e-8) {
        return(Inf)
      }
      A * beta_gr / gap
    }
    offspring <- matrix(
      c(
        eta_channel("A_00", "alpha_m_00"),
        eta_channel("A_01", "alpha_m_01"),
        eta_channel("A_10", "alpha_m_10"),
        eta_channel("A_11", "alpha_m_11")
      ),
      nrow = 2, byrow = TRUE
    )
    if (!all(is.finite(offspring))) return(Inf)
    max(Re(eigen(offspring, only.values = TRUE)$values))
  }
  assess_boot_stability <- function(par_obj, bivariate = OK_ATE_BIVARIATE) {
    if (isTRUE(bivariate)) {
      rho <- biv_branching_radius(par_obj)
      return(list(
        law = "bivariate",
        branching_metric = rho,
        rho = rho,
        eta_ctrl = NA_real_,
        eta_treat = NA_real_,
        stable = is.finite(rho) && rho < 1,
        within_fit_margin = is.finite(rho) && rho < BOOT_BRANCHING_MAX
      ))
    }
    ctrl <- par_obj$control
    treat <- par_obj$treated
    eta_ctrl <- univ_branching_ratio(ctrl)
    eta_treat <- univ_branching_ratio(treat)
    metric <- max(c(eta_ctrl, eta_treat))
    list(
      law = "univariate",
      branching_metric = metric,
      rho = NA_real_,
      eta_ctrl = eta_ctrl,
      eta_treat = eta_treat,
      stable = all(is.finite(c(eta_ctrl, eta_treat))) && metric < 1,
      within_fit_margin = all(is.finite(c(eta_ctrl, eta_treat))) &&
        metric < BOOT_BRANCHING_MAX
    )
  }
  validate_sim_obj <- function(sim_obj, label, max_events, support = win_km) {
    n_ev <- if (is.null(sim_obj) || is.null(sim_obj$t)) 0L else length(sim_obj$t)
    if (n_ev > max_events) {
      stop(sprintf("%s produced %d events (> %d cap)", label, n_ev, max_events))
    }
    if (n_ev > 0L) {
      x_ok <- !any(!is.finite(sim_obj$x))
      y_ok <- !any(!is.finite(sim_obj$y))
      t_ok <- !any(!is.finite(sim_obj$t))
      m_ok <- !any(!is.finite(sim_obj$mag))
      if (!(x_ok && y_ok && t_ok && m_ok)) {
        stop(sprintf("%s produced non-finite values", label))
      }
      inside <- spatstat.geom::inside.owin(sim_obj$x, sim_obj$y, support)
      if (!all(inside)) {
        stop(sprintf(
          "%s produced %d event(s) outside modeled support",
          label, sum(!inside)
        ))
      }
    }
    invisible(TRUE)
  }

  simulate_boot_data <- function(ctrl_seed,
                                 treat_seed,
                                 biv_seed = NULL,
                                 bivariate = OK_ATE_BIVARIATE) {
    if (OK_BOOT_GUARD_DEGENERATE) {
      dgp_stability <- if (isTRUE(bivariate)) {
        assess_boot_stability(biv_seed, bivariate = TRUE)
      } else {
        assess_boot_stability(
          list(control = ctrl_seed, treated = treat_seed),
          bivariate = FALSE
        )
      }
      if (!isTRUE(dgp_stability$within_fit_margin)) {
        stop(sprintf(
          paste0(
            "Bootstrap DGP is outside the %.3f stability margin ",
            "(%s metric=%.3f)"
          ),
          BOOT_BRANCHING_MAX, dgp_stability$law,
          dgp_stability$branching_metric
        ))
      }
    }
    # Conditional bootstrap: hold the retained observed pre-treatment history
    # fixed and simulate only the post-treatment process from each target law.
    pre_df <- canonicalize_boot_df(pp_pre)
    pre_df$process <- "control"
    pre_df$inferred_process <- "control"
    validate_sim_obj(
      pre_df, "bootstrap conditioning history", BOOT_MAX_PRE_EVENTS,
      support = win_km
    )
    history_df <- pre_df[, c("x", "y", "t", "mag"), drop = FALSE]
    if (isTRUE(bivariate)) {
      if (is.null(biv_seed)) stop("Bivariate bootstrap requires biv_seed.")
      # The post-treatment bootstrap is conditional on the complete retained
      # observed pre-treatment history. Pre-treatment belongs to the control process,
      # irrespective of the county that later receives treatment.
      history_biv <- pre_df[, c("x", "y", "t", "mag"), drop = FALSE]
      history_biv$process_id <- 0L
      if (nrow(history_biv) > 0L &&
          any(history_biv$t >= post_window_boot[1], na.rm = TRUE)) {
        stop("Bootstrap conditioning history must be strictly pre-treatment.")
      }
      joint_post_sim <- sim_etas_bivariate(
        params = biv_seed,
        windowT = post_window_boot,
        windowS = win_km,
        state_spaces = state_spaces,
        m0 = ETAS_M0,
        beta_gr = BETA_GR,
        filtration = history_biv,
        covariate_lookup = boot_covariate_lookup,
        t_trunc = SEM_T_TRUNC_DAYS
      )
      validate_sim_obj(
        joint_post_sim,
        "bootstrap bivariate post_sim",
        2L * BOOT_MAX_POST_EVENTS_PER_PROC,
        support = win_km
      )
      joint_post_df <- as_pp_bivariate_df(joint_post_sim)
      ctrl_post_sim <- joint_post_df[
        joint_post_df$location_process == "control", , drop = FALSE
      ]
      treat_post_sim <- joint_post_df[
        joint_post_df$location_process == "treated", , drop = FALSE
      ]
    } else {
      treated_history_df <- history_df[0, , drop = FALSE]
      ctrl_post_sim <- sim_etas(
        ctrl_seed, post_window_boot, windowS = control_ss,
        m0 = ETAS_M0, beta_gr = BETA_GR, filtration = history_df,
        covariate_lookup = boot_covariate_lookup$control,
        t_trunc = SEM_T_TRUNC_DAYS
      )
      treat_post_sim <- sim_etas(
        treat_seed, post_window_boot, windowS = treated_ss,
        m0 = ETAS_M0, beta_gr = BETA_GR, filtration = treated_history_df,
        covariate_lookup = boot_covariate_lookup$treated,
        t_trunc = SEM_T_TRUNC_DAYS
      )
    }
    validate_sim_obj(
      ctrl_post_sim, "bootstrap ctrl_post_sim",
      BOOT_MAX_POST_EVENTS_PER_PROC,
      support = control_ss
    )
    validate_sim_obj(
      treat_post_sim, "bootstrap treat_post_sim",
      BOOT_MAX_POST_EVENTS_PER_PROC,
      support = treated_ss
    )
    if (isTRUE(bivariate)) {
      post_ctrl_df <- ctrl_post_sim
      post_treat_df <- treat_post_sim
    } else {
      post_ctrl_df <- as_pp_df(ctrl_post_sim, "control", "control")
      post_treat_df <- as_pp_df(treat_post_sim, "treated", "treated")
    }
    pp_post_sim <- rbind(post_ctrl_df, post_treat_df)
    pp_post_sim <- pp_post_sim[order(pp_post_sim$t), , drop = FALSE]

    post_ctrl_df <- canonicalize_boot_df(post_ctrl_df)
    post_treat_df <- canonicalize_boot_df(post_treat_df)
    bg_ctrl <- canonicalize_boot_df(
      normalize_bg_weights(post_ctrl_df, control_ss, lambda_im)$new_df
    )
    bg_treat <- canonicalize_boot_df(
      normalize_bg_weights(post_treat_df, treated_ss, lambda_im)$new_df
    )
    pp_post_bg_sim <- rbind(bg_ctrl, bg_treat)
    pp_post_bg_sim <- pp_post_bg_sim[order(pp_post_bg_sim$t), , drop = FALSE]

    pre_bg_ctrl <- canonicalize_boot_df(
      normalize_bg_weights(pre_df, control_ss, lambda_im)$new_df
    )
    pp_all_bg_sim <- rbind(pre_bg_ctrl, pp_post_bg_sim)
    pp_all_bg_sim$n <- nrow(pp_all_bg_sim)
    pp_all_bg_sim <- pp_all_bg_sim[order(pp_all_bg_sim$t), , drop = FALSE]
    if (nrow(pp_all_bg_sim) > BOOT_MAX_TOTAL_EVENTS) {
      stop(sprintf("bootstrap total simulated events %d exceeded cap %d", nrow(pp_all_bg_sim), BOOT_MAX_TOTAL_EVENTS))
    }

    list(
      pre_df = pre_df,
      post_ctrl_df = post_ctrl_df,
      post_treat_df = post_treat_df,
      pp_post_bg_sim = pp_post_bg_sim,
      pp_all_bg_sim = pp_all_bg_sim,
      simulation_law = if (isTRUE(bivariate)) {
        "bivariate_post_conditional_on_observed_pre"
      } else {
        "univariate_post_conditional_on_observed_pre"
      },
      post_conditioned_on_pre = TRUE,
      conditioning_history = "observed_retained_pre_treatment",
      n_conditioning_events = nrow(pre_df)
    )
  }

  stop_explosive_boot <- function(stability, target, rep_id) {
    if (isTRUE(stability$stable)) return(invisible(stability))
    cond <- structure(
      list(
        message = sprintf(
          "Explosive %s bootstrap refit at rep %d (%s metric=%s)",
          target, as.integer(rep_id), stability$law,
          format(stability$branching_metric, digits = 6)
        ),
        call = NULL,
        stability = stability,
        target = target,
        rep = rep_id
      ),
      class = c("oklahoma_boot_explosive", "error", "condition")
    )
    stop(cond)
  }
  boot_failure <- function(e, target, rep_id) {
    is_explosive <- inherits(e, "oklahoma_boot_explosive")
    list(
      ok = FALSE,
      rep = as.integer(rep_id),
      msg = paste0("Bootstrap ", target, " failed: ", conditionMessage(e)),
      failure_type = if (is_explosive) "explosive" else "error",
      stability = if (is_explosive) e$stability else NULL
    )
  }
  summarize_boot <- function(ate_obj,
                             rep_id,
                             pre_df,
                             post_ctrl_df,
                             post_treat_df,
                             par_obj = NULL,
                             stability,
                             simulation_law) {
    if (is.null(ate_obj) || is.null(ate_obj$all_nothing_sim)) {
      return(list(ok = FALSE, rep = rep_id, msg = "ATE estimation failed"))
    }
    total_saved <- suppressWarnings(as.numeric(
      ate_obj$all_nothing_sim$total_saved
    ))
    tile_ate <- suppressWarnings(as.numeric(ate_obj$all_nothing_sim$ATE))
    if (length(total_saved) < 1L || !all(is.finite(total_saved)) ||
        length(tile_ate) < 1L || !all(is.finite(tile_ate))) {
      return(list(
        ok = FALSE, rep = as.integer(rep_id)[1L],
        msg = "ATE estimation returned empty or non-finite values",
        failure_type = "nonfinite_ate"
      ))
    }
    list(
      ok = TRUE,
      rep = as.integer(rep_id)[1L],
      n_pre_sim = as.integer(nrow(pre_df))[1L],
      n_post_ctrl_sim = as.integer(nrow(post_ctrl_df))[1L],
      n_post_treat_sim = as.integer(nrow(post_treat_df))[1L],
      ate_total_mean = as.numeric(mean(total_saved))[1L],
      ate_total_sd = as.numeric(stats::sd(total_saved))[1L],
      ate_tile_mean = as.numeric(mean(tile_ate))[1L],
      fit_law = stability$law,
      simulation_law = simulation_law,
      branching_metric = as.numeric(stability$branching_metric),
      rho = as.numeric(stability$rho),
      eta_ctrl = as.numeric(stability$eta_ctrl),
      eta_treat = as.numeric(stability$eta_treat),
      is_stable = isTRUE(stability$stable),
      within_fit_margin = isTRUE(stability$within_fit_margin),
      params = par_obj
    )
  }

  run_boot_rep <- function(rep_id) {
    t0_rep <- proc.time()[["elapsed"]]
    cat(sprintf("    [bootstrap:rep-%d] start pid=%d mem=%s\n",
                as.integer(rep_id), Sys.getpid(), mem_snapshot()))
    seed_base <- if (!is.na(BOOT_SEED)) BOOT_SEED else if (is.finite(OK_GLOBAL_SEED) && !is.na(OK_GLOBAL_SEED)) OK_GLOBAL_SEED else NA_integer_
    if (isTRUE(OK_BOOT_IDENTICAL_RANDOMNESS)) {
      boot_rep_seed <- if (!is.na(seed_base)) as.integer(seed_base) else 1L
    } else {
      boot_rep_seed <- if (!is.na(seed_base)) as.integer(seed_base + rep_id) else as.integer(10000L + rep_id)
    }
    set.seed(boot_rep_seed)
    boot_ate_crn_seed <- as.integer(100000L + 1000L * boot_rep_seed)
    out <- list(rep = rep_id)

    if ("C" %in% boot_targets_run) {
      out$C <- tryCatch({
        sim_C <- simulate_boot_data(
          c_ctrl_seed, c_treat_seed,
          biv_seed = E_params,
          bivariate = OK_ATE_BIVARIATE
        )
        c_params_boot <- if (isTRUE(OK_ATE_BIVARIATE)) {
          E_params
        } else {
          list(control = c_ctrl_seed, treated = c_treat_seed)
        }
        if (BOOT_REFIT_SCOPE %in% c("partial", "full")) {
          if (isTRUE(OK_ATE_BIVARIATE)) {
            fit_c_boot <- fit_etas_bivariate(
              params_init = E_params, realiz = sim_C$pp_all_bg_sim,
              windowT = windowT_fit, windowS = win_km, m0 = ETAS_M0,
              control_state_space = control_ss, treated_state_space = treated_ss,
              background_rate_var = "W",
              treated_background_zero_before = 0,
              control_background_everywhere_before = 0,
              control_background_pre_mass_ratio = CTRL_BG_PRE_MASS_RATIO,
              beta_gr = BETA_GR,
              enforce_finite_trigger_moments = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
              p_lower_bound = ETAS_P_LOWER_BOUND,
              q_lower_bound = ETAS_Q_LOWER_BOUND,
              max_branching_radius = ETAS_BRANCHING_MAX,
              maxit = VANILLA_MAXIT, fixed_params = SENSITIVITY_FIXED_PARAMS, trace = 0,
              t_trunc = SEM_T_TRUNC_DAYS
            )
            if (!valid_biv_fit(fit_c_boot)) {
              stop("Bivariate naive refit did not converge to a usable fit.")
            }
            c_params_boot <- fit_c_boot$par
          } else {
            fit_c_boot <- fit_indep_pair(
              pp_data_in = sim_C$pp_all_bg_sim,
              background_rate_var_in = "W",
              fixed_params_in = SENSITIVITY_FIXED_PARAMS,
              ctrl_init_in = c_ctrl_seed,
              treat_init_in = c_treat_seed,
              fit_label = sprintf("Boot C #%d", rep_id)
            )
            if (is.null(fit_c_boot) || is.null(fit_c_boot$fit_control) ||
                is.null(fit_c_boot$fit_treated)) {
              stop("Univariate naive refit did not complete both process fits.")
            }
            c_params_boot <- list(
              control = fit_c_boot$control,
              treated = fit_c_boot$treated
            )
          }
        }
        c_stability <- assess_boot_stability(
          c_params_boot, bivariate = OK_ATE_BIVARIATE
        )
        stop_explosive_boot(c_stability, "C", rep_id)
        c_marg_boot <- if (isTRUE(OK_ATE_BIVARIATE)) {
          extract_marginals(c_params_boot)
        } else {
          list(ctrl = c_params_boot$control, treat = c_params_boot$treated)
        }
        ate_c_boot <- if (isTRUE(OK_ATE_BIVARIATE)) {
          ate_estim_bivariate(
            biv_params = c_params_boot,
            windowT = windowT_ate,
            windowS = win_km,
            state_spaces_obs = state_spaces,
            label = sprintf("Boot C #%d", rep_id),
            n_sims = ATE_N_SIMS,
            # Nested under bootstrap outer parallel: keep ATE sims sequential.
            n_cores = 1L,
            m0 = ETAS_M0,
            beta_gr = BETA_GR,
            filtration_history = if (isTRUE(OK_ATE_CONDITIONAL_ON_PRE)) sim_C$pre_df else NULL,
            t_trunc = SEM_T_TRUNC_DAYS,
            n_tiles = partition$n,
            crn_base_seed = boot_ate_crn_seed,
            use_crn = OK_ATE_USE_CRN,
            crn_pair = OK_ATE_CRN_PAIR,
            quiet = TRUE,
            contrast = OK_ATE_CONTRAST,
            covariate_lookup = boot_covariate_lookup
          )
        } else {
          ate_estim_fast(
            c_marg_boot$ctrl, c_marg_boot$treat, sim_C$pp_post_bg_sim,
            label = sprintf("Boot C #%d", rep_id),
            filtration_history = if (isTRUE(OK_ATE_CONDITIONAL_ON_PRE)) sim_C$pre_df else NULL,
            crn_base_seed = boot_ate_crn_seed,
            phase = "bootstrap",
            n_tiles_used = partition$n,
            treated_idx_used = treated_idx,
            quiet = TRUE,
            contrast = OK_ATE_CONTRAST,
            covariate_lookup = boot_covariate_lookup,
            t_trunc = SEM_T_TRUNC_DAYS
          )
        }
        summarize_boot(
          ate_c_boot, rep_id,
          sim_C$pre_df, sim_C$post_ctrl_df, sim_C$post_treat_df,
          par_obj = c_params_boot,
          stability = c_stability,
          simulation_law = sim_C$simulation_law
        )
      }, error = function(e) {
        boot_failure(e, "C", rep_id)
      })
    }

    if ("D" %in% boot_targets_run) {
      out$D <- tryCatch({
        sim_D <- simulate_boot_data(
          d_ctrl_seed, d_treat_seed,
          biv_seed = F_params,
          bivariate = OK_ATE_BIVARIATE
        )
        d_params_boot <- if (isTRUE(OK_ATE_BIVARIATE)) {
          F_params
        } else {
          list(control = d_ctrl_seed, treated = d_treat_seed)
        }
        pp_post_d_boot <- sim_D$pp_post_bg_sim
        if (BOOT_REFIT_SCOPE %in% c("partial", "full")) {
          sem_boot <- run_sem_fit(
            pp_data_in = sim_D$pp_all_bg_sim,
            partition_in = partition,
            partition_processes_in = partition_processes,
            state_spaces_in = state_spaces,
            init_params_in = if (isTRUE(OK_ATE_BIVARIATE)) F_params else NULL,
            init_ctrl_params_in = d_ctrl_seed,
            init_treat_params_in = d_treat_seed,
            model_type_in = if (isTRUE(OK_ATE_BIVARIATE)) "etas_bivariate" else "etas",
            fixed_params_in = SENSITIVITY_FIXED_PARAMS,
            background_rate_var_in = "W",
            control_background_pre_mass_ratio_in = CTRL_BG_PRE_MASS_RATIO,
            sem_inner_iter_in = BOOT_SEM_INNER_ITER,
            verbose_in = FALSE,
            label = sprintf("Boot D #%d", rep_id)
          )
          if (is.null(sem_boot)) stop("SEM refit returned NULL.")
          if (isTRUE(OK_ATE_BIVARIATE)) {
            if (!valid_biv_sem(sem_boot)) {
              stop("Bivariate SEM refit did not converge to a usable fit.")
            }
            d_params_boot <- sem_boot$etas_bivariate_params
          } else {
            if (is.null(sem_boot$hawkes_params_control) ||
                is.null(sem_boot$hawkes_params_treated)) {
              stop("Univariate SEM refit returned incomplete process parameters.")
            }
            d_params_boot <- list(
              control = sem_boot$hawkes_params_control,
              treated = sem_boot$hawkes_params_treated
            )
          }
          if (!is.null(sem_boot$adaptive$adaptive_labelling)) {
            pp_post_d_boot <- sem_boot$adaptive$adaptive_labelling
            pp_post_d_boot <- pp_post_d_boot[pp_post_d_boot$t >= 0, , drop = FALSE]
          }
        }
        d_stability <- assess_boot_stability(
          d_params_boot, bivariate = OK_ATE_BIVARIATE
        )
        stop_explosive_boot(d_stability, "D", rep_id)
        d_marg_boot <- if (isTRUE(OK_ATE_BIVARIATE)) {
          extract_marginals(d_params_boot)
        } else {
          list(ctrl = d_params_boot$control, treat = d_params_boot$treated)
        }
        ate_d_boot <- if (isTRUE(OK_ATE_BIVARIATE)) {
          ate_estim_bivariate(
            biv_params = d_params_boot,
            windowT = windowT_ate,
            windowS = win_km,
            state_spaces_obs = state_spaces,
            label = sprintf("Boot D #%d", rep_id),
            n_sims = ATE_N_SIMS,
            # Nested under bootstrap outer parallel: keep ATE sims sequential.
            n_cores = 1L,
            m0 = ETAS_M0,
            beta_gr = BETA_GR,
            filtration_history = if (isTRUE(OK_ATE_CONDITIONAL_ON_PRE)) sim_D$pre_df else NULL,
            t_trunc = SEM_T_TRUNC_DAYS,
            n_tiles = partition$n,
            crn_base_seed = boot_ate_crn_seed,
            use_crn = OK_ATE_USE_CRN,
            crn_pair = OK_ATE_CRN_PAIR,
            quiet = TRUE,
            contrast = OK_ATE_CONTRAST,
            covariate_lookup = boot_covariate_lookup
          )
        } else {
          ate_estim_fast(
            d_marg_boot$ctrl, d_marg_boot$treat, pp_post_d_boot,
            label = sprintf("Boot D #%d", rep_id),
            filtration_history = if (isTRUE(OK_ATE_CONDITIONAL_ON_PRE)) sim_D$pre_df else NULL,
            crn_base_seed = boot_ate_crn_seed,
            phase = "bootstrap",
            n_tiles_used = partition$n,
            treated_idx_used = treated_idx,
            quiet = TRUE,
            contrast = OK_ATE_CONTRAST,
            covariate_lookup = boot_covariate_lookup,
            t_trunc = SEM_T_TRUNC_DAYS
          )
        }
        summarize_boot(
          ate_d_boot, rep_id,
          sim_D$pre_df, sim_D$post_ctrl_df, sim_D$post_treat_df,
          par_obj = d_params_boot,
          stability = d_stability,
          simulation_law = sim_D$simulation_law
        )
      }, error = function(e) {
        boot_failure(e, "D", rep_id)
      })
    }
    rm_vars <- intersect(
      c("sim_C", "fit_c_boot", "c_params_boot", "c_marg_boot", "ate_c_boot",
        "sim_D", "sem_boot", "d_params_boot", "d_marg_boot", "ate_d_boot", "pp_post_d_boot"),
      ls()
    )
    if (length(rm_vars) > 0L) rm(list = rm_vars)
    invisible(gc(verbose = FALSE))
    elapsed_rep <- proc.time()[["elapsed"]] - t0_rep
    ok_targets <- names(out)[names(out) %in% c("C", "D")]
    ok_count <- if (length(ok_targets) > 0) {
      sum(vapply(out[ok_targets], function(x) !is.null(x) && isTRUE(x$ok), logical(1)))
    } else {
      0L
    }
    cat(sprintf("    [bootstrap:rep-%d] done in %.1fs successful_targets=%d/%d mem=%s\n",
                as.integer(rep_id), elapsed_rep, as.integer(ok_count), as.integer(length(ok_targets)),
                mem_snapshot()))
    out
  }

  boot_results <- if (BOOT_OUTER_CORES > 1L && BOOT_N_REPS > 1L) {
    run_parallel(
      as.list(seq_len(BOOT_N_REPS)), run_boot_rep,
      cores = min(BOOT_OUTER_CORES, BOOT_N_REPS),
      label = "bootstrap"
    )
  } else {
    lapply(seq_len(BOOT_N_REPS), run_boot_rep)
  }

  make_boot_block <- function(key) {
    n_attempted <- length(boot_results)
    model_res <- lapply(boot_results, function(z) {
      if (is.null(z) || is.null(z[[key]])) return(NULL)
      z[[key]]
    })
    model_res <- Filter(Negate(is.null), model_res)
    n_missing <- n_attempted - length(model_res)
    ok <- Filter(function(x) isTRUE(x$ok), model_res)
    fail <- Filter(function(x) !isTRUE(x$ok), model_res)
    scalar1 <- function(x, default = NA_real_) {
      if (is.null(x) || length(x) < 1L) return(default)
      out <- x[[1L]]
      if (length(out) != 1L) out <- default
      out
    }
    rep_df <- if (length(ok) > 0) {
      do.call(rbind, lapply(ok, function(z) {
        data.frame(
          rep = as.integer(scalar1(z$rep, NA_integer_)),
          n_pre_sim = as.integer(scalar1(z$n_pre_sim, NA_integer_)),
          n_post_ctrl_sim = as.integer(scalar1(z$n_post_ctrl_sim, NA_integer_)),
          n_post_treat_sim = as.integer(scalar1(z$n_post_treat_sim, NA_integer_)),
          ate_total_mean = as.numeric(scalar1(z$ate_total_mean, NA_real_)),
          ate_total_sd = as.numeric(scalar1(z$ate_total_sd, NA_real_)),
          ate_tile_mean = as.numeric(scalar1(z$ate_tile_mean, NA_real_)),
          fit_law = as.character(scalar1(z$fit_law, NA_character_)),
          simulation_law = as.character(scalar1(z$simulation_law, NA_character_)),
          branching_metric = as.numeric(scalar1(z$branching_metric, NA_real_)),
          rho = as.numeric(scalar1(z$rho, NA_real_)),
          eta_ctrl = as.numeric(scalar1(z$eta_ctrl, NA_real_)),
          eta_treat = as.numeric(scalar1(z$eta_treat, NA_real_)),
          is_stable = isTRUE(scalar1(z$is_stable, FALSE)),
          within_fit_margin = isTRUE(scalar1(z$within_fit_margin, FALSE)),
          stringsAsFactors = FALSE
        )
      }))
    } else {
      data.frame(
        rep = integer(),
        n_pre_sim = integer(),
        n_post_ctrl_sim = integer(),
        n_post_treat_sim = integer(),
        ate_total_mean = numeric(),
        ate_total_sd = numeric(),
        ate_tile_mean = numeric(),
        fit_law = character(),
        simulation_law = character(),
        branching_metric = numeric(),
        rho = numeric(),
        eta_ctrl = numeric(),
        eta_treat = numeric(),
        is_stable = logical(),
        within_fit_margin = logical(),
        stringsAsFactors = FALSE
      )
    }
    failure_types <- vapply(
      fail,
      function(z) as.character(scalar1(z$failure_type, "error")),
      character(1)
    )
    list(
      replicate_summary = rep_df,
      params = lapply(ok, function(z) z$params),
      n_attempted = n_attempted,
      n_retained = length(ok),
      n_success = length(ok),
      n_fail = length(fail) + n_missing,
      n_explosive = sum(failure_types == "explosive"),
      n_error = sum(failure_types != "explosive") + n_missing,
      n_missing = n_missing,
      failures = fail
    )
  }

  bootstrap_ate <- list(
    config = list(
      reps = BOOT_N_REPS,
      targets = boot_targets_run,
      targets_requested = BOOT_TARGETS,
      targets_run = boot_targets_run,
      refit_scope = BOOT_REFIT_SCOPE,
      outer_cores = BOOT_OUTER_CORES,
      sem_inner_iter = BOOT_SEM_INNER_ITER,
      seed = BOOT_SEED,
      simulation_law = if (isTRUE(OK_ATE_BIVARIATE)) {
        "bivariate_post_conditional_on_observed_pre"
      } else {
        "univariate_post_conditional_on_observed_pre"
      },
      conditioning_history = "observed_retained_pre_treatment",
      n_conditioning_events = nrow(pp_pre),
      refit_law = if (isTRUE(OK_ATE_BIVARIATE)) "bivariate" else "univariate",
      stability_fit_margin = BOOT_BRANCHING_MAX,
      explosive_drop_threshold = 1,
      bias_correction_stage = "after_stability_filter"
    )
  )
  if ("C" %in% boot_targets_run) bootstrap_ate$fit_C <- make_boot_block("C")
  if ("D" %in% boot_targets_run) bootstrap_ate$fit_D <- make_boot_block("D")
  # Legacy aliases: older report/scenario code looked for fit_E/fit_F as primary.
  if (!is.null(bootstrap_ate$fit_C) && is.null(bootstrap_ate$fit_E)) bootstrap_ate$fit_E <- bootstrap_ate$fit_C
  if (!is.null(bootstrap_ate$fit_D) && is.null(bootstrap_ate$fit_F)) bootstrap_ate$fit_F <- bootstrap_ate$fit_D
  zero_retained <- boot_targets_run[vapply(
    boot_targets_run,
    function(target) {
      block <- bootstrap_ate[[paste0("fit_", target)]]
      is.null(block) || block$n_retained < 1L
    },
    logical(1)
  )]
  if (length(zero_retained) > 0L) {
    stop(
      "Bootstrap retained zero successful replicates for requested target(s): ",
      paste(zero_retained, collapse = ", "),
      ". Publication artifacts were not written."
    )
  }
rm(boot_results)
invisible(gc(verbose = FALSE))

  if (!is.null(bootstrap_ate$fit_C)) {
    bC <- bootstrap_ate$fit_C$replicate_summary
    if (nrow(bC) > 0) {
      cat(sprintf("  Bootstrap C complete: retained=%d/%d, explosive=%d, other_fail=%d, mean total ATE=%.1f, SD(rep means)=%.1f\n",
                  nrow(bC), bootstrap_ate$fit_C$n_attempted,
                  bootstrap_ate$fit_C$n_explosive,
                  bootstrap_ate$fit_C$n_error,
                  mean(bC$ate_total_mean, na.rm = TRUE),
                  stats::sd(bC$ate_total_mean, na.rm = TRUE)))
    } else {
      cat(sprintf("  Bootstrap C complete: success=%d, fail=%d\n",
                  bootstrap_ate$fit_C$n_success, bootstrap_ate$fit_C$n_fail))
    }
  }
  if (!is.null(bootstrap_ate$fit_D)) {
    bD <- bootstrap_ate$fit_D$replicate_summary
    if (nrow(bD) > 0) {
      cat(sprintf("  Bootstrap D complete: retained=%d/%d, explosive=%d, other_fail=%d, mean total ATE=%.1f, SD(rep means)=%.1f\n",
                  nrow(bD), bootstrap_ate$fit_D$n_attempted,
                  bootstrap_ate$fit_D$n_explosive,
                  bootstrap_ate$fit_D$n_error,
                  mean(bD$ate_total_mean, na.rm = TRUE),
                  stats::sd(bD$ate_total_mean, na.rm = TRUE)))
    } else {
      cat(sprintf("  Bootstrap D complete: success=%d, fail=%d\n",
                  bootstrap_ate$fit_D$n_success, bootstrap_ate$fit_D$n_fail))
    }
  }

  cat("Retaining pre-bootstrap checkpoints as requested.\n")
  bootstrap_elapsed <- proc.time()[["elapsed"]] - t_bootstrap
}
if (is.finite(bootstrap_elapsed)) {
  add_timing_row("bootstrap_ate_total", bootstrap_elapsed, "ok",
                 sprintf("targets=%s reps=%d", paste(boot_targets_run, collapse = ","), BOOT_N_REPS))
} else if (RUN_BOOTSTRAP_ATE && BOOT_N_REPS > 0L) {
  add_timing_row("bootstrap_ate_total", NA_real_, "skipped",
                 "requested but no eligible targets")
} else {
  add_timing_row("bootstrap_ate_total", NA_real_, "skipped", "disabled")
}

if (isTRUE(BOOTSTRAP_ONLY)) {
  if (is.null(bootstrap_ate)) {
    stop("BOOTSTRAP_ONLY: bootstrap_ate is NULL after Step 8.")
  }
  # Assemble a publication RDS that looks like a single completed analysis run
  # (no bootstrap-only / patch provenance markers).
  base_file <- normalizePath(BOOTSTRAP_PATCH_FILE, winslash = "/", mustWork = TRUE)
  bak <- paste0(
    tools::file_path_sans_ext(base_file),
    "_pre_bootstrap_",
    format(Sys.time(), "%Y%m%d%H%M%S"),
    ".rds"
  )
  ok_bak <- file.copy(base_file, bak, overwrite = FALSE)
  cat(sprintf("Backup before writing final results: %s (copied=%s)\n", bak, ok_bak))
  results_final <- readRDS(base_file)

  # Recompute point ATEs under THIS job's ATE flags so bias-correction
  # recenters bootstrap clouds to matching scenario centers.
  cat(sprintf(
    "\n--- BOOTSTRAP_ONLY: recomputing point ATE under %s ---\n",
    OK_ATE_METHOD_LABEL
  ))
  ensure_ate_psock_pool()
  ate_crn_base_seed_bootonly <- if (is.finite(OK_ATE_CRN_BASE) && !is.na(OK_ATE_CRN_BASE)) {
    as.integer(OK_ATE_CRN_BASE)
  } else if (is.finite(OK_GLOBAL_SEED) && !is.na(OK_GLOBAL_SEED)) {
    as.integer(100000L + 1000L * OK_GLOBAL_SEED)
  } else {
    100000L
  }
  recompute_point_ate <- function(biv_params, marg, observed_data, label) {
    if (isTRUE(OK_ATE_BIVARIATE) && !is.null(biv_params)) {
      ate_estim_bivariate(
        biv_params = biv_params,
        windowT = windowT_ate,
        windowS = win_km,
        state_spaces_obs = state_spaces,
        label = label,
        n_sims = ATE_N_SIMS,
        n_cores = ATE_SIM_CORES,
        m0 = ETAS_M0,
        beta_gr = BETA_GR,
        filtration_history = if (isTRUE(OK_ATE_CONDITIONAL_ON_PRE)) pp_pre else NULL,
        t_trunc = SEM_T_TRUNC_DAYS,
        n_tiles = partition$n,
        crn_base_seed = ate_crn_base_seed_bootonly,
        use_crn = OK_ATE_USE_CRN,
        crn_pair = OK_ATE_CRN_PAIR,
        contrast = OK_ATE_CONTRAST,
        covariate_lookup = KDE_BG_LOOKUP
      )
    } else {
      ate_estim_fast(
        marg$ctrl, marg$treat, observed_data, label,
        filtration_history = if (isTRUE(OK_ATE_CONDITIONAL_ON_PRE)) pp_pre else NULL,
        crn_base_seed = ate_crn_base_seed_bootonly,
        contrast = OK_ATE_CONTRAST,
        covariate_lookup = KDE_BG_LOOKUP,
        t_trunc = SEM_T_TRUNC_DAYS
      )
    }
  }
  ate_E <- recompute_point_ate(
    E_params, E_marginals, pp_post,
    sprintf("Fit C point (%s)", OK_ATE_METHOD_LABEL)
  )
  ate_F <- recompute_point_ate(
    F_params, F_marginals, pp_post,
    sprintf("Fit D point (%s)", OK_ATE_METHOD_LABEL)
  )
  # Write into public primary slots C/D; also refresh legacy E/F-as-primary if present.
  if (!is.null(results_final$fits_named$C)) {
    results_final$fits_named$C$ate <- ate_E
    results_final$fits_named$C$params <- E_params
  } else if (!is.null(results_final$fits_named$E)) {
    results_final$fits_named$E$ate <- ate_E
    results_final$fits_named$E$params <- E_params
  }
  if (!is.null(results_final$fits_named$D)) {
    results_final$fits_named$D$ate <- ate_F
    results_final$fits_named$D$params <- F_params
  } else if (!is.null(results_final$fits_named$F)) {
    results_final$fits_named$F$ate <- ate_F
    results_final$fits_named$F$params <- F_params
  }
  if (!is.null(results_final$fitE)) results_final$fitE$ate <- ate_E
  if (!is.null(results_final$fitF)) results_final$fitF$ate <- ate_F

  # Match the bootstrap_ate$config shape written by a full analysis run.
  bootstrap_ate$config <- list(
    reps = BOOT_N_REPS,
    targets = boot_targets_run,
    targets_requested = BOOT_TARGETS,
    targets_run = boot_targets_run,
    refit_scope = BOOT_REFIT_SCOPE,
    outer_cores = BOOT_OUTER_CORES,
    sem_inner_iter = BOOT_SEM_INNER_ITER,
    seed = BOOT_SEED,
    simulation_law = if (isTRUE(OK_ATE_BIVARIATE)) {
      "bivariate_post_conditional_on_observed_pre"
    } else {
      "univariate_post_conditional_on_observed_pre"
    },
    conditioning_history = "observed_retained_pre_treatment",
    n_conditioning_events = nrow(pp_pre),
    refit_law = if (isTRUE(OK_ATE_BIVARIATE)) "bivariate" else "univariate",
    stability_fit_margin = ETAS_BRANCHING_MAX,
    explosive_drop_threshold = 1,
    bias_correction_stage = "after_stability_filter",
    ate_bivariate = isTRUE(OK_ATE_BIVARIATE),
    ate_contrast = OK_ATE_CONTRAST,
    ate_method = OK_ATE_METHOD_LABEL,
    ate_scenario = if (nzchar(OK_ATE_SCENARIO)) OK_ATE_SCENARIO else NA_character_
  )
  bootstrap_ate$note <- NULL
  results_final$bootstrap_ate <- bootstrap_ate

  if (is.null(results_final$config)) results_final$config <- list()
  results_final$config$ATE_BIVARIATE <- isTRUE(OK_ATE_BIVARIATE)
  results_final$config$ATE_CONTRAST <- OK_ATE_CONTRAST
  results_final$config$ATE_METHOD <- OK_ATE_METHOD_LABEL
  results_final$config$ATE_SCENARIO <- if (nzchar(OK_ATE_SCENARIO)) OK_ATE_SCENARIO else NULL
  results_final$config$ATE_N_SIMS <- ATE_N_SIMS
  results_final$config$BOOT_N_REPS <- BOOT_N_REPS
  results_final$config$BOOT_REFIT_SCOPE <- BOOT_REFIT_SCOPE
  results_final$config$BOOT_SEM_INNER_ITER <- BOOT_SEM_INNER_ITER
  results_final$config$BOOT_OUTER_CORES <- BOOT_OUTER_CORES
  results_final$config$BOOT_TARGETS <- BOOT_TARGETS
  results_final$config$ETAS_BRANCHING_MAX <- ETAS_BRANCHING_MAX
  results_final$config$ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS <- ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS
  results_final$config$ETAS_P_LOWER_BOUND <- ETAS_P_LOWER_BOUND
  results_final$config$ETAS_Q_LOWER_BOUND <- ETAS_Q_LOWER_BOUND
  results_final$config$BOOT_BRANCHING_MAX <- BOOT_BRANCHING_MAX
  results_final$config$RUN_BOOTSTRAP_ATE <- TRUE
  # Strip any interim recompute/bootstrap-only markers so the artifact matches
  # a one-shot paper run.
  results_final$config$BOOTSTRAP_ONLY <- NULL
  results_final$config$ATE_BIVARIATE_PATCHED <- NULL
  results_final$ate_bivariate_patch <- NULL
  if (is.list(results_final$metadata)) {
    results_final$metadata$bivariate_bootstrap <- NULL
    results_final$metadata$paper_promote <- NULL
    results_final$metadata$stage <- "final"
    results_final$metadata$saved_at <- as.character(Sys.time())
    results_final$metadata$job_id <- Sys.getenv("SLURM_JOB_ID", "")
    results_final$metadata$ate_scenario <- if (nzchar(OK_ATE_SCENARIO)) OK_ATE_SCENARIO else NA_character_
    results_final$metadata$ate_method <- OK_ATE_METHOD_LABEL
  } else {
    results_final$metadata <- list(
      stage = "final",
      saved_at = as.character(Sys.time()),
      job_id = Sys.getenv("SLURM_JOB_ID", ""),
      ate_scenario = if (nzchar(OK_ATE_SCENARIO)) OK_ATE_SCENARIO else NA_character_,
      ate_method = OK_ATE_METHOD_LABEL
    )
  }

  out_tagged <- file.path(OUT_DIR, add_file_tag("oklahoma_results.rds"))
  saveRDS(results_final, out_tagged)
  scenario_outs <- character(0)
  if (nzchar(OK_ATE_SCENARIO)) {
    scen_dir <- file.path(OUT_DIR, "ate_scenarios", OK_ATE_SCENARIO)
    dir.create(scen_dir, recursive = TRUE, showWarnings = FALSE)
    scen_paper <- file.path(scen_dir, "for_paper.rds")
    scen_alias <- file.path(OUT_DIR, paste0("for_paper_", OK_ATE_SCENARIO, ".rds"))
    saveRDS(results_final, scen_paper)
    saveRDS(results_final, scen_alias)
    writeLines(
      c(
        paste0("scenario=", OK_ATE_SCENARIO),
        paste0("ate_method=", OK_ATE_METHOD_LABEL),
        paste0("ate_bivariate=", isTRUE(OK_ATE_BIVARIATE)),
        paste0("ate_contrast=", OK_ATE_CONTRAST),
        paste0("job_id=", Sys.getenv("SLURM_JOB_ID", "")),
        paste0("results=", out_tagged)
      ),
      file.path(scen_dir, "scenario.txt")
    )
    scenario_outs <- c(scen_paper, scen_alias)
  } else {
    for_paper_out <- file.path(OUT_DIR, "for_paper.rds")
    saveRDS(results_final, for_paper_out)
    scenario_outs <- for_paper_out
    # Keep the named input path in sync when it is also under the output tree.
    if (!identical(normalizePath(base_file, winslash = "/", mustWork = FALSE),
                   normalizePath(for_paper_out, winslash = "/", mustWork = FALSE)) &&
        !identical(normalizePath(base_file, winslash = "/", mustWork = FALSE),
                   normalizePath(out_tagged, winslash = "/", mustWork = FALSE))) {
      saveRDS(results_final, base_file)
    }
  }
  cat(sprintf("Wrote final paper results (full-run shape):\n  %s\n", out_tagged))
  for (p in scenario_outs) cat(sprintf("  %s\n", p))
  point_mean <- function(ate) {
    if (is.null(ate$all_nothing_sim$total_saved)) return(NA_real_)
    mean(ate$all_nothing_sim$total_saved, na.rm = TRUE)
  }
  cat(sprintf("  Point ATE C=%.1f D=%.1f method=%s\n",
              point_mean(ate_E), point_mean(ate_F), OK_ATE_METHOD_LABEL))
  boot_c_blk <- if (!is.null(bootstrap_ate$fit_C)) bootstrap_ate$fit_C else bootstrap_ate$fit_E
  boot_d_blk <- if (!is.null(bootstrap_ate$fit_D)) bootstrap_ate$fit_D else bootstrap_ate$fit_F
  if (!is.null(boot_c_blk$replicate_summary) &&
      nrow(boot_c_blk$replicate_summary) > 0) {
    cat(sprintf("  C boot mean=%.1f sd=%.1f (n=%d)\n",
                mean(boot_c_blk$replicate_summary$ate_total_mean, na.rm = TRUE),
                stats::sd(boot_c_blk$replicate_summary$ate_total_mean, na.rm = TRUE),
                nrow(boot_c_blk$replicate_summary)))
  }
  if (!is.null(boot_d_blk$replicate_summary) &&
      nrow(boot_d_blk$replicate_summary) > 0) {
    cat(sprintf("  D boot mean=%.1f sd=%.1f (n=%d)\n",
                mean(boot_d_blk$replicate_summary$ate_total_mean, na.rm = TRUE),
                stats::sd(boot_d_blk$replicate_summary$ate_total_mean, na.rm = TRUE),
                nrow(boot_d_blk$replicate_summary)))
  }
  cat("\n=== Oklahoma bootstrap stage finished successfully ===\n")
  quit(save = "no", status = 0)
}

stop_ate_psock_pool()

# Sensitivity ATE payloads are already computed in Step 6a so they are
# included in pre-bootstrap checkpoints and reused for final results.

# ============================================================================
# 7. Summary tables
# ============================================================================
cat("\n")
cat("===========================================================================\n")
cat("                       PARAMETER COMPARISON TABLE\n")
cat("===========================================================================\n\n")

fmt <- function(x) {
  x <- suppressWarnings(as.numeric(x)[1L])
  if (!is.finite(x)) return(sprintf("%8s", "NULL"))
  sprintf("%8.4f", x)
}

biv_ctrl_val <- function(par, p) {
  if (is.null(par)) return(NA_real_)
  if (p == "mu") par[["mu_0"]] else par[[paste0(p, "_00")]]
}
biv_treat_val <- function(par, p) {
  if (is.null(par)) return(NA_real_)
  if (p == "mu") par[["mu_1"]] else par[[paste0(p, "_11")]]
}

cat(sprintf("%-12s  %8s  %8s  |  %8s  %8s  |  %8s  %8s  |  %8s  %8s\n",
  "", "A.ctrl", "A.treat", "B.ctrl", "B.treat",
  "C.ctrl", "C.treat", "D.ctrl", "D.treat"))
cat(paste(rep("-", 106), collapse = ""), "\n")

for (p in c("mu", "A", "alpha_m")) {
  cat(sprintf("%-12s  %8s  %8s  |  %8s  %8s  |  %8s  %8s  |  %8s  %8s\n",
    p,
    fmt(biv_ctrl_val(B_params, p)), fmt(biv_treat_val(B_params, p)),
    fmt(biv_ctrl_val(D_params, p)), fmt(biv_treat_val(D_params, p)),
    fmt(biv_ctrl_val(E_params, p)), fmt(biv_treat_val(E_params, p)),
    fmt(biv_ctrl_val(F_params, p)), fmt(biv_treat_val(F_params, p))))
}

cat("\nCross-excitation:\n")
for (info in list(list(par = B_params, lab = "A"),
                  list(par = D_params, lab = "B"),
                  list(par = E_params, lab = "C"),
                  list(par = F_params, lab = "D"))) {
  if (is.null(info$par)) {
    cat(sprintf("  %s: NULL (fit failed)\n", info$lab))
  } else {
    cat(sprintf("  %s: A_01=%.4f  A_10=%.4f  alpha_m_01=%.4f  alpha_m_10=%.4f\n",
      info$lab, info$par[["A_01"]], info$par[["A_10"]],
      info$par[["alpha_m_01"]], info$par[["alpha_m_10"]]))
  }
}

cat("\n")
cat("===========================================================================\n")
cat(sprintf("     EARTHQUAKES SAVED COMPARISON (%d-day horizon)\n", ATE_WINDOW_DAYS))
cat("===========================================================================\n\n")
cat(sprintf("%-30s  %12s  %12s  %12s  %12s\n",
            "", "Total Saved", "SD(saved)", "eta_ctrl", "eta_treat"))
cat(paste(rep("-", 90), collapse = ""), "\n")

ate_print_list <- list(
  list(ate = ate_B, lab = sprintf("A: Naive Biv (%s)", PRIMARY_PARTITION_LABEL)),
  list(ate = ate_D, lab = sprintf("B: SEM Biv (%s)", PRIMARY_PARTITION_LABEL)),
  list(ate = ate_E, lab = sprintf("C: Naive Biv+KDE (%s)", PRIMARY_PARTITION_LABEL)),
  list(ate = ate_F, lab = sprintf("D: SEM Biv+KDE (%s)", PRIMARY_PARTITION_LABEL)))

for (pname in names(ate_partitions)) {
  ap <- ate_partitions[[pname]]
  if (!is.null(ap$ate_E))
    ate_print_list[[length(ate_print_list) + 1]] <- list(ate = ap$ate_E, lab = sprintf("C: %s", pname))
  if (!is.null(ap$ate_F))
    ate_print_list[[length(ate_print_list) + 1]] <- list(ate = ap$ate_F, lab = sprintf("D: %s", pname))
}

for (nm in ate_print_list) {
  if (!is.null(nm$ate)) {
    total <- mean(nm$ate$all_nothing_sim$total_saved, na.rm = TRUE)
    total_sd <- stats::sd(nm$ate$all_nothing_sim$total_saved, na.rm = TRUE)
    eta_c <- if (!is.null(nm$ate$analytic)) nm$ate$analytic$eta_ctrl else NA
    eta_t <- if (!is.null(nm$ate$analytic)) nm$ate$analytic$eta_treat else NA
    cat(sprintf("%-30s  %12.0f  %12.0f  %12.3f  %12.3f\n",
      nm$lab, total, total_sd, eta_c, eta_t))
  } else {
    cat(sprintf("%-30s  %12s  %12s  %12s  %12s\n",
      nm$lab, "FAILED", "FAILED", "FAILED", "FAILED"))
  }
}
cat(sprintf("\nTotal Saved: expected all-or-nothing earthquake reduction (control minus treatment) over %d days in Oklahoma.\n",
            ATE_WINDOW_DAYS))

# ============================================================================
# 7b. Louis-method SEs for SEM-based ATE estimates
# ============================================================================
cat("\n--- Step 7b: Louis-method SE skipped for this run ---\n")
louis_D <- NULL
louis_F <- NULL

# ============================================================================
# 8. Save results
# ============================================================================
cat("\n--- Saving results ---\n")

analysis_elapsed_pre_save <- proc.time()[["elapsed"]] - analysis_start_elapsed
add_timing_row("analysis_total_pre_save", analysis_elapsed_pre_save, "ok")
timing_df <- if (length(timing_rows) > 0L) {
  do.call(rbind, timing_rows)
} else {
  data.frame(
    order = integer(0),
    stage = character(0),
    elapsed_sec = numeric(0),
    elapsed_min = numeric(0),
    status = character(0),
    detail = character(0),
    stringsAsFactors = FALSE
  )
}
timing_df <- timing_df[order(timing_df$order), , drop = FALSE]
rownames(timing_df) <- NULL

# Print timing summary to stdout so it is visible in slurm.out.
cat("\n--- Runtime timing summary (chronological) ---\n")
if (nrow(timing_df) > 0L) {
  for (i in seq_len(nrow(timing_df))) {
    row_i <- timing_df[i, , drop = FALSE]
    sec_i <- suppressWarnings(as.numeric(row_i$elapsed_sec[[1]]))
    min_i <- suppressWarnings(as.numeric(row_i$elapsed_min[[1]]))
    sec_txt <- if (is.finite(sec_i)) sprintf("%.1f", sec_i) else "NA"
    min_txt <- if (is.finite(min_i)) sprintf("%.2f", min_i) else "NA"
    status_txt <- as.character(row_i$status[[1]])
    detail_txt <- trimws(as.character(row_i$detail[[1]]))
    if (nzchar(detail_txt) && !identical(detail_txt, "NA")) {
      cat(sprintf("  %02d | %-28s | %8ss (%6sm) | %-8s | %s\n",
                  as.integer(row_i$order[[1]]), as.character(row_i$stage[[1]]),
                  sec_txt, min_txt, status_txt, detail_txt))
    } else {
      cat(sprintf("  %02d | %-28s | %8ss (%6sm) | %-8s\n",
                  as.integer(row_i$order[[1]]), as.character(row_i$stage[[1]]),
                  sec_txt, min_txt, status_txt))
    }
  }
  finite_idx <- is.finite(timing_df$elapsed_sec)
  if (any(finite_idx)) {
    top_idx <- order(timing_df$elapsed_sec[finite_idx], decreasing = TRUE)
    top_n <- min(5L, length(top_idx))
    top_rows <- timing_df[which(finite_idx)[top_idx[seq_len(top_n)]], , drop = FALSE]
    cat("\n--- Slowest stages (top 5) ---\n")
    for (i in seq_len(nrow(top_rows))) {
      cat(sprintf("  %-28s : %8.1fs (%6.2fm)\n",
                  as.character(top_rows$stage[[i]]),
                  as.numeric(top_rows$elapsed_sec[[i]]),
                  as.numeric(top_rows$elapsed_min[[i]])))
    }
  }
} else {
  cat("  No timing rows recorded.\n")
}

results <- list(
  fitB = list(params = B_params, loglik = B_loglik, fit = fitB, ate = ate_B),
  fitD = list(params = D_params, ctrl = D_ctrl, treat = D_treat, sem = semD, ate = ate_D, louis = louis_D),
  # Legacy aliases: fitE/fitF hold primary all_free content (= fits_named$C/$D).
  fitE = list(params = E_params, loglik = E_loglik, fit = fitE, ate = ate_E),
  fitF = list(params = F_params, ctrl = F_ctrl, treat = F_treat, sem = semF, ate = ate_F, louis = louis_F),
  fitI = list(params = I_params, fit = fitI, ate = ate_I),
  fitJ = list(params = J_params, fit = semJ, ctrl = J_ctrl, treat = J_treat, ate = ate_J),
  fits_named = fits_named,
  kde_variants = kde_variant_fits,
  fitG = NULL,
  fitH = NULL,
  bootstrap_ate = bootstrap_ate,
  fit_variability = fit_variability,
  partition_results = partition_results,
  ate_partitions = ate_partitions,
  kde_bandwidth_sensitivity = kde_bandwidth_sensitivity,
  t_trunc_sensitivity = t_trunc_sensitivity,
  kde_info = kde_info,
  control_snapshot_fits = control_snapshot_fits,
  sem_pilot_tuning = sem_pilot_tuning,
  plots = analysis_plots,
  pp_data = list(pp_pre = pp_pre, pp_pre_holdout = pp_pre_holdout, pp_post = pp_post),
  timing_df = timing_df,
  timing_info = list(
    analysis_start = as.character(analysis_start_time),
    analysis_end_pre_save = as.character(Sys.time())
  ),
  counties = list(
    names = counties_sf_valid$NAME,
    treated_names = county_treated_names,
    n_counties = county_partition$n,
    n_treated = sum(county_treated_idx)
  ),
  fit_name_map = list(
    A = "fits_named$A",
    B = "fits_named$B",
    C = "fits_named$C",
    D = "fits_named$D",
    E = "fits_named$E",
    F = "fits_named$F",
    G = "fits_named$G",
    H = "fits_named$H",
    I = "fits_named$I",
    J = "fits_named$J"
  ),
  config = list(
    ETAS_M0 = ETAS_M0, BETA_GR = BETA_GR,
    BETA_GR_SOURCE = BETA_GR_SOURCE,
    GAMMA_FIXED = GAMMA_FIXED,
    FIXED_STRUCTURAL = FIXED_STRUCTURAL,
    SEM_N_ITER = SEM_N_ITER, SEM_INNER_ITER = SEM_INNER_ITER,
    SENS_SEM_INNER_ITER = SENS_SEM_INNER_ITER,
    SEM_INNER_PROPS = SEM_INNER_PROPS,
    SEM_N_LABELLINGS = SEM_N_LABELLINGS,
    SEM_OUTER_MAXIT = SEM_OUTER_MAXIT,
    SEM_OUTER_MAXIT_BIV = SEM_OUTER_MAXIT_BIV,
    SEM_WARMSTART_FIXED = SEM_WARMSTART_FIXED,
    SEM_CHANGE_FACTOR = SEM_CHANGE_FACTOR,
    SEM_CHANGE_FACTOR_MIN_MULT = SEM_CHANGE_FACTOR_MIN_MULT,
    SEM_CHANGE_FACTOR_MAX_MULT = SEM_CHANGE_FACTOR_MAX_MULT,
    SEM_OPTIM_METHOD = SEM_OPTIM_METHOD,
    SEM_MONOTONE_COMPLETE_LL = SEM_MONOTONE_COMPLETE_LL,
    SEM_START_FROM_C = SEM_START_FROM_C,
    SEM_BIV_N_THREADS = SEM_BIV_N_THREADS,
    SEM_SINGLE_FLIP_FROM_ITER = SEM_SINGLE_FLIP_FROM_ITER,
    SEM_SELECTION_TEMPERATURE = SEM_SELECTION_TEMPERATURE,
    RUN_SEM_PILOT = RUN_SEM_PILOT,
    SEM_PILOT_INNER_ITER = SEM_PILOT_INNER_ITER,
    SEM_PILOT_CORES = SEM_PILOT_CORES,
    SEM_PILOT_MAX_COMBOS = SEM_PILOT_MAX_COMBOS,
    SEM_TEMPORAL_WEIGHT = SEM_TEMPORAL_WEIGHT,
    SEM_TEMPORAL_SCALE_DAYS = SEM_TEMPORAL_SCALE_DAYS,
    SEM_T_TRUNC_DAYS = SEM_T_TRUNC_DAYS,
    SEM_T_TRUNC_SOURCE = SEM_T_TRUNC_SOURCE,
    SEM_T_TRUNC_REL = SEM_T_TRUNC_REL,
    T_TRUNC_SENS_DAYS = T_TRUNC_SENS_DAYS,
    RUN_T_TRUNC_SENSITIVITY = RUN_T_TRUNC_SENSITIVITY,
    RUN_SENSITIVITY = RUN_SENSITIVITY,
    RUN_FIT_VARIABILITY = RUN_FIT_VARIABILITY,
    FIT_VARIABILITY_REPS = FIT_VARIABILITY_REPS,
    FIT_VARIABILITY_CORES = FIT_VARIABILITY_CORES,
    KDE_VARIANT_MODE = KDE_VARIANT_MODE,
    KDE_BW_METHOD = KDE_BW_METHOD,
    PRIMARY_PARTITION = PRIMARY_PARTITION,
    KDE_BW_SENS_KM = KDE_BW_SENS_KM,
    RUN_PARTITION_SENSITIVITY = RUN_PARTITION_SENSITIVITY,
    CD_ONLY = CD_ONLY,
    UNIV_KDE_ONLY = UNIV_KDE_ONLY,
    RUN_KDE_PROFILE_SWEEP = RUN_KDE_PROFILE_SWEEP,
    MEMORY_SAFE = MEMORY_SAFE,
    TRIM_SENS_OBJECTS = TRIM_SENS_OBJECTS,
    SENS_CORES = SENS_CORES,
    SENS_BIV_N_THREADS = SENS_BIV_N_THREADS,
    ATE_SIM_CORES = ATE_SIM_CORES,
    ATE_USE_CRN = OK_ATE_USE_CRN,
    ATE_CRN_PAIR = OK_ATE_CRN_PAIR,
    ATE_CONDITIONAL_ON_PRE = OK_ATE_CONDITIONAL_ON_PRE,
    ATE_CRN_BASE = OK_ATE_CRN_BASE,
    SEM_MAX_RELABEL_STEP_FRAC = SEM_MAX_RELABEL_STEP_FRAC,
    SEM_FORCE_PARAM_UPDATE_FLIP_FRAC = SEM_FORCE_PARAM_UPDATE_FLIP_FRAC,
    ATE_N_SIMS = ATE_N_SIMS, ATE_WINDOW_DAYS = ATE_WINDOW_DAYS,
    ATE_BIVARIATE = OK_ATE_BIVARIATE, ATE_CONTRAST = OK_ATE_CONTRAST,
    RUN_BOOTSTRAP_ATE = RUN_BOOTSTRAP_ATE,
    BOOT_N_REPS = BOOT_N_REPS,
    BOOT_TARGETS = BOOT_TARGETS,
    SENS_TARGETS = SENS_TARGETS,
    BOOT_REFIT_SCOPE = BOOT_REFIT_SCOPE,
    BOOT_SEM_INNER_ITER = BOOT_SEM_INNER_ITER,
    BOOT_OUTER_CORES = BOOT_OUTER_CORES,
    BOOT_SEED = BOOT_SEED,
    BOOT_IDENTICAL_RANDOMNESS = OK_BOOT_IDENTICAL_RANDOMNESS,
    BOOT_GUARD_DEGENERATE = OK_BOOT_GUARD_DEGENERATE,
    ETAS_BRANCHING_MAX = ETAS_BRANCHING_MAX,
    ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS = ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS,
    ETAS_P_LOWER_BOUND = ETAS_P_LOWER_BOUND,
    ETAS_Q_LOWER_BOUND = ETAS_Q_LOWER_BOUND,
    BOOT_BRANCHING_MAX = BOOT_BRANCHING_MAX,
    BOOT_MAX_PRE_EVENTS = BOOT_MAX_PRE_EVENTS,
    BOOT_MAX_POST_EVENTS_PER_PROC = BOOT_MAX_POST_EVENTS_PER_PROC,
    BOOT_MAX_TOTAL_EVENTS = BOOT_MAX_TOTAL_EVENTS,
    TEST_MODE = TEST_MODE,
    windowT_post = windowT_post,
    n_pre = nrow(pp_pre), n_pre_holdout = nrow(pp_pre_holdout), n_pre_total = nrow(pp_pre_all),
    n_post = nrow(pp_post),
    n_counties = county_partition$n,
    n_treated_counties = sum(county_treated_idx),
    n_tiles = partition$n,
    n_treated = sum(treated_idx),
    grid_diameters = grid_diameters,
    grid_multipliers = grid_multipliers,
    kde_bandwidth_multipliers = vapply(kde_bandwidth_specs, `[[`, numeric(1), "multiplier"),
    spatial_scale_D = D_scale
  )
)

out_file <- file.path(OUT_DIR, add_file_tag("oklahoma_results.rds"))
saveRDS(results, out_file)
cat(sprintf("Results saved to: %s\n", out_file))
cat("Plots stored in results payload under results$plots (PNG file output disabled).\n")

# Keep report outputs synchronized with the newest results.
report_file <- file.path(SCRIPT_DIR, "oklahoma_report.qmd")
if (file.exists(report_file) && length(REPORT_FORMATS) > 0L) {
  quarto_bin <- Sys.which("quarto")
  if (nzchar(quarto_bin)) {
    cat(sprintf("\n--- Rendering report (%s) ---\n", paste(REPORT_FORMATS, collapse = ", ")))
    old_wd <- getwd()
    old_results_file_env <- Sys.getenv("OK_RESULTS_FILE", unset = NA_character_)
    Sys.setenv(OK_RESULTS_FILE = normalizePath(out_file, winslash = "/", mustWork = FALSE))
    setwd(dirname(report_file))
    render_errors <- character(0)
    for (fmt in REPORT_FORMATS) {
      render_status <- tryCatch(
        system2(
          quarto_bin,
          c(
            "render", basename(report_file), "--to", fmt,
            "--output-dir", path_rel_to(REPORT_OUT_DIR, SCRIPT_DIR)
          ),
          stdout = TRUE,
          stderr = TRUE
        ),
        error = function(e) e
      )
      if (inherits(render_status, "error")) {
        render_errors <- c(render_errors, sprintf("[%s] %s", fmt, render_status$message))
      } else if (!is.null(attr(render_status, "status"))) {
        render_errors <- c(
          render_errors,
          sprintf("[%s] exit code %s", fmt, as.character(attr(render_status, "status")))
        )
      }
    }
    setwd(old_wd)
    if (is.na(old_results_file_env)) {
      Sys.unsetenv("OK_RESULTS_FILE")
    } else {
      Sys.setenv(OK_RESULTS_FILE = old_results_file_env)
    }
    if (length(render_errors) > 0L) {
      cat("Report render failed:\n")
      cat(paste0("  - ", render_errors, collapse = "\n"), "\n")
    } else {
      cat("Report render complete.\n")
      sync_stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
      result_mtime <- tryCatch(file.info(out_file)$mtime[[1]], error = function(e) NA)
      result_stamp <- if (!is.na(result_mtime)) {
        format(result_mtime, "%Y%m%d_%H%M%S")
      } else {
        format(Sys.time(), "%Y%m%d_%H%M%S")
      }
      report_html_rendered <- if ("html" %in% REPORT_FORMATS) {
        file.path(REPORT_OUT_DIR, "oklahoma_report.html")
      } else {
        NA_character_
      }
      report_pdf_rendered <- if ("pdf" %in% REPORT_FORMATS) {
        file.path(REPORT_OUT_DIR, "oklahoma_report.pdf")
      } else {
        NA_character_
      }
      report_html_target <- if (is.character(report_html_rendered) &&
                                !is.na(report_html_rendered)) {
        file.path(REPORT_OUT_DIR, add_file_tag(sprintf("oklahoma_report_%s.html", result_stamp)))
      } else {
        report_html_rendered
      }
      report_pdf_target <- if (is.character(report_pdf_rendered) &&
                               !is.na(report_pdf_rendered)) {
        file.path(REPORT_OUT_DIR, add_file_tag(sprintf("oklahoma_report_%s.pdf", result_stamp)))
      } else {
        report_pdf_rendered
      }
      if (is.character(report_html_rendered) && !is.na(report_html_rendered) &&
          is.character(report_html_target) && !is.na(report_html_target) &&
          !identical(report_html_rendered, report_html_target)) {
        file.copy(report_html_rendered, report_html_target, overwrite = TRUE)
      }
      if (is.character(report_pdf_rendered) && !is.na(report_pdf_rendered) &&
          is.character(report_pdf_target) && !is.na(report_pdf_target) &&
          !identical(report_pdf_rendered, report_pdf_target)) {
        file.copy(report_pdf_rendered, report_pdf_target, overwrite = TRUE)
      }
      report_html_path <- if ("html" %in% REPORT_FORMATS) {
        normalizePath(report_html_target, winslash = "/", mustWork = FALSE)
      } else {
        NA_character_
      }
      report_pdf_path <- if ("pdf" %in% REPORT_FORMATS) {
        normalizePath(report_pdf_target, winslash = "/", mustWork = FALSE)
      } else {
        NA_character_
      }
      stamp_lines <- c(
        paste0("sync_stamp: ", sync_stamp),
        paste0("results_file: ", normalizePath(out_file, winslash = "/", mustWork = FALSE)),
        paste0("report_html: ", report_html_path),
        paste0("report_pdf: ", report_pdf_path)
      )
      stamp_root <- file.path(OUT_DIR, add_file_tag("last_run_sync_stamp.txt"))
      writeLines(stamp_lines, stamp_root)
      cat("Sync stamp written for Google Drive change detection.\n")
      if ("html" %in% REPORT_FORMATS) {
        cat(sprintf("Timestamped report (from results mtime %s): %s\n",
                    result_stamp, normalizePath(report_html_target, winslash = "/", mustWork = FALSE)))
      }
      if ("pdf" %in% REPORT_FORMATS) {
        cat(sprintf("Timestamped PDF (from results mtime %s): %s\n",
                    result_stamp, normalizePath(report_pdf_target, winslash = "/", mustWork = FALSE)))
      }
      remove_legacy_report_outputs()
    }
  } else {
    cat("Quarto not found in PATH; skipping report render.\n")
  }
} else if (file.exists(report_file)) {
  cat("Report render skipped (OK_REPORT_FORMATS empty).\n")
}

cat("\n=== Oklahoma Analysis Complete ===\n")
