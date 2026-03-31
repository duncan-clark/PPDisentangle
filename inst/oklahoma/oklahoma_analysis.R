#!/usr/bin/env Rscript
# ============================================================================
# Oklahoma Induced Seismicity Analysis — Bivariate ETAS with Sensitivity
#
# Treatment: OCC directive AOI_20150318 (wastewater injection reduction)
# State space: Oklahoma state boundary
#
# Core fits (county partition):
#   B. Naive bivariate ETAS     — joint MLE on location-labeled data
#   D. SEM bivariate ETAS       — adaptive_SEM with model_type="etas_bivariate"
#   E. Naive bivariate + KDE bg — B with non-parametric background from control
#   F. SEM bivariate + KDE bg   — D with non-parametric background from control
#
# Partition sensitivity:
#   - County boundaries (default)
#   - Grid: 1D, 2D, 5D (D = estimated triggering range)
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
# Canonical output path at repo root:
#   output/oklahoma/
OUT_DIR  <- file.path(REPO_DIR, "output", "oklahoma")
PLOT_DIR <- file.path(OUT_DIR, "plots")
for (d in unique(c(OUT_DIR, PLOT_DIR))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
analysis_plots <- list()
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
ETAS_M0    <- 2.5
BETA_GR    <- 2.3
CRS_PROJ   <- 5070

VANILLA_MAXIT <- if (QUICK_CHECK) 120 else if (TEST_MODE) 500 else 1000
VANILLA_STARTS <- list(
  list(mu = 1.0, A = 0.2, alpha_m = 0.8, c = 0.05, p = 1.2,
       D = 5.0, gamma = 0.5, q = 1.5),
  list(mu = 2.0, A = 0.1, alpha_m = 0.5, c = 0.1,  p = 1.3,
       D = 2.0, gamma = 0.3, q = 2.0)
)

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
SEM_CHANGE_FACTOR <- 0.01
SEM_CHANGE_FACTOR_MIN_MULT <- suppressWarnings(as.numeric(Sys.getenv("OK_SEM_CHANGE_FACTOR_MIN_MULT", "0.2")))
if (!is.finite(SEM_CHANGE_FACTOR_MIN_MULT) || is.na(SEM_CHANGE_FACTOR_MIN_MULT) || SEM_CHANGE_FACTOR_MIN_MULT <= 0) {
  SEM_CHANGE_FACTOR_MIN_MULT <- 0.2
}
SEM_CHANGE_FACTOR_MAX_MULT <- suppressWarnings(as.numeric(Sys.getenv("OK_SEM_CHANGE_FACTOR_MAX_MULT", "2.0")))
if (!is.finite(SEM_CHANGE_FACTOR_MAX_MULT) || is.na(SEM_CHANGE_FACTOR_MAX_MULT) || SEM_CHANGE_FACTOR_MAX_MULT < SEM_CHANGE_FACTOR_MIN_MULT) {
  SEM_CHANGE_FACTOR_MAX_MULT <- max(2.0, SEM_CHANGE_FACTOR_MIN_MULT)
}
SEM_OPTIM_METHOD <- tolower(trimws(Sys.getenv("OK_SEM_OPTIM_METHOD", "sample_weighted")))
if (!SEM_OPTIM_METHOD %in% c("max", "sample_weighted")) SEM_OPTIM_METHOD <- "sample_weighted"
SEM_SELECTION_TEMPERATURE <- suppressWarnings(as.numeric(Sys.getenv("OK_SEM_SELECTION_TEMPERATURE", "0.08")))
if (!is.finite(SEM_SELECTION_TEMPERATURE) || is.na(SEM_SELECTION_TEMPERATURE) || SEM_SELECTION_TEMPERATURE <= 0) {
  SEM_SELECTION_TEMPERATURE <- 0.08
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
SEM_T_TRUNC_DAYS_USER <- suppressWarnings(as.numeric(Sys.getenv("OK_SEM_T_TRUNC_DAYS", "")))
if (!is.finite(SEM_T_TRUNC_DAYS_USER) || is.na(SEM_T_TRUNC_DAYS_USER) || SEM_T_TRUNC_DAYS_USER <= 0) {
  SEM_T_TRUNC_DAYS_USER <- NA_real_
}
SEM_T_TRUNC_REL <- suppressWarnings(as.numeric(Sys.getenv("OK_SEM_T_TRUNC_REL", "0.05")))
if (!is.finite(SEM_T_TRUNC_REL) || is.na(SEM_T_TRUNC_REL) || SEM_T_TRUNC_REL <= 0 || SEM_T_TRUNC_REL >= 1) {
  SEM_T_TRUNC_REL <- 0.05
}
SEM_T_TRUNC_DAYS <- SEM_T_TRUNC_DAYS_USER
SEM_T_TRUNC_SOURCE <- if (is.finite(SEM_T_TRUNC_DAYS_USER) && !is.na(SEM_T_TRUNC_DAYS_USER)) "env" else "auto_from_pre50"
SEM_PARAM_UPDATE  <- if (QUICK_CHECK) 10 else if (TEST_MODE) 10 else 25
SEM_OUTER_MAXIT       <- if (QUICK_CHECK) 40 else if (TEST_MODE) 200 else 220
SEM_OUTER_MAXIT_BIV   <- 1000
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
KDE_VARIANT_MODE <- tolower(trimws(Sys.getenv("OK_KDE_VARIANT_MODE", "triple")))
if (!KDE_VARIANT_MODE %in% c("single", "triple")) KDE_VARIANT_MODE <- "triple"
RUN_KDE_PROFILE_SWEEP <- identical(KDE_VARIANT_MODE, "triple")
OK_VERBOSE <- tolower(Sys.getenv("OK_VERBOSE", "false")) %in% c("1", "true", "yes", "y")
DF_VERBOSE <- tolower(Sys.getenv("OK_DF_VERBOSE", "false")) %in% c("1", "true", "yes", "y")
SEM_WORKER_LOGS <- tolower(Sys.getenv("OK_SEM_WORKER_LOGS", "true")) %in% c("1", "true", "yes", "y")
SEM_WORKER_LOG_VERBOSE <- tolower(Sys.getenv("OK_SEM_WORKER_LOG_VERBOSE", if (SEM_WORKER_LOGS) "true" else "false")) %in% c("1", "true", "yes", "y")
SEM_WORKER_LOG_SPLIT <- tolower(Sys.getenv("OK_SEM_WORKER_LOG_SPLIT", "false")) %in% c("1", "true", "yes", "y")

STRUCT_DEFAULTS <- list(c = 0.05, p = 1.2, D = 5.0, gamma = 0.5, q = 1.5)
# Structural terms are fixed downstream after first-half pre-treatment calibration.
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
# Default to bootstrapping the all-free KDE pair (E and F).
BOOT_TARGETS_RAW <- toupper(Sys.getenv("OK_BOOT_TARGETS", "E,F"))
BOOT_TARGETS <- unique(trimws(unlist(strsplit(BOOT_TARGETS_RAW, ","))))
BOOT_TARGETS <- BOOT_TARGETS[BOOT_TARGETS %in% c("E", "F")]
if (length(BOOT_TARGETS) < 1) BOOT_TARGETS <- c("E", "F")
# Sensitivity targets are wired directly from bootstrap targets.
SENS_TARGETS <- unique(BOOT_TARGETS)
if (!setequal(SENS_TARGETS, c("E", "F"))) {
  stop(sprintf(
    "Sensitivity/bootstrap target alignment requires E,F; got: %s",
    paste(SENS_TARGETS, collapse = ",")
  ))
}
# Safer default is no per-replicate refit (still allows explicit partial/full override).
BOOT_REFIT_SCOPE <- tolower(trimws(Sys.getenv("OK_BOOT_REFIT_SCOPE", "none")))
if (!BOOT_REFIT_SCOPE %in% c("none", "partial", "full")) BOOT_REFIT_SCOPE <- "none"
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
BOOT_BRANCHING_MAX <- suppressWarnings(as.numeric(Sys.getenv("OK_BOOT_BRANCHING_MAX", "0.98")))
if (!is.finite(BOOT_BRANCHING_MAX) || is.na(BOOT_BRANCHING_MAX) || BOOT_BRANCHING_MAX <= 0) {
  BOOT_BRANCHING_MAX <- 0.98
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
BOOT_OUTER_DEFAULT <- if (MEMORY_SAFE) max(2L, min(8L, as.integer(floor(N_CORES / 4L)))) else N_CORES
BOOT_OUTER_CORES <- suppressWarnings(as.integer(ifelse(nzchar(BOOT_OUTER_CORES_RAW), BOOT_OUTER_CORES_RAW, as.character(BOOT_OUTER_DEFAULT))))
if (!is.finite(BOOT_OUTER_CORES) || is.na(BOOT_OUTER_CORES) || BOOT_OUTER_CORES < 1L) BOOT_OUTER_CORES <- 1L
BOOT_OUTER_CORES <- max(1L, min(BOOT_OUTER_CORES, N_CORES))
if (MEMORY_SAFE && RUN_BOOTSTRAP_ATE && BOOT_N_REPS > 0L && BOOT_OUTER_CORES > 1L) {
  boot_outer_cap <- suppressWarnings(as.integer(Sys.getenv("OK_BOOT_OUTER_CAP_MEMSAFE", as.character(max(2L, min(8L, as.integer(floor(N_CORES / 4L))))))))
  if (!is.finite(boot_outer_cap) || is.na(boot_outer_cap) || boot_outer_cap < 1L) boot_outer_cap <- 2L
  if (BOOT_OUTER_CORES > boot_outer_cap) {
    cat(sprintf("Memory-safe bootstrap cap: reducing BOOT_OUTER_CORES from %d to %d\n", BOOT_OUTER_CORES, boot_outer_cap))
    BOOT_OUTER_CORES <- boot_outer_cap
  }
}
SENS_CORES_DEFAULT <- if (MEMORY_SAFE) min(2L, N_CORES) else N_CORES
SENS_CORES <- suppressWarnings(as.integer(ifelse(nzchar(SENS_CORES_RAW), SENS_CORES_RAW, as.character(SENS_CORES_DEFAULT))))
if (!is.finite(SENS_CORES) || is.na(SENS_CORES) || SENS_CORES < 1L) SENS_CORES <- 1L
SENS_CORES <- max(1L, min(SENS_CORES, N_CORES))
ATE_SIM_CORES_DEFAULT <- 1L
ATE_SIM_CORES <- suppressWarnings(as.integer(ifelse(nzchar(ATE_SIM_CORES_RAW), ATE_SIM_CORES_RAW, as.character(ATE_SIM_CORES_DEFAULT))))
if (!is.finite(ATE_SIM_CORES) || is.na(ATE_SIM_CORES) || ATE_SIM_CORES < 1L) ATE_SIM_CORES <- 1L
ATE_SIM_CORES <- max(1L, min(ATE_SIM_CORES, N_CORES))
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
    all(vapply(checks, function(x) all(isTRUE(x)), logical(1)))
  }
  if (identical(export_mode, "all")) {
    worker_globals <- ls(envir = .GlobalEnv, all.names = TRUE)
  } else {
    # Export only symbols referenced by FUN that exist in .GlobalEnv.
    globs <- tryCatch(
      codetools::findGlobals(FUN, merge = FALSE),
      error = function(e) list(variables = character(0), functions = character(0))
    )
    needed <- unique(c(globs$variables, globs$functions))
    needed <- needed[nzchar(needed)]
    worker_globals <- needed[vapply(
      needed,
      function(nm) exists(nm, envir = .GlobalEnv, inherits = FALSE),
      logical(1)
    )]
  }
  if (length(worker_globals) > 0L) {
    parallel::clusterExport(cl, varlist = worker_globals, envir = .GlobalEnv)
  }
  if (!verify_worker_symbols(worker_globals) && identical(export_mode, "minimal")) {
    warning(sprintf(
      "PSOCK preflight missing symbols for %s under minimal export; retrying with full global export.",
      label
    ))
    worker_globals <- ls(envir = .GlobalEnv, all.names = TRUE)
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
cat(sprintf("SEM selection: method=%s | temperature=%.3f\n",
            SEM_OPTIM_METHOD, SEM_SELECTION_TEMPERATURE))
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
cat(sprintf("SEM inner iters: main=%d, sensitivity=%d, bootstrap=%d\n",
            SEM_INNER_ITER, SENS_SEM_INNER_ITER, BOOT_SEM_INNER_ITER))
cat(sprintf("Targets (shared): sensitivity=%s | bootstrap=%s\n",
            paste(SENS_TARGETS, collapse = ","),
            paste(BOOT_TARGETS, collapse = ",")))
cat(sprintf("SEM warm-start fixed adaptive step: %s\n", SEM_WARMSTART_FIXED))
sem_t_trunc_banner <- if (is.finite(SEM_T_TRUNC_DAYS) && !is.na(SEM_T_TRUNC_DAYS) && SEM_T_TRUNC_DAYS > 0) {
  as.character(signif(SEM_T_TRUNC_DAYS, 4))
} else {
  sprintf("auto(from pre-50 estimate, rel=%.3f)", SEM_T_TRUNC_REL)
}
cat(sprintf("SEM proposal/event t_trunc (days): %s\n", sem_t_trunc_banner))
cat(sprintf("SEM temporal relabel weight: %.3f\n", SEM_TEMPORAL_WEIGHT))
cat(sprintf("B/D SEM verbose tracing: %s\n", DF_VERBOSE))
cat(sprintf("Verbose optimizer/SEM tracing: %s\n", OK_VERBOSE))
cat(sprintf("SEM worker logs enabled: %s | worker logs verbose: %s | split-to-main-log: %s\n",
            SEM_WORKER_LOGS, SEM_WORKER_LOG_VERBOSE, SEM_WORKER_LOG_SPLIT))
cat(sprintf("KDE variant mode: %s\n", KDE_VARIANT_MODE))

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

if (!dir.exists(DATA_DIR)) {
  cat("Data not found — running data preparation script...\n")
  source(file.path(SCRIPT_DIR, "Oklahoma_data_and_viz.R"), local = TRUE)
}

ev_all <- fread(file.path(DATA_DIR, "events_all.csv"))
meta   <- jsonlite::fromJSON(readLines(file.path(DATA_DIR, "metadata.json")))

t_star_utc <- as.POSIXct(meta$design$t_star_utc, tz = "UTC")
ev_all[, time_utc := as.POSIXct(time_utc, tz = "UTC")]
ev_all[, t_days := as.numeric(difftime(time_utc, t_star_utc, units = "days"))]
ev_all[, x_km := x_m / 1000]
ev_all[, y_km := y_m / 1000]
if (!"mag" %in% names(ev_all)) stop("No 'mag' column in events data")
ev_all <- ev_all[mag >= ETAS_M0]

post_end_days <- as.numeric(difftime(
  as.POSIXct(meta$design$post_end_utc, tz = "UTC"), t_star_utc, units = "days"))

# ============================================================================
# 2. Oklahoma county tessellation
# ============================================================================
cat("\n--- Step 2: Building county tessellation ---\n")

options(tigris_use_cache = TRUE)
counties_sf <- counties(state = "OK", cb = TRUE, year = 2022)
counties_sf <- st_transform(counties_sf, CRS_PROJ)
counties_sf <- st_make_valid(counties_sf)

ok_boundary <- st_union(counties_sf)
ok_boundary <- st_make_valid(ok_boundary)

bb <- st_bbox(ok_boundary)
win_km <- owin(xrange = c(bb["xmin"], bb["xmax"]) / 1000,
               yrange = c(bb["ymin"], bb["ymax"]) / 1000)

# Convert county polygons to spatstat owin tiles
county_owins <- lapply(seq_len(nrow(counties_sf)), function(i) {
  geom <- st_geometry(counties_sf[i, ])
  coords_list <- st_coordinates(geom)
  x_km <- coords_list[, 1] / 1000
  y_km <- coords_list[, 2] / 1000
  tryCatch(
    owin(poly = list(x = rev(x_km), y = rev(y_km))),
    error = function(e) {
      tryCatch(
        owin(poly = list(x = x_km, y = y_km)),
        error = function(e2) NULL
      )
    }
  )
})

valid_idx <- !sapply(county_owins, is.null)
if (sum(valid_idx) < 50) {
  cat("  Warning: only", sum(valid_idx), "of", nrow(counties_sf), "counties converted.\n")
}

county_owins_valid <- county_owins[valid_idx]
counties_sf_valid  <- counties_sf[valid_idx, ]

names(county_owins_valid) <- counties_sf_valid$NAME
partition <- tess(tiles = county_owins_valid, window = win_km)

cat(sprintf("  Counties in tessellation: %d / %d\n",
            partition$n, nrow(counties_sf)))

# Assign treatment: county centroid inside OCC AOI
aoi_path <- file.path(DATA_DIR, "occ_aoi_layer_2.geojson")
aoi_sf <- st_read(aoi_path, quiet = TRUE)
aoi_sf <- st_transform(aoi_sf, CRS_PROJ)
aoi_sf <- st_make_valid(aoi_sf)
aoi_union <- st_union(aoi_sf)

county_centroids <- st_centroid(counties_sf_valid)
inside_aoi <- lengths(st_within(county_centroids, aoi_union)) > 0

partition_processes <- ifelse(inside_aoi, "treated", "control")
names(partition_processes) <- counties_sf_valid$NAME

treated_idx <- partition_processes == "treated"
treated_names <- names(partition_processes)[treated_idx]
control_ss <- as.owin(partition[!treated_idx])
treated_ss <- as.owin(partition[treated_idx])
state_spaces <- list(control = control_ss, treated = treated_ss)

cat(sprintf("  Treated counties: %d (%s)\n",
            sum(treated_idx),
            paste(treated_names, collapse = ", ")))
cat(sprintf("  Control counties: %d\n", sum(!treated_idx)))

# ---- Assign events to counties ----
assign_county <- function(df) {
  ti <- as.integer(tileindex(df$x, df$y, partition))
  df$location_process <- ifelse(is.na(ti), NA_character_,
                                partition_processes[pmin(pmax(ti, 1), partition$n)])
  df$W <- 1.0
  df$n <- nrow(df)
  df$background <- TRUE
  df
}

pp_pre  <- assign_county(as.data.frame(ev_all[t_days < 0,
  .(x = x_km, y = y_km, t = t_days, mag = mag)]))
pp_post <- assign_county(as.data.frame(ev_all[t_days >= 0 & t_days <= post_end_days,
  .(x = x_km, y = y_km, t = t_days, mag = mag)]))

# Drop events outside Oklahoma
pp_pre  <- pp_pre[!is.na(pp_pre$location_process), ]
pp_post <- pp_post[!is.na(pp_post$location_process), ]

# Pre-treatment split:
# - hold out first 50% (earliest events) for KDE background estimation
# - use second 50% for model estimation / carry-over in SEM
pp_pre_all <- pp_pre[order(pp_pre$t), ]
n_pre_total <- nrow(pp_pre_all)
n_pre_holdout <- floor(n_pre_total * 0.5)
if (n_pre_total > 0 && n_pre_holdout < 1) n_pre_holdout <- 1
holdout_idx <- if (n_pre_holdout > 0) seq_len(n_pre_holdout) else integer(0)
keep_idx <- if (n_pre_total > n_pre_holdout) (n_pre_holdout + 1):n_pre_total else integer(0)
if (length(keep_idx) < 1) {
  stop("Pre-treatment split left zero estimation events; cannot continue.")
}
pp_pre_holdout <- pp_pre_all[holdout_idx, ]
pp_pre <- pp_pre_all[keep_idx, ]

# Carry-over convention: all pre-treatment events are control-process
# events regardless of eventual treated/control location status.
pp_pre$process  <- "control"
pp_post$process <- pp_post$location_process
pp_pre$inferred_process  <- "control"
pp_post$inferred_process <- pp_post$location_process
pp_all <- rbind(pp_pre, pp_post)
pp_all <- pp_all[order(pp_all$t), ]

windowT_post <- c(0, post_end_days)
windowT_fit <- c(min(pp_pre$t, na.rm = TRUE), post_end_days)
pp_pre_ctrl_all <- pp_pre_all[pp_pre_all$location_process == "control", ]
pp_pre_ctrl_all <- pp_pre_ctrl_all[order(pp_pre_ctrl_all$t), ]
n_pre_ctrl_total <- nrow(pp_pre_ctrl_all)
n_pre_ctrl_struct_init <- floor(n_pre_ctrl_total * 0.5)
if (n_pre_ctrl_total > 0 && n_pre_ctrl_struct_init < 1) n_pre_ctrl_struct_init <- 1
ctrl_struct_idx <- if (n_pre_ctrl_struct_init > 0) seq_len(n_pre_ctrl_struct_init) else integer(0)
pp_pre_ctrl_struct_init <- pp_pre_ctrl_all[ctrl_struct_idx, ]
n_pre_ctrl_struct_init <- nrow(pp_pre_ctrl_struct_init)

if (n_pre_ctrl_struct_init < 5) {
  stop("Insufficient first-half control pre-treatment events for structural parameter estimation.")
}

# Data-driven ETAS calibration from first-half control pre data.
# These values are used as robust initializations for all downstream fits.
PRE_CTRL_BOOT_PARAMS <- NULL
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
        maxit = VANILLA_MAXIT,
        fixed_params = NULL,
        zero_background_region = treated_ss
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_try) && is.finite(fit_try$value) && fit_try$value > best_val) {
      best <- fit_try
      best_val <- fit_try$value
    }
  }
  if (is.null(best) || is.null(best$par)) return(STRUCT_DEFAULTS)
  pre_full <- as.list(best$par)
  pre_needed <- c("mu", "A", "alpha_m", "c", "p", "D", "gamma", "q")
  for (nm in pre_needed) {
    if (is.null(pre_full[[nm]]) || !is.finite(pre_full[[nm]])) {
      pre_full[[nm]] <- if (nm %in% names(STRUCT_DEFAULTS)) STRUCT_DEFAULTS[[nm]] else VANILLA_STARTS[[1]][[nm]]
    }
  }
  PRE_CTRL_BOOT_PARAMS <<- pre_full[pre_needed]
  out <- as.list(best$par[c("c", "p", "D", "gamma", "q")])
  out <- out[!vapply(out, is.null, logical(1))]
  needed <- c("c", "p", "D", "gamma", "q")
  for (nm in needed) {
    if (is.null(out[[nm]]) || !is.finite(out[[nm]])) out[[nm]] <- STRUCT_DEFAULTS[[nm]]
  }
  out
}
STRUCT_INIT <- estimate_structural_init()
FIXED_STRUCTURAL <- as.list(STRUCT_INIT[c("c", "p", "D", "gamma", "q")])
if (is.null(PRE_CTRL_BOOT_PARAMS)) {
  PRE_CTRL_BOOT_PARAMS <- list(mu = 1.0, A = 0.2, alpha_m = 0.8,
                               c = STRUCT_INIT$c, p = STRUCT_INIT$p,
                               D = STRUCT_INIT$D, gamma = STRUCT_INIT$gamma, q = STRUCT_INIT$q)
}
compute_temporal_trunc_from_pre <- function(c_param, p_param, rel_level = SEM_T_TRUNC_REL) {
  c_param <- suppressWarnings(as.numeric(c_param))
  p_param <- suppressWarnings(as.numeric(p_param))
  rel_level <- suppressWarnings(as.numeric(rel_level))
  if (!is.finite(c_param) || !is.finite(p_param) || !is.finite(rel_level) ||
      c_param <= 0 || p_param <= 0 || rel_level <= 0 || rel_level >= 1) {
    return(NULL)
  }
  # Solve ((t + c) / c)^(-p) = rel_level for t.
  t_cut <- c_param * ((rel_level)^(-1 / p_param) - 1)
  if (!is.finite(t_cut) || t_cut <= 0) return(NULL)
  as.numeric(t_cut)
}
if (!is.finite(SEM_T_TRUNC_DAYS) || is.na(SEM_T_TRUNC_DAYS) || SEM_T_TRUNC_DAYS <= 0) {
  SEM_T_TRUNC_DAYS <- compute_temporal_trunc_from_pre(PRE_CTRL_BOOT_PARAMS$c, PRE_CTRL_BOOT_PARAMS$p, SEM_T_TRUNC_REL)
  if (is.null(SEM_T_TRUNC_DAYS)) {
    SEM_T_TRUNC_SOURCE <- "none"
  } else {
    SEM_T_TRUNC_SOURCE <- sprintf("auto_from_pre50(rel=%.3f)", SEM_T_TRUNC_REL)
  }
} else {
  SEM_T_TRUNC_SOURCE <- "env"
}

cat(sprintf("  Structural init (from first 50%% control pre): c=%.4f, p=%.4f, D=%.4f, gamma=%.4f, q=%.4f\n",
            STRUCT_INIT$c, STRUCT_INIT$p, STRUCT_INIT$D,
            STRUCT_INIT$gamma, STRUCT_INIT$q))
cat(sprintf("  Pre-treatment full ETAS init: mu=%.4f, A=%.4f, alpha_m=%.4f, c=%.4f, p=%.4f, D=%.4f, gamma=%.4f, q=%.4f\n",
            PRE_CTRL_BOOT_PARAMS$mu, PRE_CTRL_BOOT_PARAMS$A, PRE_CTRL_BOOT_PARAMS$alpha_m,
            PRE_CTRL_BOOT_PARAMS$c, PRE_CTRL_BOOT_PARAMS$p, PRE_CTRL_BOOT_PARAMS$D,
            PRE_CTRL_BOOT_PARAMS$gamma, PRE_CTRL_BOOT_PARAMS$q))
cat(sprintf("  SEM t_trunc resolved from %s: %s days\n",
            SEM_T_TRUNC_SOURCE,
            if (is.null(SEM_T_TRUNC_DAYS)) "none" else as.character(signif(SEM_T_TRUNC_DAYS, 5))))

apply_pre_init_etas <- function(start_par) {
  out <- start_par
  for (nm in etas_names) out[[nm]] <- PRE_CTRL_BOOT_PARAMS[[nm]]
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
        maxit = VANILLA_MAXIT,
        fixed_params = NULL
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

# ============================================================================
# 3. Partition and point pattern plots
# ============================================================================
cat("\n--- Step 3: Plots ---\n")

tryCatch({
  counties_plot <- counties_sf_valid
  counties_plot$treatment <- partition_processes

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
               fixed_params = NULL, zero_background_region = zbr),
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
fit_b <- function() {
  tryCatch({
    fit_etas_bivariate(
      params_init = biv_init, realiz = pp_all,
      windowT = windowT_fit, windowS = win_km, m0 = ETAS_M0,
      control_state_space = control_ss, treated_state_space = treated_ss,
      treated_background_zero_before = 0,
      maxit = VANILLA_MAXIT, fixed_params = FIXED_STRUCTURAL, trace = 0,
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
biv_fixed <- FIXED_STRUCTURAL

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
                        fixed_params_in = biv_fixed,
                        background_rate_var_in = NULL,
                        use_pre_history_for_biv_in = TRUE,
                        treated_background_zero_before_in = 0,
                        sem_t_trunc_in = SEM_T_TRUNC_DAYS,
                        sem_inner_iter_in = SEM_INNER_ITER,
                        sem_change_factor_in = SEM_CHANGE_FACTOR,
                        sem_change_factor_min_mult_in = SEM_CHANGE_FACTOR_MIN_MULT,
                        sem_change_factor_max_mult_in = SEM_CHANGE_FACTOR_MAX_MULT,
                        sem_optim_method_in = SEM_OPTIM_METHOD,
                        sem_selection_temperature_in = SEM_SELECTION_TEMPERATURE,
                        verbose_in = DF_VERBOSE,
                        label = "SEM",
                        sem_rng_label_in = NULL) {
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
        proposal_method = "simulation",
        outer_maxit = SEM_OUTER_MAXIT, outer_maxit_biv = SEM_OUTER_MAXIT_BIV
      ),
      m0 = ETAS_M0, beta_gr = BETA_GR,
      etas_bivariate_params = if (identical(model_type_in, "etas_bivariate")) init_params_sem else NULL,
      background_rate_var = background_rate_var_in,
      use_pre_history_for_biv = use_pre_history_for_biv_in,
      treated_background_zero_before = treated_background_zero_before_in,
      t_trunc = sem_t_trunc_in
    )
  }
  run_sem_body <- function() {
    tryCatch({
      pp_data_sem <- pp_data_in
      init_params_sem <- init_params_in
      if (identical(model_type_in, "etas_bivariate") && !is.null(init_params_sem)) {
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

B_params <- biv_init
B_loglik <- NA_real_
D_params <- biv_init_D
D_ctrl <- A_ctrl
D_treat <- A_treat

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
bw_sigma <- 2 * bw_diggle
lambda_im <- density(X_bg, sigma = bw_sigma, edge = TRUE, at = "pixels")
min_nz <- min(lambda_im$v[lambda_im$v > 0], na.rm = TRUE)
lambda_im$v[lambda_im$v <= 0] <- min_nz
cat(sprintf("  KDE bandwidth: %.2f km (from %d held-out training events)\n",
            as.numeric(bw_sigma), n_pre_holdout_ctrl))

# Data-informed triggering range from the same subset used for KDE:
# use a robust upper-local scale (90th percentile nearest-neighbor distance).
nn_training <- nndist(X_bg)
trigger_range_km <- as.numeric(stats::quantile(nn_training, probs = 0.9, na.rm = TRUE))
if (!is.finite(trigger_range_km) || trigger_range_km <= 0) {
  trigger_range_km <- STRUCT_DEFAULTS$D
}
cat(sprintf("  Estimated triggering range (Q90 NN distance): %.2f km\n", trigger_range_km))

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

  pp_pre_bg_ctrl <- normalize_bg_weights(pp_pre[pp_pre$location_process == "control", ],
                                         ctrl_ss, lambda_covariate)
  pp_pre_bg_treat <- normalize_bg_weights(pp_pre[pp_pre$location_process == "treated", ],
                                          treat_ss, lambda_covariate)
  pp_all_bg_local <- rbind(pp_pre_bg_ctrl$new_df, pp_pre_bg_treat$new_df, pp_post_bg_local)
  pp_all_bg_local <- pp_all_bg_local[order(pp_all_bg_local$t), ]
  list(pp_post_bg = pp_post_bg_local, pp_all_bg = pp_all_bg_local)
}

bg_ctrl_res <- normalize_bg_weights(pp_post[pp_post$location_process == "control", ],
                                     control_ss, lambda_im)
bg_treat_res <- normalize_bg_weights(pp_post[pp_post$location_process == "treated", ],
                                      treated_ss, lambda_im)
pp_post_bg <- rbind(bg_ctrl_res$new_df, bg_treat_res$new_df)
pp_post_bg <- pp_post_bg[order(pp_post_bg$t), ]

pp_pre_bg_ctrl <- normalize_bg_weights(pp_pre[pp_pre$location_process == "control", ],
                                        control_ss, lambda_im)
pp_pre_bg_treat <- normalize_bg_weights(pp_pre[pp_pre$location_process == "treated", ],
                                         treated_ss, lambda_im)
pp_all_bg <- rbind(pp_pre_bg_ctrl$new_df, pp_pre_bg_treat$new_df, pp_post_bg)
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
         subtitle = "Estimated from held-out first 50% of pre-treatment events",
         x = "X (km)", y = "Y (km)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  store_analysis_plot("background_rate_kde", p_bg)
  cat("  Stored KDE background-rate plot in results$plots$background_rate_kde\n")
}, error = function(e) cat("  Background-rate plot error:", e$message, "\n"))

# Parameter profiles for county KDE bivariate fits.
# - all_free: no fixed parameters
# - control_only_fixed: structural terms fixed from first-half control pre-treatment;
#   productivity terms (including mu_0 and A_00) are free.
# - productivity_free: all mu/A/alpha free; structural terms fixed from pre-treatment
kde_variant_specs <- list(
  all_free = list(
    id = "all_free",
    label = "all parameters free",
    fixed_params = NULL
  ),
  control_only_fixed = list(
    id = "control_only_fixed",
    label = "structural fixed from first-half control pre-treatment",
    fixed_params = FIXED_STRUCTURAL
  ),
  productivity_free = list(
    id = "productivity_free",
    label = "mu/A/alpha free, structural fixed",
    fixed_params = FIXED_STRUCTURAL
  )
)
kde_primary_variant_id <- "control_only_fixed"
kde_variant_ids <- if (RUN_KDE_PROFILE_SWEEP) names(kde_variant_specs) else kde_primary_variant_id
kde_variant_fits <- list(E = list(), F = list())
# Sensitivity and bootstrap refits align with the all-free KDE pair (E/F).
SENSITIVITY_FIXED_PARAMS <- kde_variant_specs$all_free$fixed_params
KDE_FIT_LETTERS <- list(
  control_only_fixed = list(E = "C", F = "D"),
  all_free = list(E = "E", F = "F"),
  productivity_free = list(E = "G", F = "H")
)
kde_fit_label <- function(fit_type, variant_id) {
  letter <- KDE_FIT_LETTERS[[variant_id]][[fit_type]]
  sprintf("Fit %s [%s]", letter, variant_id)
}
# KDE bandwidth SEM (all-free F) uses this RNG key so stochastic SEM matches main Fit F
# when sigma and data match (i.e. digglex2 vs main county KDE).
OK_BW_SEM_RNG_LABEL <- kde_fit_label("F", "all_free")
cat(sprintf(
  "  KDE county fit variants to run: %s\n",
  paste(kde_variant_ids, collapse = ", ")
))
cat(sprintf(
  "  Sensitivity fixed profile: %s\n",
  if (is.null(SENSITIVITY_FIXED_PARAMS) || length(SENSITIVITY_FIXED_PARAMS) < 1) "<none (all-free)>" else paste(names(SENSITIVITY_FIXED_PARAMS), collapse = ", ")
))

# ============================================================================
# 4C. Fit C: Naive bivariate ETAS with KDE background
# ============================================================================
cat("\n--- Fit C: Naive bivariate ETAS with KDE background ---\n")

biv_init_E <- apply_pre_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
fit_e <- function(init_params = biv_init_E,
                  fixed_params = FIXED_STRUCTURAL,
                  fit_label = "Fit C") {
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
                  fixed_params = FIXED_STRUCTURAL,
                  fit_label = "Fit D") {
  tryCatch({
    init_params_use <- enforce_control_decay_start(init_params, A_ctrl)
    cat(sprintf(
      "  [%s] start decays from control fit: c=%.4f p=%.4f D=%.4f gamma=%.4f q=%.4f\n",
      fit_label,
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
      fixed_params_in = fixed_params,
      background_rate_var_in = "W",
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
if (RUN_SEM_PILOT) {
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

cat("\n--- Step 4 unified dispatch: running all county fits in parallel ---\n")
fit_jobs_all <- list(
  list(kind = "A_hom_naive", variant_id = NA_character_),
  list(kind = "B_hom_sem", variant_id = NA_character_)
)
for (vid in kde_variant_ids) {
  fit_jobs_all[[length(fit_jobs_all) + 1L]] <- list(kind = "C_kde_naive", variant_id = vid)
  fit_jobs_all[[length(fit_jobs_all) + 1L]] <- list(kind = "D_kde_sem", variant_id = vid)
}
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
    out <- fit_f(
      init_params = biv_init_F,
      fixed_params = spec$fixed_params,
      fit_label = fit_label
    )
  }
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
if (N_CORES > 1L && length(fit_jobs_all) > 1L) {
  fit_all_out <- run_parallel(
    fit_jobs_all, run_all_fit_job,
    cores = min(length(fit_jobs_all), N_CORES),
    label = "fit-all-county-jobs"
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

B_params <- if (!is.null(fitB)) fitB$par else biv_init
B_loglik <- if (!is.null(fitB)) fitB$value else NA_real_
cat("  Fit A params:", paste(biv_names, round(B_params, 4), sep = "=", collapse = ", "), "\n")
if (!is.null(semD)) {
  D_params <- semD$etas_bivariate_params
  D_ctrl <- semD$hawkes_params_control
  D_treat <- semD$hawkes_params_treated
  cat("  Fit B params:", paste(biv_names, round(D_params, 4), sep = "=", collapse = ", "), "\n")
} else {
  D_params <- biv_init_D
  D_ctrl <- A_ctrl
  D_treat <- A_treat
  cat("  Fit B failed, falling back to naive initialization.\n")
}

fitE <- NULL
semF <- NULL
for (vid in kde_variant_ids) {
  spec <- kde_variant_specs[[vid]]
  row_e <- get_fit_out("C_kde_naive", vid)
  row_f <- get_fit_out("D_kde_sem", vid)
  fitE_var <- if (!is.null(row_e)) row_e$obj else NULL
  semF_var <- if (!is.null(row_f)) row_f$obj else NULL
  E_var_params <- if (!is.null(fitE_var)) fitE_var$par else biv_init_E
  E_var_loglik <- if (!is.null(fitE_var)) fitE_var$value else NA_real_
  F_var_params <- if (!is.null(semF_var)) semF_var$etas_bivariate_params else biv_init_F
  F_var_ctrl <- if (!is.null(semF_var)) semF_var$hawkes_params_control else A_ctrl
  F_var_treat <- if (!is.null(semF_var)) semF_var$hawkes_params_treated else A_treat

  kde_variant_fits$E[[vid]] <- list(
    id = spec$id,
    label = spec$label,
    fixed_params = spec$fixed_params,
    fit = fitE_var,
    params = E_var_params,
    objective = E_var_loglik
  )
  kde_variant_fits$F[[vid]] <- list(
    id = spec$id,
    label = spec$label,
    fixed_params = spec$fixed_params,
    fit = semF_var,
    params = F_var_params,
    hawkes_params_control = F_var_ctrl,
    hawkes_params_treated = F_var_treat
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
    cat("  Fit C params:", paste(biv_names, round(E_var_params, 4), sep = "=", collapse = ", "), "\n")
    cat("  Fit D params:", paste(biv_names, round(F_var_params, 4), sep = "=", collapse = ", "), "\n")
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
    cat(sprintf("  %s params: %s\n", kde_fit_label("E", spec$id),
                paste(biv_names, round(E_var_params, 4), sep = "=", collapse = ", ")))
    cat(sprintf("  %s params: %s\n", kde_fit_label("F", spec$id),
                paste(biv_names, round(F_var_params, 4), sep = "=", collapse = ", ")))
  }
}
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
# 4I/4J. Univariate ETAS comparison (no cross-excitation)
# ============================================================================
cat("\n--- Step 4I/4J: Univariate ETAS (I/J) dispatch ---\n")

fit_indep_pair <- function(pp_data_in,
                           background_rate_var_in = NULL,
                           fixed_params_in = NULL,
                           ctrl_init_in = A_ctrl,
                           treat_init_in = A_treat,
                           fit_label = "Fit I") {
  dat <- as.data.frame(pp_data_in)
  dat <- dat[dat$t >= 0, , drop = FALSE]
  if (nrow(dat) < 1L) return(NULL)
  ctrl_dat <- dat[dat$location_process == "control", , drop = FALSE]
  treat_dat <- dat[dat$location_process == "treated", , drop = FALSE]
  if (nrow(ctrl_dat) < 1L || nrow(treat_dat) < 1L) return(NULL)
  ctrl_init <- apply_pre_init_etas(ctrl_init_in)
  treat_init <- apply_pre_init_etas(treat_init_in)
  fit_ctrl <- tryCatch(
    fit_etas(
      params_init = ctrl_init, realiz = ctrl_dat,
      windowT = windowT_post, windowS = win_km, m0 = ETAS_M0,
      maxit = VANILLA_MAXIT, fixed_params = fixed_params_in,
      zero_background_region = treated_ss,
      background_rate_var = background_rate_var_in,
      t_trunc = SEM_T_TRUNC_DAYS
    ),
    error = function(e) {
      cat(sprintf("  [%s] control independent fit error: %s\n", fit_label, e$message))
      NULL
    }
  )
  fit_treat <- tryCatch(
    fit_etas(
      params_init = treat_init, realiz = treat_dat,
      windowT = windowT_post, windowS = win_km, m0 = ETAS_M0,
      maxit = VANILLA_MAXIT, fixed_params = fixed_params_in,
      zero_background_region = control_ss,
      background_rate_var = background_rate_var_in,
      t_trunc = SEM_T_TRUNC_DAYS
    ),
    error = function(e) {
      cat(sprintf("  [%s] treated independent fit error: %s\n", fit_label, e$message))
      NULL
    }
  )
  ctrl_par <- if (!is.null(fit_ctrl) && !is.null(fit_ctrl$par)) as.list(fit_ctrl$par) else ctrl_init
  treat_par <- if (!is.null(fit_treat) && !is.null(fit_treat$par)) as.list(fit_treat$par) else treat_init
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
      # Keep structural terms free so finite-moment constraints are enforced
      # by the ETAS likelihood's soft-to-hard boundary penalties.
      fixed_params_in = NULL,
      ctrl_init_in = A_ctrl,
      treat_init_in = A_treat,
      fit_label = "Fit I"
    ),
    error = function(e) {
      cat("  Fit I (naive univariate) error:", e$message, "\n")
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
      # Match other SEM fits: enforce finite moments via likelihood penalties,
      # not by hard-fixing structural terms.
      fixed_params_in = NULL,
      background_rate_var_in = "W",
      sem_inner_iter_in = SEM_INNER_ITER,
      verbose_in = DF_VERBOSE,
      label = "Fit J"
    ),
    error = function(e) {
      cat("  Fit J (SEM univariate) error:", e$message, "\n")
      NULL
    }
  )
}

univ_jobs <- c("I", "J")
run_one_univ_job <- function(tag) {
  t0 <- proc.time()[["elapsed"]]
  out_obj <- NULL
  cat(sprintf("    [univ-job:%s] start pid=%d mem=%s\n", tag, Sys.getpid(), mem_snapshot()))
  if (tag == "I") out_obj <- fit_i()
  if (tag == "J") out_obj <- fit_j()
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
fitI <- if (!is.null(row_i)) row_i$obj else NULL
semJ <- if (!is.null(row_j)) row_j$obj else NULL
fit_I_elapsed <- if (!is.null(row_i)) row_i$elapsed else NA_real_
fit_J_elapsed <- if (!is.null(row_j)) row_j$elapsed else NA_real_
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
I_ctrl <- if (!is.null(fitI)) fitI$control else A_ctrl
I_treat <- if (!is.null(fitI)) fitI$treated else A_treat
J_ctrl <- if (!is.null(semJ)) semJ$hawkes_params_control else A_ctrl
J_treat <- if (!is.null(semJ)) semJ$hawkes_params_treated else A_treat
I_params <- list(control = I_ctrl, treated = I_treat)
J_params <- list(control = J_ctrl, treated = J_treat)
pp_post_sem_J <- if (!is.null(semJ) && !is.null(semJ$adaptive$adaptive_labelling)) {
  semJ$adaptive$adaptive_labelling
} else {
  pp_post_bg
}
pp_post_sem_J <- pp_post_sem_J[pp_post_sem_J$t >= 0, , drop = FALSE]

# ============================================================================
# 4I. KDE bandwidth sensitivity (county partition, inhomogeneous E/F)
# ============================================================================
cat("\n--- Step 4I: KDE bandwidth sensitivity (county, inhom E/F) ---\n")
if (RUN_SENSITIVITY) {
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
  sigma_local <- bw_mult * bw_diggle
  cat(sprintf("  [BW %s] start (sigma=%.3f km)\n", bw_label, sigma_local))

  lambda_local <- density(X_bg, sigma = sigma_local, edge = TRUE, at = "pixels")
  min_nz_local <- min(lambda_local$v[lambda_local$v > 0], na.rm = TRUE)
  lambda_local$v[lambda_local$v <= 0] <- min_nz_local

  bg_data <- build_background_weighted_data(lambda_local, control_ss, treated_ss)
  pp_post_bg_local <- bg_data$pp_post_bg
  pp_all_bg_local <- bg_data$pp_all_bg

  biv_init_local <- apply_pre_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
  fitE_local <- tryCatch({
    fit_etas_bivariate(
      params_init = biv_init_local, realiz = pp_all_bg_local,
      windowT = windowT_fit, windowS = win_km, m0 = ETAS_M0,
      control_state_space = control_ss, treated_state_space = treated_ss,
      background_rate_var = "W",
      treated_background_zero_before = 0,
      maxit = VANILLA_MAXIT, fixed_params = SENSITIVITY_FIXED_PARAMS, trace = 0,
      t_trunc = SEM_T_TRUNC_DAYS
    )
  }, error = function(e) {
    cat(sprintf("  [BW %s] Fit C error: %s\n", bw_label, e$message))
    NULL
  })
  E_params_local <- if (!is.null(fitE_local)) fitE_local$par else biv_init_local
  E_loglik_local <- if (!is.null(fitE_local)) fitE_local$value else NA_real_

  semF_local <- tryCatch({
    run_sem_fit(
      pp_data_in = pp_all_bg_local,
      partition_in = partition,
      partition_processes_in = partition_processes,
      state_spaces_in = state_spaces,
      init_params_in = biv_init_local,
      fixed_params_in = SENSITIVITY_FIXED_PARAMS,
      background_rate_var_in = "W",
      # Match main Fit F (not SENS_SEM_INNER_ITER) so digglex2 reproduces F when sigma matches.
      sem_inner_iter_in = SEM_INNER_ITER,
      verbose_in = FALSE,
      label = sprintf("BW %s Fit D", bw_label),
      sem_rng_label_in = OK_BW_SEM_RNG_LABEL
    )
  }, error = function(e) {
    cat(sprintf("  [BW %s] Fit D error: %s\n", bw_label, e$message))
    NULL
  })
  F_params_local <- if (!is.null(semF_local)) semF_local$etas_bivariate_params else biv_init_local

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
if (QUICK_CHECK) {
  grid_multipliers <- c(1.0)
  grid_diameters <- c(NA_real_)
  grid_max_tiles <- c(2000L)
  cat("  QUICK_CHECK mode: using one coarse grid (~100 tiles)\n")
} else {
  grid_multipliers <- c(1.0, 2.0, 5.0)
  grid_diameters <- grid_multipliers * trigger_range_km
  # Use distinct tile caps so the three radii produce meaningfully different grids.
  grid_max_tiles <- if (TEST_MODE) c(2000L, 1600L, 1200L) else c(5000L, 3000L, 1500L)
  cat(sprintf("  Triggering range = %.2f km; grid diameters: %s km\n",
              trigger_range_km,
              paste(round(grid_diameters, 2), collapse = ", ")))
}

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
  grid_part <- quadrats(win, nx = nx, ny = ny)

  n_tiles <- grid_part$n
  grid_procs <- character(n_tiles)
  tile_list <- tiles(grid_part)
  for (i in seq_len(n_tiles)) {
    tile_cent <- c(mean(tile_list[[i]]$xrange), mean(tile_list[[i]]$yrange))
    grid_procs[i] <- if (inside.owin(tile_cent[1], tile_cent[2], aoi_owin)) {
      "treated"
    } else {
      "control"
    }
  }
  names(grid_procs) <- tilenames(grid_part)

  g_treated_idx <- grid_procs == "treated"
  g_ctrl_ss <- tryCatch(as.owin(grid_part[!g_treated_idx]), error = function(e) control_ss)
  g_treat_ss <- tryCatch(as.owin(grid_part[g_treated_idx]), error = function(e) treated_ss)
  g_state_spaces <- list(control = g_ctrl_ss, treated = g_treat_ss)

  cat(sprintf("  [%s] %d tiles (%d treated, %d control)\n",
              label, n_tiles, sum(g_treated_idx), sum(!g_treated_idx)))

  list(partition = grid_part, processes = grid_procs,
       state_spaces = g_state_spaces, label = label,
       treated_idx = g_treated_idx)
}

aoi_owin <- tryCatch({
  # Use the precomputed AOI union and avoid an additional costly st_union call.
  aoi_coords <- st_coordinates(aoi_union)
  x_km <- aoi_coords[, 1] / 1000
  y_km <- aoi_coords[, 2] / 1000
  tryCatch(
    owin(poly = list(x = rev(x_km), y = rev(y_km))),
    error = function(e) owin(poly = list(x = x_km, y = y_km))
  )
}, error = function(e) {
  cat("  AOI owin construction failed, using treated_ss as fallback\n")
  treated_ss
})

grid_partitions <- lapply(seq_along(grid_diameters), function(i) {
  build_grid_partition(
    diameter = grid_diameters[i],
    win = win_km,
    aoi_owin = aoi_owin,
    label = if (QUICK_CHECK) "grid_coarse" else sprintf("grid_%.1fR", grid_multipliers[i]),
    max_tiles = grid_max_tiles[i]
  )
})

cat("\n  Building AOI partition (AOI = treated, rest = control)...\n")
aoi_treat_ss <- intersect.owin(aoi_owin, win_km)
aoi_ctrl_ss <- setminus.owin(win_km, aoi_treat_ss)
aoi_partition <- tess(tiles = list(control = aoi_ctrl_ss, treated = aoi_treat_ss),
                      window = win_km)
aoi_procs <- c(control = "control", treated = "treated")
aoi_state_spaces <- list(control = aoi_ctrl_ss, treated = aoi_treat_ss)
aoi_treated_idx <- c(FALSE, TRUE)

aoi_scheme <- list(partition = aoi_partition, processes = aoi_procs,
                   state_spaces = aoi_state_spaces, label = "aoi_region",
                   treated_idx = aoi_treated_idx)

all_partitions <- c(
  list(county = list(partition = partition, processes = partition_processes,
                     state_spaces = state_spaces, label = "county",
                     treated_idx = treated_idx)),
  setNames(grid_partitions, sapply(grid_partitions, `[[`, "label")),
  list(aoi_region = aoi_scheme)
)

if (!RUN_SENSITIVITY) {
  all_partitions <- list()
}
cat(sprintf("  Total partition schemes scheduled: %d\n", length(all_partitions)))

analysis_plots$partition_map_ids <- vapply(all_partitions, `[[`, character(1), "label")
cat("  Partition map files are disabled; render partition map panels from results payload.\n")

# ============================================================================
# 5b. Joint sensitivity dispatch: bandwidth + partition fits (single parallel layer)
# ============================================================================
cat("\n--- Step 5b: Joint sensitivity dispatch (bandwidth + partition) ---\n")

assign_to_partition <- function(df, part_info) {
  ti <- as.integer(tileindex(df$x, df$y, part_info$partition))
  df$location_process <- ifelse(is.na(ti), NA_character_,
                                part_info$processes[pmin(pmax(ti, 1),
                                  part_info$partition$n)])
  df$process <- df$location_process
  df$inferred_process <- df$location_process
  df$W <- 1.0
  df$background <- TRUE
  df[!is.na(df$location_process), ]
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
  pp_pre_p_bg <- rbind(
    normalize_bg_weights(pp_pre_p[pp_pre_p$location_process == "control", ],
                         p_ctrl_ss, lambda_im)$new_df,
    normalize_bg_weights(pp_pre_p[pp_pre_p$location_process == "treated", ],
                         p_treat_ss, lambda_im)$new_df
  )
  pp_pre_p_bg <- pp_pre_p_bg[order(pp_pre_p_bg$t), ]
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
      maxit = VANILLA_MAXIT, fixed_params = SENSITIVITY_FIXED_PARAMS, trace = 0,
      t_trunc = SEM_T_TRUNC_DAYS
    )
  }, error = function(e) { cat(sprintf("    [%s] Inhom naive fit error: %s\n", label, e$message)); NULL })

  E_params_p <- if (!is.null(fitE_p)) fitE_p$par else biv_init_p

  # SEM bivariate + KDE background
  biv_init_sem_p <- biv_init_p
  semF_p <- tryCatch({
    run_sem_fit(
      pp_data_in = pp_all_p,
      partition_in = part_info$partition,
      partition_processes_in = part_info$processes,
      state_spaces_in = p_state_spaces,
      init_params_in = biv_init_sem_p,
      fixed_params_in = SENSITIVITY_FIXED_PARAMS,
      background_rate_var_in = "W",
      sem_inner_iter_in = SENS_SEM_INNER_ITER,
      verbose_in = FALSE,
      label = sprintf("Partition %s Fit D", label)
    )
  }, error = function(e) { cat(sprintf("    [%s] Inhom SEM error: %s\n", label, e$message)); NULL })

  F_params_p <- if (!is.null(semF_p)) semF_p$etas_bivariate_params else biv_init_sem_p
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

  cat(sprintf("    [%s] E params: %s\n", label,
              paste(biv_names, round(E_params_p, 4), sep = "=", collapse = ", ")))
  cat(sprintf("    [%s] F params: %s\n", label,
              paste(biv_names, round(F_params_p, 4), sep = "=", collapse = ", ")))

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
cat(sprintf("  Sensitivity jobs: %d bandwidth + %d partition = %d total\n",
            length(kde_bandwidth_specs), length(all_partitions), length(sensitivity_jobs)))

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

  kde_bandwidth_fits_local <- lapply(vapply(kde_bandwidth_specs, `[[`, character(1), "label"), function(lbl) {
    idx <- which(vapply(sens_out, function(z) identical(z$type, "bandwidth") && identical(z$id, lbl), logical(1)))
    if (length(idx) == 0) return(NULL)
    sens_out[[idx[1]]]$out
  })
  names(kde_bandwidth_fits_local) <- vapply(kde_bandwidth_specs, `[[`, character(1), "label")

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

windowT_ate <- c(0, ATE_WINDOW_DAYS)
cat(sprintf("  ATE window: [0, %d] days\n", ATE_WINDOW_DAYS))

extract_marginals <- function(biv_par) {
  if (is.null(biv_par)) return(NULL)
  c_val <- biv_par[["c"]]; if (!is.finite(c_val)) c_val <- STRUCT_DEFAULTS$c
  p_val <- biv_par[["p"]]; if (!is.finite(p_val)) p_val <- STRUCT_DEFAULTS$p
  D_val <- biv_par[["D"]]; if (!is.finite(D_val)) D_val <- STRUCT_DEFAULTS$D
  g_val <- biv_par[["gamma"]]; if (!is.finite(g_val)) g_val <- STRUCT_DEFAULTS$gamma
  q_val <- biv_par[["q"]]; if (!is.finite(q_val)) q_val <- STRUCT_DEFAULTS$q
  structural <- c(c = c_val, p = p_val, D = D_val, gamma = g_val, q = q_val)
  ctrl <- as.list(c(mu = biv_par[["mu_0"]], A = biv_par[["A_00"]],
    alpha_m = biv_par[["alpha_m_00"]], structural))
  treat <- as.list(c(mu = biv_par[["mu_1"]], A = biv_par[["A_11"]],
    alpha_m = biv_par[["alpha_m_11"]], structural))
  list(ctrl = ctrl, treat = treat)
}

extract_marginals_indep <- function(indep_par) {
  if (is.null(indep_par) || is.null(indep_par$control) || is.null(indep_par$treated)) return(NULL)
  list(ctrl = indep_par$control, treat = indep_par$treated)
}

B_marginals <- extract_marginals(B_params)
D_marginals <- extract_marginals(D_params)
E_marginals <- extract_marginals(E_params)
F_marginals <- extract_marginals(F_params)
G_marginals <- extract_marginals(G_params)
H_marginals <- extract_marginals(H_params)
I_marginals <- extract_marginals_indep(I_params)
J_marginals <- extract_marginals_indep(J_params)

ate_estim_fast <- function(ctrl_pp, treat_pp, observed_data, label,
                           n_tiles_used = partition$n,
                           treated_idx_used = treated_idx,
                           filtration_history = NULL,
                           crn_base_seed = NA_integer_,
                           phase = "main_fit",
                           quiet = FALSE) {
  phase_tag <- tolower(gsub("[^A-Za-z0-9]+", "_", as.character(phase)))
  if (!nzchar(phase_tag)) phase_tag <- "main_fit"
  phase_label <- switch(
    phase_tag,
    main_fit = "main-fit",
    sensitivity = "sensitivity",
    bootstrap = "bootstrap",
    phase_tag
  )
  ate_prefix <- sprintf("[ATE:%s]", phase_label)
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
    ate_approx_branching_ratio <- function(par_obj, beta_gr = BETA_GR) {
      if (is.null(par_obj)) return(NA_real_)
      A_val <- suppressWarnings(as.numeric(par_obj$A))
      alpha_val <- suppressWarnings(as.numeric(par_obj$alpha_m))
      if (!is.finite(A_val) || !is.finite(alpha_val)) return(NA_real_)
      A_val * exp(alpha_val / beta_gr)
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
    ctrl_br <- ate_approx_branching_ratio(ctrl_pp)
    treat_br <- ate_approx_branching_ratio(treat_pp)
    ctrl_diag <- summarize_process_stability(ctrl_pp)
    treat_diag <- summarize_process_stability(treat_pp)
    if (!quiet) {
      cat(sprintf("    %s sim config: sims=%d cores=%d backend=%s branch_proxy(ctrl=%.3f, treat=%.3f)\n",
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
    if ((is.finite(ctrl_br) && ctrl_br > BOOT_BRANCHING_MAX) ||
        (is.finite(treat_br) && treat_br > BOOT_BRANCHING_MAX)) {
      cat(sprintf("    [warning:ate-explosive] %s: branching proxy exceeds %.3f (ctrl=%.3f, treat=%.3f)\n",
                  label, BOOT_BRANCHING_MAX, ctrl_br, treat_br))
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
    run_one_sim <- function(s) {
      s_int <- suppressWarnings(as.integer(s))
      if (!is.finite(s_int) || is.na(s_int)) s_int <- 1L
      if (isTRUE(OK_ATE_USE_CRN)) {
        seed_s <- as.integer(crn_base_seed + s_int)
        if (isTRUE(OK_ATE_CRN_PAIR)) {
          set.seed(seed_s)
          c_sim <- sim_etas(ctrl_pp, windowT_ate, windowS = win_km,
                            m0 = ETAS_M0, beta_gr = BETA_GR,
                            filtration = pre_history)
          set.seed(seed_s)
          t_sim <- sim_etas(treat_pp, windowT_ate, windowS = win_km,
                            m0 = ETAS_M0, beta_gr = BETA_GR,
                            filtration = pre_history)
        } else {
          set.seed(seed_s)
          c_sim <- sim_etas(ctrl_pp, windowT_ate, windowS = win_km,
                            m0 = ETAS_M0, beta_gr = BETA_GR,
                            filtration = pre_history)
          set.seed(as.integer(seed_s + 1000000L))
          t_sim <- sim_etas(treat_pp, windowT_ate, windowS = win_km,
                            m0 = ETAS_M0, beta_gr = BETA_GR,
                            filtration = pre_history)
        }
      } else {
        c_sim <- sim_etas(ctrl_pp, windowT_ate, windowS = win_km,
                          m0 = ETAS_M0, beta_gr = BETA_GR,
                          filtration = pre_history)
        t_sim <- sim_etas(treat_pp, windowT_ate, windowS = win_km,
                          m0 = ETAS_M0, beta_gr = BETA_GR,
                          filtration = pre_history)
      }
      c(c_count = length(c_sim$t), t_count = length(t_sim$t))
    }
    sim_results <- if (N_CORES > 1 && ATE_N_SIMS > 1) {
      if (!is.null(ate_cl_reuse)) {
        run_parallel_on_cluster(
          ate_cl_reuse,
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
    # Estimand of interest: earthquakes saved by treatment regime versus control-everywhere.
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
         analytic_saved = analytic_saved)
  }, error = function(e) { cat("    Error:", e$message, "\n"); NULL })
}

ate_cl_reuse <- NULL
ate_cl_reuse_cores <- max(1L, min(ATE_SIM_CORES, ATE_N_SIMS))
if (PARALLEL_BACKEND == "psock" && N_CORES > 1 && ATE_N_SIMS > 1 && ate_cl_reuse_cores > 1L) {
  cat(sprintf("  [parallel:ate-sim] initializing reusable PSOCK pool: cores=%d\n", ate_cl_reuse_cores))
  ate_cl_reuse <- parallel::makePSOCKcluster(ate_cl_reuse_cores, outfile = "")
  on.exit({
    if (!is.null(ate_cl_reuse)) {
      try(parallel::stopCluster(ate_cl_reuse), silent = TRUE)
      ate_cl_reuse <<- NULL
    }
  }, add = TRUE)
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

t_ate_B <- proc.time()[["elapsed"]]
ate_B <- ate_estim_fast(B_marginals$ctrl, B_marginals$treat, pp_post,
                        "Fit A (naive bivariate)")
ate_B_elapsed <- proc.time()[["elapsed"]] - t_ate_B
add_timing_row("ate_B", ate_B_elapsed, if (!is.null(ate_B)) "ok" else "failed")
t_ate_D <- proc.time()[["elapsed"]]
ate_D <- ate_estim_fast(D_marginals$ctrl, D_marginals$treat, pp_post_sem_D,
                        "Fit B (SEM bivariate)")
ate_D_elapsed <- proc.time()[["elapsed"]] - t_ate_D
add_timing_row("ate_D", ate_D_elapsed, if (!is.null(ate_D)) "ok" else "failed")
t_ate_E <- proc.time()[["elapsed"]]
ate_E <- ate_estim_fast(E_marginals$ctrl, E_marginals$treat, pp_post_bg,
                        "Fit C (naive biv+KDE, control-only fixed)")
ate_E_elapsed <- proc.time()[["elapsed"]] - t_ate_E
add_timing_row("ate_E", ate_E_elapsed, if (!is.null(ate_E)) "ok" else "failed")
t_ate_F <- proc.time()[["elapsed"]]
ate_F <- ate_estim_fast(F_marginals$ctrl, F_marginals$treat, pp_post_sem_F,
                        "Fit D (SEM biv+KDE, control-only fixed)")
ate_F_elapsed <- proc.time()[["elapsed"]] - t_ate_F
add_timing_row("ate_F", ate_F_elapsed, if (!is.null(ate_F)) "ok" else "failed")

# Compute ATE/saved estimands for all KDE fix-profile pairs (C/D, E/F, G/H).
for (vid in names(kde_variant_specs)) {
  if (!is.null(kde_variant_fits$E[[vid]])) {
    p_e <- kde_variant_fits$E[[vid]]$params
    m_e <- extract_marginals(p_e)
    letter_e <- KDE_FIT_LETTERS[[vid]]$E
    kde_variant_fits$E[[vid]]$marginals <- m_e
    kde_variant_fits$E[[vid]]$ate <- ate_estim_fast(
      m_e$ctrl, m_e$treat, pp_post_bg,
      label = sprintf("Fit %s (naive biv+KDE, %s)", letter_e, vid),
      quiet = FALSE
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
    kde_variant_fits$F[[vid]]$ate <- ate_estim_fast(
      m_f$ctrl, m_f$treat, post_sem_local,
      label = sprintf("Fit %s (SEM biv+KDE, %s)", letter_f, vid),
      quiet = FALSE
    )
  }
}
ate_I <- if (!is.null(I_marginals)) {
  t_ate_I <- proc.time()[["elapsed"]]
  ate_estim_fast(I_marginals$ctrl, I_marginals$treat, pp_post_bg,
                 "Fit I (Naive univariate ETAS)")
} else NULL
if (exists("t_ate_I", inherits = FALSE)) {
  ate_I_elapsed <- proc.time()[["elapsed"]] - t_ate_I
  add_timing_row("ate_I", ate_I_elapsed, if (!is.null(ate_I)) "ok" else "failed")
}
ate_J <- if (!is.null(J_marginals)) {
  t_ate_J <- proc.time()[["elapsed"]]
  ate_estim_fast(J_marginals$ctrl, J_marginals$treat, pp_post_sem_J,
                 "Fit J (SEM univariate ETAS)")
} else NULL
if (exists("t_ate_J", inherits = FALSE)) {
  ate_J_elapsed <- proc.time()[["elapsed"]] - t_ate_J
  add_timing_row("ate_J", ate_J_elapsed, if (!is.null(ate_J)) "ok" else "failed")
}

# Save checkpoint with all main-fit ATE payloads, but without sensitivity payloads.
cat("\n--- Step 6a checkpoint: saving pre-sensitivity results ---\n")
pre_sens_saved_at <- as.character(Sys.time())
fits_named_pre_sensitivity <- list(
  A = list(
    letter = "A", label = "Naive bivariate (county)",
    method = "county_bivariate", algorithm = "naive",
    params = B_params, fit = fitB, ate = ate_B
  ),
  B = list(
    letter = "B", label = "SEM bivariate (county)",
    method = "county_bivariate", algorithm = "sem",
    params = D_params, fit = semD, ctrl = D_ctrl, treat = D_treat, ate = ate_D
  ),
  C = c(list(letter = "C", label = "Naive biv+KDE (control_only_fixed)",
             method = "kde_control_only_fixed", algorithm = "naive"),
        kde_variant_fits$E$control_only_fixed),
  D = c(list(letter = "D", label = "SEM biv+KDE (control_only_fixed)",
             method = "kde_control_only_fixed", algorithm = "sem"),
        kde_variant_fits$F$control_only_fixed),
  E = c(list(letter = "E", label = "Naive biv+KDE (all_free)",
             method = "kde_all_free", algorithm = "naive"),
        kde_variant_fits$E$all_free),
  F = c(list(letter = "F", label = "SEM biv+KDE (all_free)",
             method = "kde_all_free", algorithm = "sem"),
        kde_variant_fits$F$all_free),
  G = c(list(letter = "G", label = "Naive biv+KDE (productivity_free)",
             method = "kde_productivity_free", algorithm = "naive"),
        kde_variant_fits$E$productivity_free),
  H = c(list(letter = "H", label = "SEM biv+KDE (productivity_free)",
             method = "kde_productivity_free", algorithm = "sem"),
        kde_variant_fits$F$productivity_free),
  I = list(
    letter = "I", label = "Naive univariate ETAS + KDE (county)",
    method = "county_univariate_kde", algorithm = "naive",
    params = I_params, fit = fitI, ate = ate_I
  ),
  J = list(
    letter = "J", label = "SEM univariate ETAS + KDE (county)",
    method = "county_univariate_kde", algorithm = "sem",
    params = J_params, fit = semJ, ctrl = J_ctrl, treat = J_treat, ate = ate_J
  )
)
results_pre_sensitivity <- list(
  fitB = list(params = B_params, loglik = B_loglik, fit = fitB, ate = ate_B),
  fitD = list(params = D_params, ctrl = D_ctrl, treat = D_treat, sem = semD, ate = ate_D),
  fitE = list(params = E_params, loglik = E_loglik, fit = fitE, ate = ate_E),
  fitF = list(params = F_params, ctrl = F_ctrl, treat = F_treat, sem = semF, ate = ate_F),
  fitI = list(params = I_params, fit = fitI, ate = ate_I),
  fitJ = list(params = J_params, fit = semJ, ctrl = J_ctrl, treat = J_treat, ate = ate_J),
  fits_named = fits_named_pre_sensitivity,
  kde_variants = kde_variant_fits,
  fitG = NULL,
  fitH = NULL,
  bootstrap_ate = NULL,
  partition_results = NULL,
  ate_partitions = NULL,
  kde_bandwidth_sensitivity = NULL,
  pp_data = list(pp_pre = pp_pre, pp_pre_holdout = pp_pre_holdout, pp_post = pp_post),
  kde_info = list(
    bw_sigma = as.numeric(bw_sigma), n_training = n_pre_holdout_ctrl,
    n_pre_ctrl_holdout = n_pre_holdout_ctrl,
    n_pre_holdout = nrow(pp_pre_holdout),
    n_pre_used = nrow(pp_pre),
    trigger_range_km = trigger_range_km
  ),
  control_snapshot_fits = control_snapshot_fits,
  sem_pilot_tuning = sem_pilot_tuning,
  plots = analysis_plots,
  counties = list(
    names = counties_sf_valid$NAME,
    treated_names = treated_names,
    n_counties = partition$n,
    n_treated = sum(treated_idx)
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
    FIXED_STRUCTURAL = FIXED_STRUCTURAL,
    SEM_N_ITER = SEM_N_ITER, SEM_INNER_ITER = SEM_INNER_ITER,
    SEM_INNER_PROPS = SEM_INNER_PROPS,
    SEM_N_LABELLINGS = SEM_N_LABELLINGS,
    SEM_CHANGE_FACTOR = SEM_CHANGE_FACTOR,
    SEM_CHANGE_FACTOR_MIN_MULT = SEM_CHANGE_FACTOR_MIN_MULT,
    SEM_CHANGE_FACTOR_MAX_MULT = SEM_CHANGE_FACTOR_MAX_MULT,
    SEM_OPTIM_METHOD = SEM_OPTIM_METHOD,
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
    RUN_SENSITIVITY = RUN_SENSITIVITY,
    KDE_VARIANT_MODE = KDE_VARIANT_MODE,
    RUN_KDE_PROFILE_SWEEP = RUN_KDE_PROFILE_SWEEP,
    MEMORY_SAFE = MEMORY_SAFE,
    TRIM_SENS_OBJECTS = TRIM_SENS_OBJECTS,
    SENS_CORES = SENS_CORES,
    ATE_SIM_CORES = ATE_SIM_CORES,
    ATE_USE_CRN = OK_ATE_USE_CRN,
    ATE_CRN_PAIR = OK_ATE_CRN_PAIR,
    ATE_CONDITIONAL_ON_PRE = OK_ATE_CONDITIONAL_ON_PRE,
    ATE_CRN_BASE = OK_ATE_CRN_BASE,
    SEM_MAX_RELABEL_STEP_FRAC = SEM_MAX_RELABEL_STEP_FRAC,
    SEM_FORCE_PARAM_UPDATE_FLIP_FRAC = SEM_FORCE_PARAM_UPDATE_FLIP_FRAC,
    ATE_N_SIMS = ATE_N_SIMS, ATE_WINDOW_DAYS = ATE_WINDOW_DAYS,
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
    BOOT_BRANCHING_MAX = BOOT_BRANCHING_MAX,
    BOOT_MAX_PRE_EVENTS = BOOT_MAX_PRE_EVENTS,
    BOOT_MAX_POST_EVENTS_PER_PROC = BOOT_MAX_POST_EVENTS_PER_PROC,
    BOOT_MAX_TOTAL_EVENTS = BOOT_MAX_TOTAL_EVENTS,
    TEST_MODE = TEST_MODE,
    windowT_post = windowT_post,
    n_pre = nrow(pp_pre), n_pre_holdout = nrow(pp_pre_holdout), n_pre_total = nrow(pp_pre_all),
    n_post = nrow(pp_post),
    n_counties = partition$n, n_treated = sum(treated_idx),
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

# ATE sensitivity by KDE bandwidth (county only, inhomogeneous E/F)
kde_bandwidth_sensitivity <- lapply(kde_bandwidth_fits, function(kf) {
  if (is.null(kf) || is.null(kf$E_params) || is.null(kf$F_params)) return(NULL)
  em <- extract_marginals(kf$E_params)
  fm <- extract_marginals(kf$F_params)
  list(
    label = kf$label,
    multiplier = kf$multiplier,
    sigma = kf$sigma,
    E_params = kf$E_params,
    F_params = kf$F_params,
    n_relabel = kf$n_relabel,
    ate_E = tryCatch(
      ate_estim_fast(em$ctrl, em$treat, kf$pp_post_bg,
                     sprintf("BW %s naive biv+KDE", kf$label),
                     phase = "sensitivity",
                     n_tiles_used = partition$n,
                     treated_idx_used = treated_idx),
      error = function(e) NULL
    ),
    ate_F = tryCatch(
      ate_estim_fast(fm$ctrl, fm$treat, kf$pp_post_sem,
                     sprintf("BW %s SEM biv+KDE", kf$label),
                     phase = "sensitivity",
                     n_tiles_used = partition$n,
                     treated_idx_used = treated_idx),
      error = function(e) NULL
    )
  )
})

# ATE for alternative partitions
ate_partitions <- lapply(partition_results, function(pr) {
  if (pr$label == "county") return(NULL)
  em <- extract_marginals(pr$E_params)
  fm <- extract_marginals(pr$F_params)
  part_info <- all_partitions[[pr$label]]
  pp_post_part <- assign_to_partition(pp_post[, c("x", "y", "t", "mag")], part_info)
  pp_post_part_bg <- rbind(
    normalize_bg_weights(pp_post_part[pp_post_part$location_process == "control", ],
                         part_info$state_spaces$control, lambda_im)$new_df,
    normalize_bg_weights(pp_post_part[pp_post_part$location_process == "treated", ],
                         part_info$state_spaces$treated, lambda_im)$new_df
  )
  pp_post_part_bg <- pp_post_part_bg[order(pp_post_part_bg$t), ]
  pp_post_sem_part <- if (!is.null(pr$pp_post_sem)) {
    pr$pp_post_sem
  } else if (!is.null(pr$semF) && !is.null(pr$semF$adaptive$adaptive_labelling)) {
    pr$semF$adaptive$adaptive_labelling
  } else {
    pp_post_part_bg
  }
  pp_post_sem_part <- pp_post_sem_part[pp_post_sem_part$t >= 0, ]
  list(
    ate_E = tryCatch(
      ate_estim_fast(em$ctrl, em$treat, pp_post_part_bg,
                     sprintf("%s naive biv+KDE", pr$label),
                     phase = "sensitivity",
                     n_tiles_used = pr$n_tiles,
                     treated_idx_used = part_info$treated_idx),
      error = function(e) NULL),
    ate_F = tryCatch(
      ate_estim_fast(fm$ctrl, fm$treat, pp_post_sem_part,
                     sprintf("%s SEM biv+KDE", pr$label),
                     phase = "sensitivity",
                     n_tiles_used = pr$n_tiles,
                     treated_idx_used = part_info$treated_idx),
      error = function(e) NULL)
  )
})

# Save a checkpoint before bootstrap so long runs retain core fit outputs
# even if bootstrap gets interrupted or OOM-killed.
cat("\n--- Step 7b checkpoint: saving pre-bootstrap results ---\n")
# Include report-facing fields (pp_data, counties, kde_info, full config) so Quarto can
# render from this checkpoint alone — same keys as final oklahoma_results.rds.
pre_boot_saved_at <- as.character(Sys.time())

fits_named <- list(
  A = list(
    letter = "A", label = "Naive bivariate (county)",
    method = "county_bivariate", algorithm = "naive",
    params = B_params, fit = fitB, ate = ate_B
  ),
  B = list(
    letter = "B", label = "SEM bivariate (county)",
    method = "county_bivariate", algorithm = "sem",
    params = D_params, fit = semD, ctrl = D_ctrl, treat = D_treat, ate = ate_D
  ),
  C = c(list(letter = "C", label = "Naive biv+KDE (control_only_fixed)",
             method = "kde_control_only_fixed", algorithm = "naive"),
        kde_variant_fits$E$control_only_fixed),
  D = c(list(letter = "D", label = "SEM biv+KDE (control_only_fixed)",
             method = "kde_control_only_fixed", algorithm = "sem"),
        kde_variant_fits$F$control_only_fixed),
  E = c(list(letter = "E", label = "Naive biv+KDE (all_free)",
             method = "kde_all_free", algorithm = "naive"),
        kde_variant_fits$E$all_free),
  F = c(list(letter = "F", label = "SEM biv+KDE (all_free)",
             method = "kde_all_free", algorithm = "sem"),
        kde_variant_fits$F$all_free),
  G = c(list(letter = "G", label = "Naive biv+KDE (productivity_free)",
             method = "kde_productivity_free", algorithm = "naive"),
        kde_variant_fits$E$productivity_free),
  H = c(list(letter = "H", label = "SEM biv+KDE (productivity_free)",
             method = "kde_productivity_free", algorithm = "sem"),
        kde_variant_fits$F$productivity_free),
  I = list(
    letter = "I", label = "Naive univariate ETAS + KDE (county)",
    method = "county_univariate_kde", algorithm = "naive",
    params = I_params, fit = fitI, ate = ate_I
  ),
  J = list(
    letter = "J", label = "SEM univariate ETAS + KDE (county)",
    method = "county_univariate_kde", algorithm = "sem",
    params = J_params, fit = semJ, ctrl = J_ctrl, treat = J_treat, ate = ate_J
  )
)

results_pre_bootstrap <- list(
  fitB = list(params = B_params, loglik = B_loglik, fit = fitB, ate = ate_B),
  fitD = list(params = D_params, ctrl = D_ctrl, treat = D_treat, sem = semD, ate = ate_D),
  fitE = list(params = E_params, loglik = E_loglik, fit = fitE, ate = ate_E),
  fitF = list(params = F_params, ctrl = F_ctrl, treat = F_treat, sem = semF, ate = ate_F),
  fitI = list(params = I_params, fit = fitI, ate = ate_I),
  fitJ = list(params = J_params, fit = semJ, ctrl = J_ctrl, treat = J_treat, ate = ate_J),
  fits_named = fits_named,
  kde_variants = kde_variant_fits,
  fitG = NULL,
  fitH = NULL,
  bootstrap_ate = NULL,
  partition_results = partition_results,
  ate_partitions = ate_partitions,
  kde_bandwidth_sensitivity = kde_bandwidth_sensitivity,
  pp_data = list(pp_pre = pp_pre, pp_pre_holdout = pp_pre_holdout, pp_post = pp_post),
  kde_info = list(
    bw_sigma = as.numeric(bw_sigma), n_training = n_pre_holdout_ctrl,
    n_pre_ctrl_holdout = n_pre_holdout_ctrl,
    n_pre_holdout = nrow(pp_pre_holdout),
    n_pre_used = nrow(pp_pre),
    trigger_range_km = trigger_range_km
  ),
  control_snapshot_fits = control_snapshot_fits,
  sem_pilot_tuning = sem_pilot_tuning,
  plots = analysis_plots,
  counties = list(
    names = counties_sf_valid$NAME,
    treated_names = treated_names,
    n_counties = partition$n,
    n_treated = sum(treated_idx)
  ),
  metadata = list(stage = "pre_bootstrap", saved_at = pre_boot_saved_at),
  checkpoint = list(
    stage = "pre_bootstrap",
    saved_at = pre_boot_saved_at,
    boot_targets_requested = BOOT_TARGETS,
    boot_targets_run = intersect(BOOT_TARGETS, c("E", "F")),
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
    FIXED_STRUCTURAL = FIXED_STRUCTURAL,
    SEM_N_ITER = SEM_N_ITER, SEM_INNER_ITER = SEM_INNER_ITER,
    SEM_INNER_PROPS = SEM_INNER_PROPS,
    SEM_N_LABELLINGS = SEM_N_LABELLINGS,
    SEM_CHANGE_FACTOR = SEM_CHANGE_FACTOR,
    SEM_CHANGE_FACTOR_MIN_MULT = SEM_CHANGE_FACTOR_MIN_MULT,
    SEM_CHANGE_FACTOR_MAX_MULT = SEM_CHANGE_FACTOR_MAX_MULT,
    SEM_OPTIM_METHOD = SEM_OPTIM_METHOD,
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
    RUN_SENSITIVITY = RUN_SENSITIVITY,
    KDE_VARIANT_MODE = KDE_VARIANT_MODE,
    RUN_KDE_PROFILE_SWEEP = RUN_KDE_PROFILE_SWEEP,
    MEMORY_SAFE = MEMORY_SAFE,
    TRIM_SENS_OBJECTS = TRIM_SENS_OBJECTS,
    SENS_CORES = SENS_CORES,
    ATE_SIM_CORES = ATE_SIM_CORES,
    ATE_USE_CRN = OK_ATE_USE_CRN,
    ATE_CRN_PAIR = OK_ATE_CRN_PAIR,
    ATE_CONDITIONAL_ON_PRE = OK_ATE_CONDITIONAL_ON_PRE,
    ATE_CRN_BASE = OK_ATE_CRN_BASE,
    SEM_MAX_RELABEL_STEP_FRAC = SEM_MAX_RELABEL_STEP_FRAC,
    SEM_FORCE_PARAM_UPDATE_FLIP_FRAC = SEM_FORCE_PARAM_UPDATE_FLIP_FRAC,
    ATE_N_SIMS = ATE_N_SIMS, ATE_WINDOW_DAYS = ATE_WINDOW_DAYS,
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
    BOOT_BRANCHING_MAX = BOOT_BRANCHING_MAX,
    BOOT_MAX_PRE_EVENTS = BOOT_MAX_PRE_EVENTS,
    BOOT_MAX_POST_EVENTS_PER_PROC = BOOT_MAX_POST_EVENTS_PER_PROC,
    BOOT_MAX_TOTAL_EVENTS = BOOT_MAX_TOTAL_EVENTS,
    TEST_MODE = TEST_MODE,
    windowT_post = windowT_post,
    n_pre = nrow(pp_pre), n_pre_holdout = nrow(pp_pre_holdout), n_pre_total = nrow(pp_pre_all),
    n_post = nrow(pp_post),
    n_counties = partition$n, n_treated = sum(treated_idx),
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

# Parametric bootstrap ATEs for the all-free KDE pair (E/F).
bootstrap_ate <- NULL
boot_targets_run <- intersect(BOOT_TARGETS, c("E", "F"))
bootstrap_elapsed <- NA_real_
if (RUN_BOOTSTRAP_ATE && BOOT_N_REPS > 0L && length(boot_targets_run) > 0L) {
  t_bootstrap <- proc.time()[["elapsed"]]
  cat(sprintf("\n--- Step 8: Parametric bootstrap ATEs (targets=%s, reps=%d, scope=%s, outer cores=%d, boot SEM inner=%d) ---\n",
              paste(boot_targets_run, collapse = ","), BOOT_N_REPS, BOOT_REFIT_SCOPE, BOOT_OUTER_CORES, BOOT_SEM_INNER_ITER))

  if (BOOT_REFIT_SCOPE == "full") {
    cat("  Full scope selected: for current targets (E/F), this runs per-replicate refits before ATE estimation.\n")
  }

  as_pp_df <- function(sim_obj, location_process_value, inferred_process_value = NULL) {
    if (is.null(sim_obj) || is.null(sim_obj$t) || length(sim_obj$t) < 1) {
      out <- pp_post[0, c("x", "y", "t", "mag", "location_process", "process", "inferred_process"), drop = FALSE]
      out$location_process <- character(0)
      out$process <- character(0)
      out$inferred_process <- character(0)
      return(out)
    }
    out <- data.frame(
      x = sim_obj$x, y = sim_obj$y, t = sim_obj$t, mag = sim_obj$mag,
      stringsAsFactors = FALSE
    )
    out$location_process <- location_process_value
    out$process <- location_process_value
    out$inferred_process <- if (is.null(inferred_process_value)) location_process_value else inferred_process_value
    out
  }

  pre_window_boot <- c(min(pp_pre$t, na.rm = TRUE), 0)
  post_window_boot <- windowT_post
  pre_ctrl_seed <- PRE_CTRL_BOOT_PARAMS
  # Bootstrap targets E/F correspond to the all-free KDE pair.
  c_ctrl_seed <- E_marginals$ctrl
  c_treat_seed <- E_marginals$treat
  d_ctrl_seed <- F_marginals$ctrl
  d_treat_seed <- F_marginals$treat

  get_num <- function(obj, nm, default = NA_real_) {
    if (is.null(obj) || is.null(obj[[nm]])) return(default)
    v <- suppressWarnings(as.numeric(obj[[nm]]))
    if (length(v) < 1L || !is.finite(v[[1]])) return(default)
    v[[1]]
  }
  approx_branching_ratio <- function(par_obj, beta_gr = BETA_GR) {
    A <- get_num(par_obj, "A", NA_real_)
    alpha_m <- get_num(par_obj, "alpha_m", 0)
    if (!is.finite(A) || A < 0) return(Inf)
    if (!is.finite(alpha_m)) alpha_m <- 0
    # Rough ETAS branching proxy under exponential magnitudes.
    A * exp(alpha_m / beta_gr)
  }
  validate_sim_obj <- function(sim_obj, label, max_events) {
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
    }
    invisible(TRUE)
  }

  simulate_boot_data <- function(ctrl_seed, treat_seed) {
    if (OK_BOOT_GUARD_DEGENERATE) {
      pre_br <- approx_branching_ratio(pre_ctrl_seed)
      ctrl_br <- approx_branching_ratio(ctrl_seed)
      treat_br <- approx_branching_ratio(treat_seed)
      if (any(!is.finite(c(pre_br, ctrl_br, treat_br))) ||
          pre_br > BOOT_BRANCHING_MAX || ctrl_br > BOOT_BRANCHING_MAX || treat_br > BOOT_BRANCHING_MAX) {
        stop(sprintf(
          "Degenerate bootstrap guard: branching proxy too high (pre=%.3f ctrl=%.3f treat=%.3f, limit=%.3f)",
          pre_br, ctrl_br, treat_br, BOOT_BRANCHING_MAX
        ))
      }
    }
    pre_sim <- sim_etas(pre_ctrl_seed, pre_window_boot, windowS = control_ss,
                        m0 = ETAS_M0, beta_gr = BETA_GR)
    validate_sim_obj(pre_sim, "bootstrap pre_sim", BOOT_MAX_PRE_EVENTS)
    pre_df <- as_pp_df(pre_sim, "control", "control")
    history_df <- pre_df[, c("x", "y", "t", "mag"), drop = FALSE]
    ctrl_post_sim <- sim_etas(ctrl_seed, post_window_boot, windowS = control_ss,
                              m0 = ETAS_M0, beta_gr = BETA_GR, filtration = history_df)
    validate_sim_obj(ctrl_post_sim, "bootstrap ctrl_post_sim", BOOT_MAX_POST_EVENTS_PER_PROC)
    treat_post_sim <- sim_etas(treat_seed, post_window_boot, windowS = treated_ss,
                               m0 = ETAS_M0, beta_gr = BETA_GR, filtration = history_df)
    validate_sim_obj(treat_post_sim, "bootstrap treat_post_sim", BOOT_MAX_POST_EVENTS_PER_PROC)
    post_ctrl_df <- as_pp_df(ctrl_post_sim, "control", "control")
    post_treat_df <- as_pp_df(treat_post_sim, "treated", "treated")
    pp_post_sim <- rbind(post_ctrl_df, post_treat_df)
    pp_post_sim <- pp_post_sim[order(pp_post_sim$t), , drop = FALSE]

    bg_ctrl <- normalize_bg_weights(post_ctrl_df, control_ss, lambda_im)$new_df
    bg_treat <- normalize_bg_weights(post_treat_df, treated_ss, lambda_im)$new_df
    pp_post_bg_sim <- rbind(bg_ctrl, bg_treat)
    pp_post_bg_sim <- pp_post_bg_sim[order(pp_post_bg_sim$t), , drop = FALSE]

    pre_bg_ctrl <- normalize_bg_weights(pre_df, control_ss, lambda_im)$new_df
    pre_bg_treat <- normalize_bg_weights(pre_df[0, , drop = FALSE], treated_ss, lambda_im)$new_df
    pp_all_bg_sim <- rbind(pre_bg_ctrl, pre_bg_treat, pp_post_bg_sim)
    pp_all_bg_sim <- pp_all_bg_sim[order(pp_all_bg_sim$t), , drop = FALSE]
    if (nrow(pp_all_bg_sim) > BOOT_MAX_TOTAL_EVENTS) {
      stop(sprintf("bootstrap total simulated events %d exceeded cap %d", nrow(pp_all_bg_sim), BOOT_MAX_TOTAL_EVENTS))
    }

    list(
      pre_df = pre_df,
      post_ctrl_df = post_ctrl_df,
      post_treat_df = post_treat_df,
      pp_post_bg_sim = pp_post_bg_sim,
      pp_all_bg_sim = pp_all_bg_sim
    )
  }

  summarize_boot <- function(ate_obj, rep_id, pre_df, post_ctrl_df, post_treat_df, par_obj = NULL) {
    if (is.null(ate_obj) || is.null(ate_obj$all_nothing_sim)) {
      return(list(ok = FALSE, rep = rep_id, msg = "ATE estimation failed"))
    }
    list(
      ok = TRUE,
      rep = rep_id,
      n_pre_sim = nrow(pre_df),
      n_post_ctrl_sim = nrow(post_ctrl_df),
      n_post_treat_sim = nrow(post_treat_df),
      ate_total_mean = mean(ate_obj$all_nothing_sim$total_saved, na.rm = TRUE),
      ate_total_sd = stats::sd(ate_obj$all_nothing_sim$total_saved, na.rm = TRUE),
      ate_tile_mean = mean(ate_obj$all_nothing_sim$ATE, na.rm = TRUE),
      eta_ctrl = ate_obj$analytic$eta_ctrl,
      eta_treat = ate_obj$analytic$eta_treat,
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

    if ("E" %in% boot_targets_run) {
      out$E <- tryCatch({
        sim_C <- simulate_boot_data(c_ctrl_seed, c_treat_seed)
        c_params_boot <- E_params
        if (BOOT_REFIT_SCOPE %in% c("partial", "full")) {
          fit_c_boot <- tryCatch({
            fit_etas_bivariate(
              params_init = E_params, realiz = sim_C$pp_all_bg_sim,
              windowT = windowT_fit, windowS = win_km, m0 = ETAS_M0,
              control_state_space = control_ss, treated_state_space = treated_ss,
              background_rate_var = "W",
              treated_background_zero_before = 0,
              maxit = VANILLA_MAXIT, fixed_params = SENSITIVITY_FIXED_PARAMS, trace = 0,
              t_trunc = SEM_T_TRUNC_DAYS
            )
          }, error = function(e) NULL)
          if (!is.null(fit_c_boot) && !is.null(fit_c_boot$par)) c_params_boot <- fit_c_boot$par
        }
        c_marg_boot <- extract_marginals(c_params_boot)
        ate_c_boot <- ate_estim_fast(
          c_marg_boot$ctrl, c_marg_boot$treat, sim_C$pp_post_bg_sim,
          label = sprintf("Boot E #%d", rep_id),
          filtration_history = sim_C$pre_df,
          crn_base_seed = boot_ate_crn_seed,
          phase = "bootstrap",
          n_tiles_used = partition$n,
          treated_idx_used = treated_idx,
          quiet = TRUE
        )
        summarize_boot(ate_c_boot, rep_id, sim_C$pre_df, sim_C$post_ctrl_df, sim_C$post_treat_df, c_params_boot)
      }, error = function(e) {
        list(ok = FALSE, rep = rep_id, msg = paste0("Bootstrap E failed: ", conditionMessage(e)))
      })
    }

    if ("F" %in% boot_targets_run) {
      out$F <- tryCatch({
        sim_D <- simulate_boot_data(d_ctrl_seed, d_treat_seed)
        d_params_boot <- F_params
        pp_post_d_boot <- sim_D$pp_post_bg_sim
        if (BOOT_REFIT_SCOPE %in% c("partial", "full")) {
          sem_boot <- run_sem_fit(
            pp_data_in = sim_D$pp_all_bg_sim,
            partition_in = partition,
            partition_processes_in = partition_processes,
            state_spaces_in = state_spaces,
            init_params_in = F_params,
            fixed_params_in = SENSITIVITY_FIXED_PARAMS,
            background_rate_var_in = "W",
            sem_inner_iter_in = BOOT_SEM_INNER_ITER,
            verbose_in = FALSE,
            label = sprintf("Boot F #%d", rep_id)
          )
          if (!is.null(sem_boot) && !is.null(sem_boot$etas_bivariate_params)) {
            d_params_boot <- sem_boot$etas_bivariate_params
          }
          if (!is.null(sem_boot) && !is.null(sem_boot$adaptive$adaptive_labelling)) {
            pp_post_d_boot <- sem_boot$adaptive$adaptive_labelling
            pp_post_d_boot <- pp_post_d_boot[pp_post_d_boot$t >= 0, , drop = FALSE]
          }
        }
        d_marg_boot <- extract_marginals(d_params_boot)
        ate_d_boot <- ate_estim_fast(
          d_marg_boot$ctrl, d_marg_boot$treat, pp_post_d_boot,
          label = sprintf("Boot F #%d", rep_id),
          filtration_history = sim_D$pre_df,
          crn_base_seed = boot_ate_crn_seed,
          phase = "bootstrap",
          n_tiles_used = partition$n,
          treated_idx_used = treated_idx,
          quiet = TRUE
        )
        summarize_boot(ate_d_boot, rep_id, sim_D$pre_df, sim_D$post_ctrl_df, sim_D$post_treat_df, d_params_boot)
      }, error = function(e) {
        list(ok = FALSE, rep = rep_id, msg = paste0("Bootstrap F failed: ", conditionMessage(e)))
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
    ok_targets <- names(out)[names(out) %in% c("E", "F")]
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
    model_res <- lapply(boot_results, function(z) z[[key]])
    model_res <- Filter(Negate(is.null), model_res)
    ok <- Filter(function(x) isTRUE(x$ok), model_res)
    fail <- Filter(function(x) !isTRUE(x$ok), model_res)
    rep_df <- if (length(ok) > 0) {
      do.call(rbind, lapply(ok, function(z) data.frame(
        rep = z$rep,
        n_pre_sim = z$n_pre_sim,
        n_post_ctrl_sim = z$n_post_ctrl_sim,
        n_post_treat_sim = z$n_post_treat_sim,
        ate_total_mean = z$ate_total_mean,
        ate_total_sd = z$ate_total_sd,
        ate_tile_mean = z$ate_tile_mean,
        eta_ctrl = z$eta_ctrl,
        eta_treat = z$eta_treat
      )))
    } else data.frame()
    list(
      replicate_summary = rep_df,
      params = lapply(ok, function(z) z$params),
      n_success = length(ok),
      n_fail = length(fail),
      failures = fail
    )
  }

  bootstrap_ate <- list(
    config = list(
      reps = BOOT_N_REPS,
      targets = BOOT_TARGETS,
      refit_scope = BOOT_REFIT_SCOPE,
      outer_cores = BOOT_OUTER_CORES,
      sem_inner_iter = BOOT_SEM_INNER_ITER,
      seed = BOOT_SEED
    )
  )
  if ("E" %in% boot_targets_run) bootstrap_ate$fit_E <- make_boot_block("E")
  if ("F" %in% boot_targets_run) bootstrap_ate$fit_F <- make_boot_block("F")
  # Backward-compatibility aliases expected by older report code.
  if (!is.null(bootstrap_ate$fit_E) && is.null(bootstrap_ate$fit_C)) bootstrap_ate$fit_C <- bootstrap_ate$fit_E
  if (!is.null(bootstrap_ate$fit_F) && is.null(bootstrap_ate$fit_D)) bootstrap_ate$fit_D <- bootstrap_ate$fit_F
rm(boot_results)
invisible(gc(verbose = FALSE))

  if (!is.null(bootstrap_ate$fit_E)) {
    bC <- bootstrap_ate$fit_E$replicate_summary
    if (nrow(bC) > 0) {
      cat(sprintf("  Bootstrap E complete: success=%d, fail=%d, mean total ATE=%.1f, SD(rep means)=%.1f\n",
                  nrow(bC), bootstrap_ate$fit_E$n_fail,
                  mean(bC$ate_total_mean, na.rm = TRUE),
                  stats::sd(bC$ate_total_mean, na.rm = TRUE)))
    } else {
      cat(sprintf("  Bootstrap E complete: success=%d, fail=%d\n",
                  bootstrap_ate$fit_E$n_success, bootstrap_ate$fit_E$n_fail))
    }
  }
  if (!is.null(bootstrap_ate$fit_F)) {
    bD <- bootstrap_ate$fit_F$replicate_summary
    if (nrow(bD) > 0) {
      cat(sprintf("  Bootstrap F complete: success=%d, fail=%d, mean total ATE=%.1f, SD(rep means)=%.1f\n",
                  nrow(bD), bootstrap_ate$fit_F$n_fail,
                  mean(bD$ate_total_mean, na.rm = TRUE),
                  stats::sd(bD$ate_total_mean, na.rm = TRUE)))
    } else {
      cat(sprintf("  Bootstrap F complete: success=%d, fail=%d\n",
                  bootstrap_ate$fit_F$n_success, bootstrap_ate$fit_F$n_fail))
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

# Sensitivity ATE payloads are already computed in Step 6a so they are
# included in pre-bootstrap checkpoints and reused for final results.

# ============================================================================
# 7. Summary tables
# ============================================================================
cat("\n")
cat("===========================================================================\n")
cat("                       PARAMETER COMPARISON TABLE\n")
cat("===========================================================================\n\n")

fmt <- function(x) sprintf("%8.4f", x)

biv_ctrl_val <- function(par, p) {
  if (p == "mu") par[["mu_0"]] else par[[paste0(p, "_00")]]
}
biv_treat_val <- function(par, p) {
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
  cat(sprintf("  %s: A_01=%.4f  A_10=%.4f  alpha_m_01=%.4f  alpha_m_10=%.4f\n",
    info$lab, info$par[["A_01"]], info$par[["A_10"]],
    info$par[["alpha_m_01"]], info$par[["alpha_m_10"]]))
}

cat("\n")
cat("===========================================================================\n")
cat(sprintf("     EARTHQUAKES SAVED COMPARISON (%d-day horizon)\n", ATE_WINDOW_DAYS))
cat("===========================================================================\n\n")
cat(sprintf("%-30s  %12s  %12s  %12s  %12s\n",
            "", "Total Saved", "SD(saved)", "eta_ctrl", "eta_treat"))
cat(paste(rep("-", 90), collapse = ""), "\n")

ate_print_list <- list(
  list(ate = ate_B, lab = "A: Naive Biv (county)"),
  list(ate = ate_D, lab = "B: SEM Biv (county)"),
  list(ate = ate_E, lab = "C: Naive Biv+KDE (county)"),
  list(ate = ate_F, lab = "D: SEM Biv+KDE (county)"))

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
  fitE = list(params = E_params, loglik = E_loglik, fit = fitE, ate = ate_E),
  fitF = list(params = F_params, ctrl = F_ctrl, treat = F_treat, sem = semF, ate = ate_F, louis = louis_F),
  fitI = list(params = I_params, fit = fitI, ate = ate_I),
  fitJ = list(params = J_params, fit = semJ, ctrl = J_ctrl, treat = J_treat, ate = ate_J),
  fits_named = fits_named,
  kde_variants = kde_variant_fits,
  fitG = NULL,
  fitH = NULL,
  bootstrap_ate = bootstrap_ate,
  partition_results = partition_results,
  ate_partitions = ate_partitions,
  kde_bandwidth_sensitivity = kde_bandwidth_sensitivity,
  kde_info = list(bw_sigma = as.numeric(bw_sigma), n_training = n_pre_holdout_ctrl,
                  n_pre_ctrl_holdout = n_pre_holdout_ctrl,
                  n_pre_holdout = nrow(pp_pre_holdout),
                  n_pre_used = nrow(pp_pre),
                  trigger_range_km = trigger_range_km),
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
    treated_names = treated_names,
    n_counties = partition$n,
    n_treated = sum(treated_idx)
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
    FIXED_STRUCTURAL = FIXED_STRUCTURAL,
    SEM_N_ITER = SEM_N_ITER, SEM_INNER_ITER = SEM_INNER_ITER,
    SEM_INNER_PROPS = SEM_INNER_PROPS,
    SEM_N_LABELLINGS = SEM_N_LABELLINGS,
    SEM_CHANGE_FACTOR = SEM_CHANGE_FACTOR,
    SEM_CHANGE_FACTOR_MIN_MULT = SEM_CHANGE_FACTOR_MIN_MULT,
    SEM_CHANGE_FACTOR_MAX_MULT = SEM_CHANGE_FACTOR_MAX_MULT,
    SEM_OPTIM_METHOD = SEM_OPTIM_METHOD,
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
    RUN_SENSITIVITY = RUN_SENSITIVITY,
    KDE_VARIANT_MODE = KDE_VARIANT_MODE,
    RUN_KDE_PROFILE_SWEEP = RUN_KDE_PROFILE_SWEEP,
    MEMORY_SAFE = MEMORY_SAFE,
    TRIM_SENS_OBJECTS = TRIM_SENS_OBJECTS,
    SENS_CORES = SENS_CORES,
    ATE_SIM_CORES = ATE_SIM_CORES,
    ATE_USE_CRN = OK_ATE_USE_CRN,
    ATE_CRN_PAIR = OK_ATE_CRN_PAIR,
    ATE_CONDITIONAL_ON_PRE = OK_ATE_CONDITIONAL_ON_PRE,
    ATE_CRN_BASE = OK_ATE_CRN_BASE,
    SEM_MAX_RELABEL_STEP_FRAC = SEM_MAX_RELABEL_STEP_FRAC,
    SEM_FORCE_PARAM_UPDATE_FLIP_FRAC = SEM_FORCE_PARAM_UPDATE_FLIP_FRAC,
    ATE_N_SIMS = ATE_N_SIMS, ATE_WINDOW_DAYS = ATE_WINDOW_DAYS,
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
    BOOT_BRANCHING_MAX = BOOT_BRANCHING_MAX,
    BOOT_MAX_PRE_EVENTS = BOOT_MAX_PRE_EVENTS,
    BOOT_MAX_POST_EVENTS_PER_PROC = BOOT_MAX_POST_EVENTS_PER_PROC,
    BOOT_MAX_TOTAL_EVENTS = BOOT_MAX_TOTAL_EVENTS,
    TEST_MODE = TEST_MODE,
    windowT_post = windowT_post,
    n_pre = nrow(pp_pre), n_pre_holdout = nrow(pp_pre_holdout), n_pre_total = nrow(pp_pre_all),
    n_post = nrow(pp_post),
    n_counties = partition$n, n_treated = sum(treated_idx),
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
          c("render", basename(report_file), "--to", fmt),
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
        file.path(SCRIPT_DIR, "oklahoma_report.html")
      } else {
        NA_character_
      }
      report_pdf_rendered <- if ("pdf" %in% REPORT_FORMATS) {
        file.path(SCRIPT_DIR, "oklahoma_report.pdf")
      } else {
        NA_character_
      }
      report_html_target <- if (is.character(report_html_rendered) &&
                                !is.na(report_html_rendered)) {
        file.path(SCRIPT_DIR, add_file_tag(sprintf("oklahoma_report_%s.html", result_stamp)))
      } else {
        report_html_rendered
      }
      report_pdf_target <- if (is.character(report_pdf_rendered) &&
                               !is.na(report_pdf_rendered)) {
        file.path(SCRIPT_DIR, add_file_tag(sprintf("oklahoma_report_%s.pdf", result_stamp)))
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
    }
  } else {
    cat("Quarto not found in PATH; skipping report render.\n")
  }
} else if (file.exists(report_file)) {
  cat("Report render skipped (OK_REPORT_FORMATS empty).\n")
}

cat("\n=== Oklahoma Analysis Complete ===\n")
