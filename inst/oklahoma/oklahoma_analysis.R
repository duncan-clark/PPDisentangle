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
# Canonical output path at repo root, with mirror to legacy script-local output.
OUT_DIR  <- file.path(REPO_DIR, "output")
LEGACY_OUT_DIR <- file.path(REPO_DIR, "inst", "oklahoma", "output")
PLOT_DIR <- file.path(OUT_DIR, "plots")
for (d in c(OUT_DIR, PLOT_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
if (!dir.exists(LEGACY_OUT_DIR)) dir.create(LEGACY_OUT_DIR, recursive = TRUE)
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

SEM_N_LABELLINGS  <- if (QUICK_CHECK) 3 else if (TEST_MODE) 5  else 30
SEM_N_ITER        <- 1
SEM_INNER_ITER    <- if (QUICK_CHECK) 2 else if (TEST_MODE) 5 else 1000
env_sem_inner_iter <- suppressWarnings(as.integer(Sys.getenv("OK_SEM_INNER_ITER", "")))
if (!is.na(env_sem_inner_iter) && env_sem_inner_iter > 0L) {
  SEM_INNER_ITER <- env_sem_inner_iter
}
SENS_SEM_INNER_ITER <- suppressWarnings(as.integer(Sys.getenv("OK_SENS_SEM_INNER_ITER", as.character(SEM_INNER_ITER))))
if (!is.finite(SENS_SEM_INNER_ITER) || is.na(SENS_SEM_INNER_ITER) || SENS_SEM_INNER_ITER < 1L) {
  SENS_SEM_INNER_ITER <- SEM_INNER_ITER
}
SEM_INNER_PROPS   <- if (QUICK_CHECK) 3 else if (TEST_MODE) 5  else 20
SEM_CHANGE_FACTOR <- 0.01
SEM_TEMPORAL_WEIGHT <- if (QUICK_CHECK) 0.2 else if (TEST_MODE) 0.4 else 0.6
SEM_TEMPORAL_SCALE_DAYS <- 15
SEM_PARAM_UPDATE  <- if (QUICK_CHECK) 10 else if (TEST_MODE) 10 else 25
SEM_OUTER_MAXIT       <- if (QUICK_CHECK) 40 else if (TEST_MODE) 200 else 220
SEM_OUTER_MAXIT_BIV   <- if (QUICK_CHECK) 30 else if (TEST_MODE) 100 else 60
# Keep Decode iterations separate from SEM inner iterations.
DECODE_ITER <- if (QUICK_CHECK) 2 else if (TEST_MODE) 5 else 200
env_decode_iter <- suppressWarnings(as.integer(Sys.getenv("OK_DECODE_ITER", "")))
if (!is.na(env_decode_iter) && env_decode_iter > 0L) {
  DECODE_ITER <- env_decode_iter
}
# Leave Decode capability intact, but default off for long runs unless explicitly enabled.
RUN_DECODE <- tolower(Sys.getenv("OK_RUN_DECODE", "false")) %in% c("1", "true", "yes", "y")
# Optional speed mode for proof-of-concept runs.
RUN_SENSITIVITY <- tolower(Sys.getenv("OK_RUN_SENSITIVITY", "true")) %in% c("1", "true", "yes", "y")
OK_VERBOSE <- tolower(Sys.getenv("OK_VERBOSE", "false")) %in% c("1", "true", "yes", "y")

STRUCT_DEFAULTS <- list(c = 0.05, p = 1.2, D = 5.0, gamma = 0.5, q = 1.5)
# Oklahoma fits now estimate all ETAS parameters freely.
FIXED_STRUCTURAL <- list()

ATE_N_SIMS    <- if (QUICK_CHECK) 5 else if (TEST_MODE) 20 else 40
ATE_WINDOW_DAYS <- 100

RUN_BOOTSTRAP_ATE <- tolower(Sys.getenv("OK_RUN_BOOTSTRAP_ATE", "false")) %in% c("1", "true", "yes", "y")
BOOT_N_REPS <- suppressWarnings(as.integer(Sys.getenv("OK_BOOT_N_REPS", "0")))
if (!is.finite(BOOT_N_REPS) || is.na(BOOT_N_REPS) || BOOT_N_REPS < 0L) BOOT_N_REPS <- 0L
BOOT_TARGETS_RAW <- toupper(Sys.getenv("OK_BOOT_TARGETS", "E,F"))
BOOT_TARGETS <- unique(trimws(unlist(strsplit(BOOT_TARGETS_RAW, ","))))
BOOT_TARGETS <- BOOT_TARGETS[BOOT_TARGETS %in% c("E", "F")]
if (length(BOOT_TARGETS) < 1) BOOT_TARGETS <- c("E", "F")
BOOT_REFIT_SCOPE <- tolower(trimws(Sys.getenv("OK_BOOT_REFIT_SCOPE", "partial")))
if (!BOOT_REFIT_SCOPE %in% c("none", "partial", "full")) BOOT_REFIT_SCOPE <- "partial"
BOOT_SEM_INNER_ITER <- suppressWarnings(as.integer(Sys.getenv("OK_BOOT_SEM_INNER_ITER", "200")))
if (!is.finite(BOOT_SEM_INNER_ITER) || is.na(BOOT_SEM_INNER_ITER) || BOOT_SEM_INNER_ITER < 1L) {
  BOOT_SEM_INNER_ITER <- 200L
}
BOOT_OUTER_CORES_RAW <- Sys.getenv("OK_BOOT_OUTER_CORES", "")
BOOT_SEED <- suppressWarnings(as.integer(Sys.getenv("OK_BOOT_SEED", "")))
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
  N_CORES <- min(N_CORES, 4L)
}
BOOT_OUTER_DEFAULT <- if (MEMORY_SAFE) 1L else N_CORES
BOOT_OUTER_CORES <- suppressWarnings(as.integer(ifelse(nzchar(BOOT_OUTER_CORES_RAW), BOOT_OUTER_CORES_RAW, as.character(BOOT_OUTER_DEFAULT))))
if (!is.finite(BOOT_OUTER_CORES) || is.na(BOOT_OUTER_CORES) || BOOT_OUTER_CORES < 1L) BOOT_OUTER_CORES <- 1L
BOOT_OUTER_CORES <- max(1L, min(BOOT_OUTER_CORES, N_CORES))
if (MEMORY_SAFE && RUN_BOOTSTRAP_ATE && BOOT_N_REPS > 0L && BOOT_OUTER_CORES > 1L) {
  cat(sprintf("Memory-safe bootstrap override: forcing BOOT_OUTER_CORES from %d to 1\n", BOOT_OUTER_CORES))
  BOOT_OUTER_CORES <- 1L
}
SENS_CORES_DEFAULT <- if (MEMORY_SAFE) min(2L, N_CORES) else N_CORES
SENS_CORES <- suppressWarnings(as.integer(ifelse(nzchar(SENS_CORES_RAW), SENS_CORES_RAW, as.character(SENS_CORES_DEFAULT))))
if (!is.finite(SENS_CORES) || is.na(SENS_CORES) || SENS_CORES < 1L) SENS_CORES <- 1L
SENS_CORES <- max(1L, min(SENS_CORES, N_CORES))
ATE_SIM_CORES_DEFAULT <- if (MEMORY_SAFE) min(2L, N_CORES) else N_CORES
ATE_SIM_CORES <- suppressWarnings(as.integer(ifelse(nzchar(ATE_SIM_CORES_RAW), ATE_SIM_CORES_RAW, as.character(ATE_SIM_CORES_DEFAULT))))
if (!is.finite(ATE_SIM_CORES) || is.na(ATE_SIM_CORES) || ATE_SIM_CORES < 1L) ATE_SIM_CORES <- 1L
ATE_SIM_CORES <- max(1L, min(ATE_SIM_CORES, N_CORES))
PARALLEL_BACKEND <- tolower(trimws(Sys.getenv(
  "OK_PARALLEL_BACKEND",
  if (MEMORY_SAFE) "psock" else "fork"
)))
if (!PARALLEL_BACKEND %in% c("fork", "psock", "sequential")) PARALLEL_BACKEND <- "psock"

run_parallel <- function(X, FUN, cores, label = "job") {
  n <- length(X)
  cores_use <- max(1L, min(as.integer(cores), as.integer(n)))
  if (n <= 1L || cores_use <= 1L || identical(PARALLEL_BACKEND, "sequential")) {
    return(lapply(X, FUN))
  }
  if (identical(PARALLEL_BACKEND, "fork")) {
    return(parallel::mclapply(X, FUN, mc.cores = cores_use))
  }
  cl <- parallel::makePSOCKcluster(cores_use)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
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
  tryCatch(
    parallel::parLapply(cl, X, FUN),
    error = function(e) {
      warning(sprintf("PSOCK failed for %s; falling back to sequential: %s", label, e$message))
      lapply(X, FUN)
    }
  )
}

mode_label <- if (QUICK_CHECK) "QUICK_CHECK" else if (TEST_MODE) "TEST" else "FULL"
cat("=== Oklahoma County-Based ETAS Analysis ===\n")
cat(sprintf("Mode: %s | SEM iters: %d | Change factor: %.3f | Cores: %d\n",
            mode_label, SEM_N_ITER, SEM_CHANGE_FACTOR,
            N_CORES))
cat(sprintf("Memory safe: %s | Sens cores: %d | ATE sim cores: %d | Boot outer cores: %d | Trim sensitivity objects: %s\n",
            MEMORY_SAFE, SENS_CORES, ATE_SIM_CORES, BOOT_OUTER_CORES, TRIM_SENS_OBJECTS))
cat(sprintf("Parallel backend: %s\n", PARALLEL_BACKEND))
cat(sprintf("SEM inner iters: main=%d, sensitivity=%d, bootstrap=%d\n",
            SEM_INNER_ITER, SENS_SEM_INNER_ITER, BOOT_SEM_INNER_ITER))
cat(sprintf("Verbose optimizer/SEM tracing: %s\n", OK_VERBOSE))

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

# Data-driven structural parameter calibration from first-half control pre data.
# These values are used as starting points (not fixed) for downstream fits.
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
FIXED_STRUCTURAL <- list()
if (is.null(PRE_CTRL_BOOT_PARAMS)) {
  PRE_CTRL_BOOT_PARAMS <- list(mu = 1.0, A = 0.2, alpha_m = 0.8,
                               c = STRUCT_INIT$c, p = STRUCT_INIT$p,
                               D = STRUCT_INIT$D, gamma = STRUCT_INIT$gamma, q = STRUCT_INIT$q)
}

cat(sprintf("  Structural init (from first 50%% control pre): c=%.4f, p=%.4f, D=%.4f, gamma=%.4f, q=%.4f\n",
            STRUCT_INIT$c, STRUCT_INIT$p, STRUCT_INIT$D,
            STRUCT_INIT$gamma, STRUCT_INIT$q))

apply_structural_init_etas <- function(start_par) {
  out <- start_par
  for (nm in c("c", "p", "D", "gamma", "q")) out[[nm]] <- STRUCT_INIT[[nm]]
  out
}

apply_structural_init_biv <- function(par_vec) {
  out <- par_vec
  out[c("c", "p", "D", "gamma", "q")] <- unlist(STRUCT_INIT[c("c", "p", "D", "gamma", "q")])
  out
}

cat(sprintf("  Events in OK: pre=%d, post=%d, total=%d\n",
            nrow(pp_pre), nrow(pp_post), nrow(pp_all)))
cat(sprintf("  Pre-treatment split: holdout=%d (first 50%%), estimation=%d (second 50%%)\n",
            nrow(pp_pre_holdout), nrow(pp_pre)))
cat(sprintf("  Post-treatment: %d control, %d treated\n",
            sum(pp_post$location_process == "control"),
            sum(pp_post$location_process == "treated")))

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
  partition_map_file <- add_file_tag("partition_map.png")
  ggsave(file.path(PLOT_DIR, partition_map_file), p_partition,
         width = 12, height = 7, dpi = 150)
  cat(sprintf("  Saved %s\n", partition_map_file))
}, error = function(e) cat("  Partition plot error:", e$message, "\n"))

tryCatch({
  p_pre <- ggplot(pp_pre, aes(x = x, y = y, alpha = t)) +
    geom_point(color = "#377eb8", size = 1.5, shape = 21, fill = "#377eb8") +
    scale_alpha_continuous(name = "Time (days)", range = c(0.15, 0.8)) +
    labs(title = "Pre-treatment Events", x = "X (km)", y = "Y (km)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  pre_treatment_file <- add_file_tag("pp_pre_treatment.png")
  ggsave(file.path(PLOT_DIR, pre_treatment_file), p_pre,
         width = 12, height = 7, dpi = 150)
  cat(sprintf("  Saved %s\n", pre_treatment_file))

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
  post_treatment_file <- add_file_tag("pp_post_treatment.png")
  ggsave(file.path(PLOT_DIR, post_treatment_file), p_post,
         width = 12, height = 7, dpi = 150)
  cat(sprintf("  Saved %s\n", post_treatment_file))
}, error = function(e) cat("  PP plot error:", e$message, "\n"))

# ============================================================================
# 4A. Fit A: Naive independent ETAS  [DISABLED — kept for reference]
# ============================================================================
# Univariate fits are no longer run. The bivariate model is the focus.

fit_best_indep <- function(realiz, zbr, starts, maxit) {
  best_fit <- NULL; best_val <- -Inf
  starts_use <- lapply(starts, apply_structural_init_etas)
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
A_ctrl <- apply_structural_init_etas(A_ctrl)
A_treat <- apply_structural_init_etas(A_treat)
names(A_ctrl) <- etas_names; names(A_treat) <- etas_names

cat("  Control:", paste(etas_names, round(unlist(A_ctrl), 4), sep = "=", collapse = ", "), "\n")
cat("  Treated:", paste(etas_names, round(unlist(A_treat), 4), sep = "=", collapse = ", "), "\n")

# ============================================================================
# 4B. Fit B: Naive bivariate ETAS
# ============================================================================
cat("\n--- Fit B: Naive bivariate ETAS ---\n")

biv_init <- apply_structural_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
biv_names <- names(biv_init)
fit_b <- function() {
  tryCatch({
    fit_etas_bivariate(
      params_init = biv_init, realiz = pp_all,
      windowT = windowT_fit, windowS = win_km, m0 = ETAS_M0,
      control_state_space = control_ss, treated_state_space = treated_ss,
      treated_background_zero_before = 0,
      maxit = VANILLA_MAXIT, fixed_params = NULL, trace = if (OK_VERBOSE) 1 else 0
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
# 4D. Fit D: SEM bivariate ETAS
# ============================================================================
cat("\n--- Fit D: SEM bivariate ETAS ---\n")

biv_init_D <- apply_structural_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
biv_fixed <- NULL

run_sem_fit <- function(pp_data_in,
                        partition_in,
                        partition_processes_in,
                        state_spaces_in,
                        init_params_in,
                        background_rate_var_in = NULL,
                        use_pre_history_for_biv_in = TRUE,
                        treated_background_zero_before_in = 0,
                        sem_inner_iter_in = SEM_INNER_ITER,
                        verbose_in = FALSE,
                        label = "SEM") {
  t0 <- proc.time()[["elapsed"]]
  cat(sprintf("  [%s] start (n=%d, pid=%d)\n", label, nrow(pp_data_in), Sys.getpid()))
  out <- tryCatch({
    adaptive_SEM(
      pp_data = pp_data_in, partition = partition_in,
      partition_processes = partition_processes_in,
      statespace = win_km, time_window = windowT_post, treatment_time = 0,
      hawkes_params_control = A_ctrl, hawkes_params_treated = A_treat,
      N_labellings = SEM_N_LABELLINGS, N_iter = SEM_N_ITER, verbose = verbose_in,
      model_type = "etas_bivariate",
      adaptive_control = list(
        param_update_cadence = SEM_PARAM_UPDATE,
        proposal_update_cadence = 1,
        state_spaces = state_spaces_in,
        iter = sem_inner_iter_in, n_props = SEM_INNER_PROPS,
        change_factor = SEM_CHANGE_FACTOR, verbose = FALSE,
        temporal_weight = SEM_TEMPORAL_WEIGHT,
        temporal_scale_days = SEM_TEMPORAL_SCALE_DAYS,
        update_starting_data = TRUE, include_starting_data = TRUE,
        update_control_params = TRUE, fixed_params = biv_fixed,
        proposal_method = "simulation",
        outer_maxit = SEM_OUTER_MAXIT, outer_maxit_biv = SEM_OUTER_MAXIT_BIV
      ),
      m0 = ETAS_M0, beta_gr = BETA_GR,
      etas_bivariate_params = init_params_in,
      background_rate_var = background_rate_var_in,
      use_pre_history_for_biv = use_pre_history_for_biv_in,
      treated_background_zero_before = treated_background_zero_before_in
    )
  }, error = function(e) {
    cat(sprintf("  [%s] error: %s\n", label, e$message))
    NULL
  })
  t1 <- proc.time()[["elapsed"]]
  cat(sprintf("  [%s] done in %.1fs\n", label, t1 - t0))
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
      verbose_in = OK_VERBOSE,
      label = "Fit D"
    )
  }, error = function(e) { cat("  SEM-biv error:", e$message, "\n"); NULL })
}

t_fit_B <- proc.time()[["elapsed"]]
fitB <- fit_b()
fit_B_elapsed <- proc.time()[["elapsed"]] - t_fit_B
add_timing_row(
  stage = "fit_B_naive_bivariate",
  elapsed_sec = fit_B_elapsed,
  status = if (!is.null(fitB)) "ok" else "failed"
)
semD <- NULL

B_params <- if (!is.null(fitB)) fitB$par else biv_init
B_loglik <- if (!is.null(fitB)) fitB$value else NA
cat("  Params:", paste(biv_names, round(B_params, 4), sep = "=", collapse = ", "), "\n")
D_params <- biv_init_D
D_ctrl <- A_ctrl
D_treat <- A_treat

# ============================================================================
# 4E. Non-parametric background rate from first 50% of control pre-treatment
# ============================================================================
cat("\n--- Step 4E: KDE background rate from held-out control pre-treatment data ---\n")

pp_pre_holdout_ctrl <- pp_pre_holdout[pp_pre_holdout$location_process == "control", ]
pp_pre_holdout_ctrl <- pp_pre_holdout_ctrl[order(pp_pre_holdout_ctrl$t), ]
n_pre_holdout_ctrl <- nrow(pp_pre_holdout_ctrl)
if (n_pre_holdout_ctrl < 2) {
  stop("Insufficient held-out control pre-treatment events for KDE background estimation.")
}

cat(sprintf("  Held-out control pre-treatment events for KDE: %d\n",
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
         subtitle = "Estimated from held-out first 50% of pre-treatment (control only)",
         x = "X (km)", y = "Y (km)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  background_rate_file <- add_file_tag("background_rate_kde.png")
  ggsave(file.path(PLOT_DIR, background_rate_file), p_bg,
         width = 10, height = 6, dpi = 150)
  cat(sprintf("  Saved %s\n", background_rate_file))
}, error = function(e) cat("  Background-rate plot error:", e$message, "\n"))

# ============================================================================
# 4F. Fit E: Naive bivariate ETAS with KDE background
# ============================================================================
cat("\n--- Fit E: Naive bivariate ETAS with KDE background ---\n")

biv_init_E <- apply_structural_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
fit_e <- function() {
  tryCatch({
    fit_etas_bivariate(
      params_init = biv_init_E, realiz = pp_all_bg,
      windowT = windowT_fit, windowS = win_km, m0 = ETAS_M0,
      control_state_space = control_ss, treated_state_space = treated_ss,
      background_rate_var = "W",
      treated_background_zero_before = 0,
      maxit = VANILLA_MAXIT, fixed_params = NULL, trace = if (OK_VERBOSE) 1 else 0
    )
  }, error = function(e) { cat("  Bivariate+KDE fit error:", e$message, "\n"); NULL })
}

# ============================================================================
# 4G. Fit F: SEM bivariate ETAS with KDE background
# ============================================================================
cat("\n--- Fit F: SEM bivariate ETAS with KDE background ---\n")

biv_init_F <- apply_structural_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
fit_f <- function() {
  tryCatch({
    run_sem_fit(
      pp_data_in = pp_all_bg,
      partition_in = partition,
      partition_processes_in = partition_processes,
      state_spaces_in = state_spaces,
      init_params_in = biv_init_F,
      background_rate_var_in = "W",
      verbose_in = OK_VERBOSE,
      label = "Fit F"
    )
  }, error = function(e) { cat("  SEM-biv+KDE error:", e$message, "\n"); NULL })
}

t_fit_E <- proc.time()[["elapsed"]]
fitE <- fit_e()
fit_E_elapsed <- proc.time()[["elapsed"]] - t_fit_E
add_timing_row(
  stage = "fit_E_naive_bivariate_kde",
  elapsed_sec = fit_E_elapsed,
  status = if (!is.null(fitE)) "ok" else "failed"
)
semF <- NULL

E_params <- if (!is.null(fitE)) fitE$par else biv_init_E
E_loglik <- if (!is.null(fitE)) fitE$value else NA
cat("  Params:", paste(biv_names, round(E_params, 4), sep = "=", collapse = ", "), "\n")
F_params <- biv_init_F
F_ctrl <- A_ctrl
F_treat <- A_treat

# ============================================================================
# 4H. Decode: hard-EM / coordinate-ascent relabelling (county only)
# ============================================================================
cat("\n--- Fit G/H: Decode hard-EM relabelling (county only) ---\n")
if (!RUN_DECODE) {
  cat("  Decode disabled (set OK_RUN_DECODE=true to enable).\n")
}

decode_hard_em_bivariate <- function(start_labelling,
                                     init_params,
                                     label = "Decode",
                                     background_rate_var = NULL,
                                     treated_background_zero_before = 0,
                                     decode_iter = SEM_INNER_ITER,
                                     refit_cadence = 10L,
                                     max_target_cells = 12L,
                                     max_flips_per_cell = 8L,
                                     icm_passes = 3L,
                                     icm_cadence = 25L,
                                     target_update_cadence = 5L,
                                     implied_sims = 2L) {
  if (is.null(start_labelling) || nrow(start_labelling) == 0 || is.null(init_params)) {
    cat(sprintf("  [%s] skipped (missing start labels/params)\n", label))
    return(NULL)
  }
  dat <- start_labelling[order(start_labelling$t), , drop = FALSE]
  if (is.null(dat$inferred_process)) dat$inferred_process <- dat$location_process
  pre <- dat[dat$t < 0, , drop = FALSE]
  post <- dat[dat$t >= 0, , drop = FALSE]
  if (nrow(post) == 0) {
    cat(sprintf("  [%s] skipped (no post-treatment points)\n", label))
    return(NULL)
  }
  pre$location_process <- "control"
  pre$inferred_process <- "control"
  history_window_start <- if (nrow(pre) > 0) min(pre$t, na.rm = TRUE) else windowT_post[1]
  history_window <- c(history_window_start, windowT_post[2])
  n_pre <- nrow(pre)
  n_post <- nrow(post)

  fix_names <- names(FIXED_STRUCTURAL)
  free_idx <- which(!names(init_params) %in% fix_names)
  fixed_idx <- which(names(init_params) %in% fix_names)

  rr_full_template <- rbind(pre, post)
  rr_full_template$..post_idx <- c(rep(NA_integer_, n_pre), seq_len(n_post))
  rr_full_template <- rr_full_template[order(rr_full_template$t), , drop = FALSE]
  post_rows <- which(!is.na(rr_full_template$..post_idx))
  post_label_idx <- as.integer(rr_full_template$..post_idx[post_rows])

  eval_obj <- function(labels, params) {
    rr_full <- rr_full_template
    rr_full$inferred_process[post_rows] <- labels[post_label_idx]
    rr_full$..post_idx <- NULL
    loglik_etas_bivariate(
      params = params, realiz = rr_full,
      windowT = history_window, windowS = win_km,
      control_state_space = control_ss, treated_state_space = treated_ss,
      background_rate_var = background_rate_var,
      treated_background_zero_before = treated_background_zero_before
    )
  }

  obj_post <- function(labels, params) {
    eval_obj(labels, params)
  }

  refit_params <- function(labels, params_start) {
    obj <- function(p15) {
      eval_obj(labels, p15)
    }
    p_cur <- params_start
    if (length(fixed_idx) > 0) p_cur[fixed_idx] <- params_start[fixed_idx]
    if (length(free_idx) == 0) return(p_cur)
    wrap <- function(pf) {
      p <- p_cur
      p[free_idx] <- pf
      obj(p)
    }
    fit <- tryCatch(
      optim(par = p_cur[free_idx], fn = wrap, method = "Nelder-Mead",
            control = list(fnscale = -1, trace = 0, maxit = 250)),
      error = function(e) NULL
    )
    if (is.null(fit)) return(p_cur)
    p_out <- p_cur
    p_out[free_idx] <- fit$par
    p_out
  }

  implied_treated_counts <- function(params, n_sims = implied_sims) {
    nb <- partition$n
    agg <- numeric(nb)
    for (s in seq_len(max(1L, n_sims))) {
      sim <- tryCatch(
        generate_inhomogeneous_etas_bivariate(
          Omega = win_km, partition = partition, time_window = windowT_post,
          partition_processes = partition_processes,
          etas_bivariate_params = params, m0 = ETAS_M0, beta_gr = BETA_GR,
          state_spaces = state_spaces, filtration = pre
        ),
        error = function(e) NULL
      )
      if (is.null(sim) || nrow(sim) == 0) next
      sim <- sim[sim$t >= 0, , drop = FALSE]
      if (nrow(sim) == 0) next
      ti <- as.integer(tileindex(sim$x, sim$y, partition))
      treated_rows <- if ("process" %in% names(sim)) sim$process == "treated" else sim$process_id == 1
      agg <- agg + tabulate(ti[treated_rows], nbins = nb)
    }
    agg / max(1L, n_sims)
  }

  labels_cur <- as.character(post$inferred_process)
  params_cur <- refit_params(labels_cur, init_params)
  obj_cur <- obj_post(labels_cur, params_cur)
  if (!is.finite(obj_cur)) obj_cur <- obj_post(labels_cur, init_params)
  if (!is.finite(obj_cur)) {
    cat(sprintf("  [%s] skipped (non-finite objective)\n", label))
    return(NULL)
  }

  tile_idx <- as.integer(tileindex(post$x, post$y, partition))
  obs_treat <- tabulate(tile_idx[post$location_process == "treated"], nbins = partition$n)
  n_flips <- 0L
  trace_avg_flips <- numeric(0)
  trace_max_flips <- numeric(0)
  trace_objective <- numeric(0)
  trace_phase <- character(0)
  add_trace <- function(avg_flips, max_flips, objective, phase) {
    trace_avg_flips <<- c(trace_avg_flips, avg_flips)
    trace_max_flips <<- c(trace_max_flips, max_flips)
    trace_objective <<- c(trace_objective, objective)
    trace_phase <<- c(trace_phase, phase)
  }
  eval_one_flip <- function(j, labels, params, obj0) {
    prop <- labels
    prop[j] <- if (prop[j] == "control") "treated" else "control"
    obj1 <- obj_post(prop, params)
    c(delta = obj1 - obj0, obj = obj1)
  }

  target_tiles <- integer(0)
  refit_cad <- max(1L, as.integer(refit_cadence))
  for (iter_idx in seq_len(max(1L, decode_iter))) {
    if (iter_idx == 1L || (iter_idx %% refit_cad) == 0L) {
      params_new <- refit_params(labels_cur, params_cur)
      obj_new <- obj_post(labels_cur, params_new)
      if (is.finite(obj_new) && obj_new >= obj_cur) {
        params_cur <- params_new
        obj_cur <- obj_new
      }
    }

    if (iter_idx == 1L || (iter_idx %% max(1L, as.integer(target_update_cadence))) == 0L) {
      imp_treat <- implied_treated_counts(params_cur)
      discrepancy <- abs(obs_treat - imp_treat)
      target_tiles <- which(discrepancy > 0)
      if (length(target_tiles) > 0) {
        ord <- order(discrepancy[target_tiles], decreasing = TRUE)
        target_tiles <- target_tiles[ord][seq_len(min(length(target_tiles), max_target_cells))]
      }
    }

    flips_iter <- 0L

    for (tile in target_tiles) {
      idx <- which(tile_idx == tile)
      if (length(idx) == 0) next
      for (k in seq_len(max_flips_per_cell)) {
        improved <- FALSE
        idx_scan <- sample(idx, length(idx), replace = FALSE)
        for (j in idx_scan) {
          res <- eval_one_flip(j, labels_cur, params_cur, obj_cur)
          if (is.finite(res["delta"]) && res["delta"] > 0) {
            labels_cur[j] <- if (labels_cur[j] == "control") "treated" else "control"
            obj_cur <- res["obj"]
            flips_iter <- flips_iter + 1L
            n_flips <- n_flips + 1L
            improved <- TRUE
            break
          }
        }
        if (!improved) break
      }
    }

    do_icm <- (iter_idx %% max(1L, as.integer(icm_cadence))) == 0L || iter_idx == max(1L, decode_iter)
    if (do_icm) {
      for (pass in seq_len(max(1L, icm_passes))) {
        improved <- FALSE
        for (j in seq_len(length(labels_cur))) {
          res <- eval_one_flip(j, labels_cur, params_cur, obj_cur)
          if (is.finite(res["delta"]) && res["delta"] > 0) {
            labels_cur[j] <- if (labels_cur[j] == "control") "treated" else "control"
            obj_cur <- res["obj"]
            improved <- TRUE
            flips_iter <- flips_iter + 1L
            n_flips <- n_flips + 1L
          }
        }
        if (!improved) break
      }
    }

    add_trace(flips_iter, flips_iter, obj_cur, sprintf("decode_iter_%d", iter_idx))
    if (iter_idx == 1L || iter_idx %% 10L == 0L || iter_idx == max(1L, decode_iter)) {
      cat(sprintf("  [%s] iter %d/%d | flips_this_iter=%d | objective=%.3f\n",
                  label, iter_idx, max(1L, decode_iter), flips_iter, obj_cur))
    }
  }

  post_dec <- post
  post_dec$inferred_process <- labels_cur
  final_lab <- rbind(pre, post_dec)
  final_lab <- final_lab[order(final_lab$t), , drop = FALSE]

  cat(sprintf("  [%s] complete: %d flips over %d iterations, objective=%.3f\n",
              label, n_flips, max(1L, decode_iter), obj_cur))
  list(
    labelling = final_lab, params = params_cur, n_flips = n_flips, objective = obj_cur,
    trace = list(
      average_flips = trace_avg_flips,
      max_metric_flips = trace_max_flips,
      objective = trace_objective,
      phase = trace_phase
    )
  )
}

decode_D <- NULL
decode_F <- NULL
G_params <- NULL
H_params <- NULL
pp_post_decode_G <- NULL
pp_post_decode_H <- NULL
cat("\n--- Step 4H/4J: SEM (D/F) + Decode (G/H) joint fit dispatch ---\n")
fit_jobs <- c("D", "F")
if (RUN_DECODE) fit_jobs <- c(fit_jobs, "G", "H")

run_one_fit_job <- function(tag) {
  t0 <- proc.time()[["elapsed"]]
  out_obj <- NULL
  if (tag == "D") {
    out_obj <- fit_d()
  }
  if (tag == "F" && is.null(out_obj)) {
    out_obj <- fit_f()
  }
  if (tag == "G" && is.null(out_obj)) {
    out_obj <- decode_hard_em_bivariate(
      start_labelling = pp_all,
      init_params = B_params,
      label = "Decode-SEM Biv",
      background_rate_var = NULL,
      decode_iter = DECODE_ITER,
      icm_passes = 1L,
      icm_cadence = 25L,
      target_update_cadence = 5L,
      implied_sims = 1L
    )
  }
  if (tag == "H" && is.null(out_obj)) {
    out_obj <- decode_hard_em_bivariate(
      start_labelling = pp_all_bg,
      init_params = E_params,
      label = "Decode-SEM Biv+KDE",
      background_rate_var = "W",
      decode_iter = DECODE_ITER,
      icm_passes = 1L,
      icm_cadence = 25L,
      target_update_cadence = 5L,
      implied_sims = 1L
    )
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  list(tag = tag, obj = out_obj, elapsed = elapsed)
}

if (!TEST_MODE && !QUICK_CHECK && N_CORES > 1 && length(fit_jobs) > 1) {
  # Single outer parallel layer only.
  fit_out <- run_parallel(
    fit_jobs, run_one_fit_job,
    cores = min(length(fit_jobs), N_CORES),
    label = "fit-jobs"
  )
} else {
  fit_out <- lapply(fit_jobs, run_one_fit_job)
}

get_fit_job <- function(tag) {
  idx <- which(vapply(fit_out, function(z) identical(z$tag, tag), logical(1)))
  if (length(idx) == 0) return(NULL)
  fit_out[[idx[1]]]$obj
}

get_fit_elapsed <- function(tag) {
  idx <- which(vapply(fit_out, function(z) identical(z$tag, tag), logical(1)))
  if (length(idx) == 0) return(NA_real_)
  fit_out[[idx[1]]]$elapsed
}

semD <- get_fit_job("D")
semF <- get_fit_job("F")
if (RUN_DECODE) {
  decode_D <- get_fit_job("G")
  decode_F <- get_fit_job("H")
}

if (!is.null(semD)) {
  D_params <- semD$etas_bivariate_params
  D_ctrl <- semD$hawkes_params_control
  D_treat <- semD$hawkes_params_treated
  cat("  SEM-biv completed in", round(semD$time, 1), "s\n")
  cat("  Bivariate params:", paste(biv_names, round(D_params, 4), sep = "=", collapse = ", "), "\n")
} else {
  D_params <- biv_init_D
  D_ctrl <- A_ctrl
  D_treat <- A_treat
  cat("  SEM-biv failed, falling back to naive.\n")
}
fit_D_elapsed <- get_fit_elapsed("D")
add_timing_row(
  stage = "fit_D_sem_bivariate",
  elapsed_sec = fit_D_elapsed,
  status = if (!is.null(semD)) "ok" else "failed",
  detail = "elapsed from joint fit dispatch"
)

if (!is.null(semF)) {
  F_params <- semF$etas_bivariate_params
  F_ctrl <- semF$hawkes_params_control
  F_treat <- semF$hawkes_params_treated
  cat("  SEM-biv+KDE completed in", round(semF$time, 1), "s\n")
  cat("  Bivariate params:", paste(biv_names, round(F_params, 4), sep = "=", collapse = ", "), "\n")
} else {
  F_params <- biv_init_F
  F_ctrl <- A_ctrl
  F_treat <- A_treat
  cat("  SEM-biv+KDE failed, falling back to naive.\n")
}
fit_F_elapsed <- get_fit_elapsed("F")
add_timing_row(
  stage = "fit_F_sem_bivariate_kde",
  elapsed_sec = fit_F_elapsed,
  status = if (!is.null(semF)) "ok" else "failed",
  detail = "elapsed from joint fit dispatch"
)

if (RUN_DECODE) {
  decode_D_time <- get_fit_elapsed("G")
  decode_F_time <- get_fit_elapsed("H")
  if (is.finite(decode_D_time)) cat("  Decode-biv completed in", round(decode_D_time, 1), "s\n")
  if (is.finite(decode_F_time)) cat("  Decode-biv+KDE completed in", round(decode_F_time, 1), "s\n")
  G_params <- if (!is.null(decode_D)) decode_D$params else NULL
  H_params <- if (!is.null(decode_F)) decode_F$params else NULL
  pp_post_decode_G <- if (!is.null(decode_D)) decode_D$labelling[decode_D$labelling$t >= 0, ] else NULL
  pp_post_decode_H <- if (!is.null(decode_F)) decode_F$labelling[decode_F$labelling$t >= 0, ] else NULL
  add_timing_row(
    stage = "fit_G_decode_bivariate",
    elapsed_sec = decode_D_time,
    status = if (!is.null(decode_D)) "ok" else "failed",
    detail = "elapsed from joint fit dispatch"
  )
  add_timing_row(
    stage = "fit_H_decode_bivariate_kde",
    elapsed_sec = decode_F_time,
    status = if (!is.null(decode_F)) "ok" else "failed",
    detail = "elapsed from joint fit dispatch"
  )
}

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

  biv_init_local <- apply_structural_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
  fitE_local <- tryCatch({
    fit_etas_bivariate(
      params_init = biv_init_local, realiz = pp_all_bg_local,
      windowT = windowT_fit, windowS = win_km, m0 = ETAS_M0,
      control_state_space = control_ss, treated_state_space = treated_ss,
      background_rate_var = "W",
      treated_background_zero_before = 0,
      maxit = VANILLA_MAXIT, fixed_params = NULL, trace = 0
    )
  }, error = function(e) {
    cat(sprintf("  [BW %s] Fit E error: %s\n", bw_label, e$message))
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
      background_rate_var_in = "W",
      sem_inner_iter_in = SENS_SEM_INNER_ITER,
      verbose_in = FALSE,
      label = sprintf("BW %s Fit F", bw_label)
    )
  }, error = function(e) {
    cat(sprintf("  [BW %s] Fit F error: %s\n", bw_label, e$message))
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

# Save partition maps for all schemes.
tryCatch({
  for (part_nm in names(all_partitions)) {
    part_info <- all_partitions[[part_nm]]
    tile_cols <- ifelse(part_info$processes == "treated", "#fc9272", "#deebf7")
    outfile_name <- add_file_tag(sprintf("partition_map_%s.png", part_info$label))
    outfile <- file.path(PLOT_DIR, outfile_name)
    png(outfile, width = 1400, height = 850, res = 150)
    plot(part_info$partition, col = tile_cols, border = "grey45", main = paste("Partition:", part_info$label))
    points(pp_post$x, pp_post$y, pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.35))
    legend("bottomleft", legend = c("Control", "Treated"),
           fill = c("#deebf7", "#fc9272"), bty = "n")
    dev.off()
    cat(sprintf("  Saved %s\n", outfile_name))
  }
}, error = function(e) cat("  Partition plotting error:", e$message, "\n"))

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
  biv_init_p <- apply_structural_init_biv(init_bivariate_from_independent(A_ctrl, A_treat))
  fitE_p <- tryCatch({
    fit_etas_bivariate(
      params_init = biv_init_p, realiz = pp_all_p,
      windowT = windowT_p, windowS = win_km, m0 = ETAS_M0,
      control_state_space = p_ctrl_ss, treated_state_space = p_treat_ss,
      background_rate_var = "W",
      treated_background_zero_before = 0,
      maxit = VANILLA_MAXIT, fixed_params = NULL, trace = 0
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
      background_rate_var_in = "W",
      sem_inner_iter_in = SENS_SEM_INNER_ITER,
      verbose_in = FALSE,
      label = sprintf("Partition %s Fit F", label)
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
  if (identical(job$type, "bandwidth")) {
    out <- run_kde_bandwidth_fit(job$payload)
  } else if (identical(job$type, "partition")) {
    out <- run_biv_for_partition(job$payload)
  } else {
    out <- NULL
  }
  list(type = job$type, id = job$id, out = out)
}

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

kde_bandwidth_fits <- lapply(vapply(kde_bandwidth_specs, `[[`, character(1), "label"), function(lbl) {
  idx <- which(vapply(sens_out, function(z) identical(z$type, "bandwidth") && identical(z$id, lbl), logical(1)))
  if (length(idx) == 0) return(NULL)
  sens_out[[idx[1]]]$out
})
names(kde_bandwidth_fits) <- vapply(kde_bandwidth_specs, `[[`, character(1), "label")

partition_results <- lapply(sapply(all_partitions, `[[`, "label"), function(lbl) {
  idx <- which(vapply(sens_out, function(z) identical(z$type, "partition") && identical(z$id, lbl), logical(1)))
  if (length(idx) == 0) return(NULL)
  sens_out[[idx[1]]]$out
})
names(partition_results) <- sapply(all_partitions, `[[`, "label")
rm(sens_out)
invisible(gc(verbose = FALSE))

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
        sem_post_file <- add_file_tag(sprintf("pp_post_sem_%s.png", nm$label))
        ggsave(file.path(PLOT_DIR, sem_post_file),
               p_s, width = 12, height = 7, dpi = 150)
        cat(sprintf("  Saved %s\n", sem_post_file))
      }
      p_f <- plot_flips(nm$res)
      sem_flips_file <- add_file_tag(sprintf("sem_flips_%s.png", nm$label))
      ggsave(file.path(PLOT_DIR, sem_flips_file),
             p_f, width = 8, height = 5, dpi = 150)
      cat(sprintf("  Saved %s\n", sem_flips_file))
    }
  }, error = function(e) cat(sprintf("  %s plot error: %s\n", nm$title, e$message)))
}

# ============================================================================
# 6. ATE estimation (post-treatment horizon)
# ============================================================================
cat(sprintf("\n--- Step 6: ATE estimation (%d-day simulation horizon) ---\n", ATE_WINDOW_DAYS))

pp_post_sem_D <- if (!is.null(semD)) semD$adaptive$adaptive_labelling else NULL
pp_post_sem_D <- if (!is.null(pp_post_sem_D)) pp_post_sem_D[pp_post_sem_D$t >= 0, ] else pp_post
pp_post_sem_F <- if (!is.null(semF)) semF$adaptive$adaptive_labelling else NULL
pp_post_sem_F <- if (!is.null(pp_post_sem_F)) pp_post_sem_F[pp_post_sem_F$t >= 0, ] else pp_post

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

B_marginals <- extract_marginals(B_params)
D_marginals <- extract_marginals(D_params)
E_marginals <- extract_marginals(E_params)
F_marginals <- extract_marginals(F_params)
G_marginals <- extract_marginals(G_params)
H_marginals <- extract_marginals(H_params)

ate_estim_fast <- function(ctrl_pp, treat_pp, observed_data, label,
                           n_tiles_used = partition$n,
                           treated_idx_used = treated_idx,
                           quiet = FALSE) {
  if (!quiet) cat(sprintf("  Computing ATE for %s...\n", label))
  tryCatch({
    pre_history <- pp_pre[, c("x", "y", "t", "mag"), drop = FALSE]
    run_one_sim <- function(s) {
      c_sim <- sim_etas(ctrl_pp, windowT_ate, windowS = win_km,
                        m0 = ETAS_M0, beta_gr = BETA_GR,
                        filtration = pre_history)
      t_sim <- sim_etas(treat_pp, windowT_ate, windowS = win_km,
                        m0 = ETAS_M0, beta_gr = BETA_GR,
                        filtration = pre_history)
      c(c_count = length(c_sim$t), t_count = length(t_sim$t))
    }
    sim_results <- if (N_CORES > 1 && ATE_N_SIMS > 1) {
      run_parallel(
        as.list(seq_len(ATE_N_SIMS)), run_one_sim,
        cores = min(ATE_SIM_CORES, ATE_N_SIMS),
        label = "ate-sim"
      )
    } else {
      lapply(seq_len(ATE_N_SIMS), run_one_sim)
    }
    c_counts <- vapply(sim_results, function(z) as.numeric(z[["c_count"]]), numeric(1))
    t_counts <- vapply(sim_results, function(z) as.numeric(z[["t_count"]]), numeric(1))
    total_effect <- t_counts - c_counts
    all_nothing_sim <- data.frame(
      c_total = c_counts,
      t_total = t_counts,
      total_effect = total_effect,
      c_mean = c_counts / n_tiles_used,
      t_mean = t_counts / n_tiles_used,
      ATE    = total_effect / n_tiles_used
    )
    n_ctrl_loc <- sum(observed_data$location_process == "control", na.rm = TRUE)
    n_treat_loc <- sum(observed_data$location_process == "treated", na.rm = TRUE)
    n_ctrl_tiles <- sum(!treated_idx_used)
    n_treat_tiles <- sum(treated_idx_used)
    ATE_naive <- NA_real_
    ATE_spillover <- NA_real_
    total_naive <- NA_real_
    total_spillover <- NA_real_
    if (n_ctrl_tiles > 0 && n_treat_tiles > 0) {
      ATE_naive <- n_treat_loc / n_treat_tiles - n_ctrl_loc / n_ctrl_tiles
      total_naive <- ATE_naive * n_tiles_used
      ATE_spillover <- if ("inferred_process" %in% names(observed_data)) {
        n_ctrl_loc / n_ctrl_tiles -
          sum(observed_data$inferred_process == "control" &
              observed_data$location_process == "control", na.rm = TRUE) / n_ctrl_tiles
      } else { 0 }
      total_spillover <- ATE_spillover * n_tiles_used
    }
    analytic <- ATE_analytic_etas(ctrl_pp, treat_pp,
                                  windowT = windowT_ate, n_tiles = n_tiles_used,
                                  beta_gr = BETA_GR, m0 = ETAS_M0)
    list(all_nothing_sim = all_nothing_sim,
         ATE_naive = ATE_naive, ATE_spillover = ATE_spillover,
         total_naive = total_naive, total_spillover = total_spillover,
         n_tiles_used = n_tiles_used,
         treated_pp = treat_pp, control_pp = ctrl_pp,
         analytic = analytic)
  }, error = function(e) { cat("    Error:", e$message, "\n"); NULL })
}

t_ate_B <- proc.time()[["elapsed"]]
ate_B <- ate_estim_fast(B_marginals$ctrl, B_marginals$treat, pp_post,
                        "Fit B (naive bivariate)")
ate_B_elapsed <- proc.time()[["elapsed"]] - t_ate_B
add_timing_row("ate_B", ate_B_elapsed, if (!is.null(ate_B)) "ok" else "failed")
t_ate_D <- proc.time()[["elapsed"]]
ate_D <- ate_estim_fast(D_marginals$ctrl, D_marginals$treat, pp_post_sem_D,
                        "Fit D (SEM bivariate)")
ate_D_elapsed <- proc.time()[["elapsed"]] - t_ate_D
add_timing_row("ate_D", ate_D_elapsed, if (!is.null(ate_D)) "ok" else "failed")
t_ate_E <- proc.time()[["elapsed"]]
ate_E <- ate_estim_fast(E_marginals$ctrl, E_marginals$treat, pp_post_bg,
                        "Fit E (naive biv+KDE)")
ate_E_elapsed <- proc.time()[["elapsed"]] - t_ate_E
add_timing_row("ate_E", ate_E_elapsed, if (!is.null(ate_E)) "ok" else "failed")
t_ate_F <- proc.time()[["elapsed"]]
ate_F <- ate_estim_fast(F_marginals$ctrl, F_marginals$treat, pp_post_sem_F,
                        "Fit F (SEM biv+KDE)")
ate_F_elapsed <- proc.time()[["elapsed"]] - t_ate_F
add_timing_row("ate_F", ate_F_elapsed, if (!is.null(ate_F)) "ok" else "failed")
ate_G <- if (RUN_DECODE && !is.null(G_marginals) && !is.null(pp_post_decode_G)) {
  t_ate_G <- proc.time()[["elapsed"]]
  ate_estim_fast(G_marginals$ctrl, G_marginals$treat, pp_post_decode_G,
                 "Fit G (Decode biv)")
} else NULL
if (exists("t_ate_G", inherits = FALSE)) {
  ate_G_elapsed <- proc.time()[["elapsed"]] - t_ate_G
  add_timing_row("ate_G", ate_G_elapsed, if (!is.null(ate_G)) "ok" else "failed")
} else {
  add_timing_row("ate_G", NA_real_, "skipped", "decode disabled or unavailable")
}
ate_H <- if (RUN_DECODE && !is.null(H_marginals) && !is.null(pp_post_decode_H)) {
  t_ate_H <- proc.time()[["elapsed"]]
  ate_estim_fast(H_marginals$ctrl, H_marginals$treat, pp_post_decode_H,
                 "Fit H (Decode biv+KDE)")
} else NULL
if (exists("t_ate_H", inherits = FALSE)) {
  ate_H_elapsed <- proc.time()[["elapsed"]] - t_ate_H
  add_timing_row("ate_H", ate_H_elapsed, if (!is.null(ate_H)) "ok" else "failed")
} else {
  add_timing_row("ate_H", NA_real_, "skipped", "decode disabled or unavailable")
}

# Save a checkpoint before bootstrap so long runs retain core fit outputs
# even if bootstrap gets interrupted or OOM-killed.
cat("\n--- Step 6a checkpoint: saving pre-bootstrap results ---\n")
# Include report-facing fields (pp_data, counties, kde_info, full config) so Quarto can
# render from this checkpoint alone — same keys as final oklahoma_results.rds.
pre_boot_saved_at <- as.character(Sys.time())
results_pre_bootstrap <- list(
  fitB = list(params = B_params, loglik = B_loglik, fit = fitB, ate = ate_B),
  fitD = list(params = D_params, ctrl = D_ctrl, treat = D_treat, sem = semD, ate = ate_D),
  fitE = list(params = E_params, loglik = E_loglik, fit = fitE, ate = ate_E),
  fitF = list(params = F_params, ctrl = F_ctrl, treat = F_treat, sem = semF, ate = ate_F),
  fitG = if (RUN_DECODE) list(params = G_params, decode = decode_D, ate = ate_G) else NULL,
  fitH = if (RUN_DECODE) list(params = H_params, decode = decode_F, ate = ate_H) else NULL,
  bootstrap_ate = NULL,
  pp_data = list(pp_pre = pp_pre, pp_pre_holdout = pp_pre_holdout, pp_post = pp_post),
  kde_info = list(
    bw_sigma = as.numeric(bw_sigma), n_training = n_pre_holdout_ctrl,
    n_pre_ctrl_holdout = n_pre_holdout_ctrl,
    n_pre_holdout = nrow(pp_pre_holdout),
    n_pre_used = nrow(pp_pre),
    trigger_range_km = trigger_range_km
  ),
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
  config = list(
    ETAS_M0 = ETAS_M0, BETA_GR = BETA_GR,
    FIXED_STRUCTURAL = FIXED_STRUCTURAL,
    SEM_N_ITER = SEM_N_ITER, SEM_INNER_ITER = SEM_INNER_ITER,
    SEM_INNER_PROPS = SEM_INNER_PROPS,
    SEM_N_LABELLINGS = SEM_N_LABELLINGS,
    SEM_CHANGE_FACTOR = SEM_CHANGE_FACTOR,
    SEM_TEMPORAL_WEIGHT = SEM_TEMPORAL_WEIGHT,
    SEM_TEMPORAL_SCALE_DAYS = SEM_TEMPORAL_SCALE_DAYS,
    RUN_DECODE = RUN_DECODE,
    RUN_SENSITIVITY = RUN_SENSITIVITY,
    MEMORY_SAFE = MEMORY_SAFE,
    TRIM_SENS_OBJECTS = TRIM_SENS_OBJECTS,
    SENS_CORES = SENS_CORES,
    ATE_SIM_CORES = ATE_SIM_CORES,
    DECODE_ITER = DECODE_ITER,
    ATE_N_SIMS = ATE_N_SIMS, ATE_WINDOW_DAYS = ATE_WINDOW_DAYS,
    RUN_BOOTSTRAP_ATE = RUN_BOOTSTRAP_ATE,
    BOOT_N_REPS = BOOT_N_REPS,
    BOOT_TARGETS = BOOT_TARGETS,
    BOOT_REFIT_SCOPE = BOOT_REFIT_SCOPE,
    BOOT_SEM_INNER_ITER = BOOT_SEM_INNER_ITER,
    BOOT_OUTER_CORES = BOOT_OUTER_CORES,
    BOOT_SEED = BOOT_SEED,
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
legacy_pre_bootstrap_out_file <- file.path(LEGACY_OUT_DIR, add_file_tag("oklahoma_results_pre_bootstrap.rds"))
if (!identical(normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE),
               normalizePath(LEGACY_OUT_DIR, winslash = "/", mustWork = FALSE))) {
  invisible(file.copy(pre_bootstrap_out_file, legacy_pre_bootstrap_out_file, overwrite = TRUE))
}
rm(results_pre_bootstrap)
invisible(gc(verbose = FALSE))

# Parametric bootstrap ATEs (supports E and F currently).
bootstrap_ate <- NULL
boot_targets_run <- intersect(BOOT_TARGETS, c("E", "F"))
bootstrap_elapsed <- NA_real_
if (RUN_BOOTSTRAP_ATE && BOOT_N_REPS > 0L && length(boot_targets_run) > 0L) {
  t_bootstrap <- proc.time()[["elapsed"]]
  cat(sprintf("\n--- Step 6b: Parametric bootstrap ATEs (targets=%s, reps=%d, scope=%s, outer cores=%d, boot SEM inner=%d) ---\n",
              paste(boot_targets_run, collapse = ","), BOOT_N_REPS, BOOT_REFIT_SCOPE, BOOT_OUTER_CORES, BOOT_SEM_INNER_ITER))

  if (BOOT_REFIT_SCOPE == "full") {
    cat("  Note: full scope currently falls back to partial refit for E/F.\n")
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
  e_ctrl_seed <- E_marginals$ctrl
  e_treat_seed <- E_marginals$treat
  f_ctrl_seed <- F_marginals$ctrl
  f_treat_seed <- F_marginals$treat

  simulate_boot_data <- function(ctrl_seed, treat_seed) {
    pre_sim <- sim_etas(pre_ctrl_seed, pre_window_boot, windowS = control_ss,
                        m0 = ETAS_M0, beta_gr = BETA_GR)
    pre_df <- as_pp_df(pre_sim, "control", "control")
    history_df <- pre_df[, c("x", "y", "t", "mag"), drop = FALSE]
    ctrl_post_sim <- sim_etas(ctrl_seed, post_window_boot, windowS = control_ss,
                              m0 = ETAS_M0, beta_gr = BETA_GR, filtration = history_df)
    treat_post_sim <- sim_etas(treat_seed, post_window_boot, windowS = treated_ss,
                               m0 = ETAS_M0, beta_gr = BETA_GR, filtration = history_df)
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
      ate_total_mean = mean(ate_obj$all_nothing_sim$total_effect, na.rm = TRUE),
      ate_total_sd = stats::sd(ate_obj$all_nothing_sim$total_effect, na.rm = TRUE),
      ate_tile_mean = mean(ate_obj$all_nothing_sim$ATE, na.rm = TRUE),
      eta_ctrl = ate_obj$analytic$eta_ctrl,
      eta_treat = ate_obj$analytic$eta_treat,
      params = par_obj
    )
  }

  run_boot_rep <- function(rep_id) {
    if (!is.na(BOOT_SEED)) set.seed(BOOT_SEED + rep_id)
    out <- list(rep = rep_id)

    if ("E" %in% boot_targets_run) {
      sim_E <- simulate_boot_data(e_ctrl_seed, e_treat_seed)
      e_params_boot <- E_params
      if (BOOT_REFIT_SCOPE %in% c("partial", "full")) {
        fit_e_boot <- tryCatch({
          fit_etas_bivariate(
            params_init = E_params, realiz = sim_E$pp_all_bg_sim,
            windowT = windowT_fit, windowS = win_km, m0 = ETAS_M0,
            control_state_space = control_ss, treated_state_space = treated_ss,
            background_rate_var = "W",
            treated_background_zero_before = 0,
            maxit = VANILLA_MAXIT, fixed_params = NULL, trace = 0
          )
        }, error = function(e) NULL)
        if (!is.null(fit_e_boot) && !is.null(fit_e_boot$par)) e_params_boot <- fit_e_boot$par
      }
      e_marg_boot <- extract_marginals(e_params_boot)
      ate_e_boot <- ate_estim_fast(
        e_marg_boot$ctrl, e_marg_boot$treat, sim_E$pp_post_bg_sim,
        label = sprintf("Boot E #%d", rep_id),
        n_tiles_used = partition$n,
        treated_idx_used = treated_idx,
        quiet = TRUE
      )
      out$E <- summarize_boot(ate_e_boot, rep_id, sim_E$pre_df, sim_E$post_ctrl_df, sim_E$post_treat_df, e_params_boot)
    }

    if ("F" %in% boot_targets_run) {
      sim_F <- simulate_boot_data(f_ctrl_seed, f_treat_seed)
      f_params_boot <- F_params
      pp_post_f_boot <- sim_F$pp_post_bg_sim
      if (BOOT_REFIT_SCOPE %in% c("partial", "full")) {
        sem_boot <- run_sem_fit(
          pp_data_in = sim_F$pp_all_bg_sim,
          partition_in = partition,
          partition_processes_in = partition_processes,
          state_spaces_in = state_spaces,
          init_params_in = F_params,
          background_rate_var_in = "W",
          sem_inner_iter_in = BOOT_SEM_INNER_ITER,
          verbose_in = FALSE,
          label = sprintf("Boot F #%d", rep_id)
        )
        if (!is.null(sem_boot) && !is.null(sem_boot$etas_bivariate_params)) {
          f_params_boot <- sem_boot$etas_bivariate_params
        }
        if (!is.null(sem_boot) && !is.null(sem_boot$adaptive$adaptive_labelling)) {
          pp_post_f_boot <- sem_boot$adaptive$adaptive_labelling
          pp_post_f_boot <- pp_post_f_boot[pp_post_f_boot$t >= 0, , drop = FALSE]
        }
      }
      f_marg_boot <- extract_marginals(f_params_boot)
      ate_f_boot <- ate_estim_fast(
        f_marg_boot$ctrl, f_marg_boot$treat, pp_post_f_boot,
        label = sprintf("Boot F #%d", rep_id),
        n_tiles_used = partition$n,
        treated_idx_used = treated_idx,
        quiet = TRUE
      )
      out$F <- summarize_boot(ate_f_boot, rep_id, sim_F$pre_df, sim_F$post_ctrl_df, sim_F$post_treat_df, f_params_boot)
    }
    rm_vars <- intersect(
      c("sim_E", "fit_e_boot", "e_params_boot", "e_marg_boot", "ate_e_boot",
        "sim_F", "sem_boot", "f_params_boot", "f_marg_boot", "ate_f_boot", "pp_post_f_boot"),
      ls()
    )
    if (length(rm_vars) > 0L) rm(list = rm_vars)
    invisible(gc(verbose = FALSE))
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
rm(boot_results)
invisible(gc(verbose = FALSE))

  if (!is.null(bootstrap_ate$fit_E)) {
    bE <- bootstrap_ate$fit_E$replicate_summary
    if (nrow(bE) > 0) {
      cat(sprintf("  Bootstrap E complete: success=%d, fail=%d, mean total ATE=%.1f, SD(rep means)=%.1f\n",
                  nrow(bE), bootstrap_ate$fit_E$n_fail,
                  mean(bE$ate_total_mean, na.rm = TRUE),
                  stats::sd(bE$ate_total_mean, na.rm = TRUE)))
    } else {
      cat(sprintf("  Bootstrap E complete: success=%d, fail=%d\n",
                  bootstrap_ate$fit_E$n_success, bootstrap_ate$fit_E$n_fail))
    }
  }
  if (!is.null(bootstrap_ate$fit_F)) {
    bF <- bootstrap_ate$fit_F$replicate_summary
    if (nrow(bF) > 0) {
      cat(sprintf("  Bootstrap F complete: success=%d, fail=%d, mean total ATE=%.1f, SD(rep means)=%.1f\n",
                  nrow(bF), bootstrap_ate$fit_F$n_fail,
                  mean(bF$ate_total_mean, na.rm = TRUE),
                  stats::sd(bF$ate_total_mean, na.rm = TRUE)))
    } else {
      cat(sprintf("  Bootstrap F complete: success=%d, fail=%d\n",
                  bootstrap_ate$fit_F$n_success, bootstrap_ate$fit_F$n_fail))
    }
  }

  # Bootstrap reached completion, so drop pre-bootstrap checkpoint.
  if (file.exists(pre_bootstrap_out_file)) {
    unlink(pre_bootstrap_out_file, force = TRUE)
    cat(sprintf("Deleted pre-bootstrap checkpoint: %s\n", pre_bootstrap_out_file))
  }
  if (exists("legacy_pre_bootstrap_out_file", inherits = FALSE) &&
      file.exists(legacy_pre_bootstrap_out_file)) {
    unlink(legacy_pre_bootstrap_out_file, force = TRUE)
    cat(sprintf("Deleted legacy pre-bootstrap checkpoint: %s\n", legacy_pre_bootstrap_out_file))
  }
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
                     n_tiles_used = partition$n,
                     treated_idx_used = treated_idx),
      error = function(e) NULL
    ),
    ate_F = tryCatch(
      ate_estim_fast(fm$ctrl, fm$treat, kf$pp_post_sem,
                     sprintf("BW %s SEM biv+KDE", kf$label),
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
                     n_tiles_used = pr$n_tiles,
                     treated_idx_used = part_info$treated_idx),
      error = function(e) NULL),
    ate_F = tryCatch(
      ate_estim_fast(fm$ctrl, fm$treat, pp_post_sem_part,
                     sprintf("%s SEM biv+KDE", pr$label),
                     n_tiles_used = pr$n_tiles,
                     treated_idx_used = part_info$treated_idx),
      error = function(e) NULL)
  )
})

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
  "", "B.ctrl", "B.treat", "D.ctrl", "D.treat",
  "E.ctrl", "E.treat", "F.ctrl", "F.treat"))
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
for (info in list(list(par = B_params, lab = "B"),
                  list(par = D_params, lab = "D"),
                  list(par = E_params, lab = "E"),
                  list(par = F_params, lab = "F"))) {
  cat(sprintf("  %s: A_01=%.4f  A_10=%.4f  alpha_m_01=%.4f  alpha_m_10=%.4f\n",
    info$lab, info$par[["A_01"]], info$par[["A_10"]],
    info$par[["alpha_m_01"]], info$par[["alpha_m_10"]]))
}

cat("\n")
cat("===========================================================================\n")
cat(sprintf("              ATE COMPARISON (%d-day horizon)\n", ATE_WINDOW_DAYS))
cat("===========================================================================\n\n")
cat(sprintf("%-30s  %12s  %12s  %12s  %12s\n",
            "", "Total Effect", "SD(total)", "eta_ctrl", "eta_treat"))
cat(paste(rep("-", 90), collapse = ""), "\n")

ate_print_list <- list(
  list(ate = ate_B, lab = "B: Naive Biv (county)"),
  list(ate = ate_D, lab = "D: SEM Biv (county)"),
  list(ate = ate_E, lab = "E: Naive Biv+KDE (county)"),
  list(ate = ate_F, lab = "F: SEM Biv+KDE (county)"))

for (pname in names(ate_partitions)) {
  ap <- ate_partitions[[pname]]
  if (!is.null(ap$ate_E))
    ate_print_list[[length(ate_print_list) + 1]] <- list(ate = ap$ate_E, lab = sprintf("E: %s", pname))
  if (!is.null(ap$ate_F))
    ate_print_list[[length(ate_print_list) + 1]] <- list(ate = ap$ate_F, lab = sprintf("F: %s", pname))
}

for (nm in ate_print_list) {
  if (!is.null(nm$ate)) {
    total <- mean(nm$ate$all_nothing_sim$total_effect, na.rm = TRUE)
    total_sd <- stats::sd(nm$ate$all_nothing_sim$total_effect, na.rm = TRUE)
    eta_c <- if (!is.null(nm$ate$analytic)) nm$ate$analytic$eta_ctrl else NA
    eta_t <- if (!is.null(nm$ate$analytic)) nm$ate$analytic$eta_treat else NA
    cat(sprintf("%-30s  %12.0f  %12.0f  %12.3f  %12.3f\n",
      nm$lab, total, total_sd, eta_c, eta_t))
  } else {
    cat(sprintf("%-30s  %12s  %12s  %12s  %12s\n",
      nm$lab, "FAILED", "FAILED", "FAILED", "FAILED"))
  }
}
cat(sprintf("\nTotal Effect: expected all-or-nothing total earthquakes over %d days in Oklahoma.\n",
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

results <- list(
  fitB = list(params = B_params, loglik = B_loglik, fit = fitB, ate = ate_B),
  fitD = list(params = D_params, ctrl = D_ctrl, treat = D_treat, sem = semD, ate = ate_D, louis = louis_D),
  fitE = list(params = E_params, loglik = E_loglik, fit = fitE, ate = ate_E),
  fitF = list(params = F_params, ctrl = F_ctrl, treat = F_treat, sem = semF, ate = ate_F, louis = louis_F),
  fitG = if (RUN_DECODE) list(params = G_params, decode = decode_D, ate = ate_G) else NULL,
  fitH = if (RUN_DECODE) list(params = H_params, decode = decode_F, ate = ate_H) else NULL,
  bootstrap_ate = bootstrap_ate,
  partition_results = partition_results,
  ate_partitions = ate_partitions,
  kde_bandwidth_sensitivity = kde_bandwidth_sensitivity,
  kde_info = list(bw_sigma = as.numeric(bw_sigma), n_training = n_pre_holdout_ctrl,
                  n_pre_ctrl_holdout = n_pre_holdout_ctrl,
                  n_pre_holdout = nrow(pp_pre_holdout),
                  n_pre_used = nrow(pp_pre),
                  trigger_range_km = trigger_range_km),
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
  config = list(
    ETAS_M0 = ETAS_M0, BETA_GR = BETA_GR,
    FIXED_STRUCTURAL = FIXED_STRUCTURAL,
    SEM_N_ITER = SEM_N_ITER, SEM_INNER_ITER = SEM_INNER_ITER,
    SEM_INNER_PROPS = SEM_INNER_PROPS,
    SEM_N_LABELLINGS = SEM_N_LABELLINGS,
    SEM_CHANGE_FACTOR = SEM_CHANGE_FACTOR,
    SEM_TEMPORAL_WEIGHT = SEM_TEMPORAL_WEIGHT,
    SEM_TEMPORAL_SCALE_DAYS = SEM_TEMPORAL_SCALE_DAYS,
    RUN_DECODE = RUN_DECODE,
    RUN_SENSITIVITY = RUN_SENSITIVITY,
    MEMORY_SAFE = MEMORY_SAFE,
    TRIM_SENS_OBJECTS = TRIM_SENS_OBJECTS,
    SENS_CORES = SENS_CORES,
    ATE_SIM_CORES = ATE_SIM_CORES,
    DECODE_ITER = DECODE_ITER,
    ATE_N_SIMS = ATE_N_SIMS, ATE_WINDOW_DAYS = ATE_WINDOW_DAYS,
    RUN_BOOTSTRAP_ATE = RUN_BOOTSTRAP_ATE,
    BOOT_N_REPS = BOOT_N_REPS,
    BOOT_TARGETS = BOOT_TARGETS,
    BOOT_REFIT_SCOPE = BOOT_REFIT_SCOPE,
    BOOT_SEM_INNER_ITER = BOOT_SEM_INNER_ITER,
    BOOT_OUTER_CORES = BOOT_OUTER_CORES,
    BOOT_SEED = BOOT_SEED,
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
cat(sprintf("Plots saved to:   %s\n", PLOT_DIR))

# Mirror latest results into legacy output path so older workflows stay fresh.
legacy_out_file <- file.path(LEGACY_OUT_DIR, add_file_tag("oklahoma_results.rds"))
if (!identical(normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE),
               normalizePath(LEGACY_OUT_DIR, winslash = "/", mustWork = FALSE))) {
  copied <- file.copy(out_file, legacy_out_file, overwrite = TRUE)
  if (copied) {
    cat(sprintf("Results mirrored to legacy path: %s\n", legacy_out_file))
  } else {
    cat(sprintf("Warning: failed to mirror results to legacy path: %s\n", legacy_out_file))
  }
}

# Keep report outputs synchronized with the newest results.
report_file <- file.path(SCRIPT_DIR, "oklahoma_report.qmd")
if (file.exists(report_file) && length(REPORT_FORMATS) > 0L) {
  quarto_bin <- Sys.which("quarto")
  if (nzchar(quarto_bin)) {
    cat(sprintf("\n--- Rendering report (%s) ---\n", paste(REPORT_FORMATS, collapse = ", ")))
    old_wd <- getwd()
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
    if (length(render_errors) > 0L) {
      cat("Report render failed:\n")
      cat(paste0("  - ", render_errors, collapse = "\n"), "\n")
    } else {
      cat("Report render complete.\n")
      sync_stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
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
                                !is.na(report_html_rendered) && nzchar(FILE_TAG)) {
        file.path(SCRIPT_DIR, add_file_tag("oklahoma_report.html"))
      } else {
        report_html_rendered
      }
      report_pdf_target <- if (is.character(report_pdf_rendered) &&
                               !is.na(report_pdf_rendered) && nzchar(FILE_TAG)) {
        file.path(SCRIPT_DIR, add_file_tag("oklahoma_report.pdf"))
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
      if (!identical(normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE),
                     normalizePath(LEGACY_OUT_DIR, winslash = "/", mustWork = FALSE))) {
        writeLines(stamp_lines, file.path(LEGACY_OUT_DIR, add_file_tag("last_run_sync_stamp.txt")))
      }
      cat("Sync stamp written for Google Drive change detection.\n")
    }
  } else {
    cat("Quarto not found in PATH; skipping report render.\n")
  }
} else if (file.exists(report_file)) {
  cat("Report render skipped (OK_REPORT_FORMATS empty).\n")
}

cat("\n=== Oklahoma Analysis Complete ===\n")
