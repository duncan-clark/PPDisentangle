#!/usr/bin/env Rscript
# Simulation study for PPDisentangle.
# Usage:
#   Rscript sim_study.R                  # local (auto-detect cores)
#   Rscript sim_study.R --test --sims 2  # quick smoke test
#   Rscript sim_study.R --small          # local, reduced
#   sbatch run_nesi.sh --sims 100        # NeSI cluster
#
# Output (canonical): ../PPDisentangle-output/sim_study/{JOB_ID}.rds

prepend_user_lib_paths <- function() {
  user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
  if (!nzchar(user_lib)) return(invisible(NULL))
  extra_libs <- strsplit(user_lib, .Platform$path.sep, fixed = TRUE)[[1]]
  extra_libs <- extra_libs[nzchar(extra_libs)]
  if (length(extra_libs) > 0L) {
    .libPaths(c(extra_libs, .libPaths()))
  }
}

wait_for_namespace <- function(pkg, timeout_s = 180L, sleep_s = 5L) {
  prepend_user_lib_paths()
  lock_dir <- file.path(Sys.getenv("R_LIBS_USER", unset = ""), paste0("00LOCK-", pkg))
  start_time <- Sys.time()
  repeat {
    prepend_user_lib_paths()
    if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
    waited_s <- as.integer(difftime(Sys.time(), start_time, units = "secs"))
    if (waited_s >= timeout_s) {
      stop(sprintf(
        "Package '%s' not visible after %ds. .libPaths()=%s",
        pkg, timeout_s, paste(.libPaths(), collapse = " | ")
      ))
    }
    if (nzchar(lock_dir) && dir.exists(lock_dir)) {
      message(sprintf("Waiting for package '%s' while lock exists (%s); waited %ds...", pkg, lock_dir, waited_s))
    } else {
      message(sprintf("Waiting for package '%s' to become visible; waited %ds...", pkg, waited_s))
    }
    Sys.sleep(sleep_s)
  }
}

wait_for_namespace("PPDisentangle")
library(PPDisentangle)
library(spatstat)
library(ggplot2)
library(dplyr)
library(data.table)
library(parallel)
library(doParallel)
normalize_hawkes_kernel <- getFromNamespace("normalize_hawkes_kernel", "PPDisentangle")
normalize_hawkes_spatial_kernel <- getFromNamespace("normalize_hawkes_spatial_kernel", "PPDisentangle")
as_hawkes_params <- getFromNamespace("as_hawkes_params", "PPDisentangle")

time_start_global <- proc.time()[3]

args <- commandArgs(trailingOnly = TRUE)
SMALL <- "--small" %in% args
FORCE_CLUSTER <- "--cluster" %in% args
TEST <- "--test" %in% args
MODE_ENV <- tolower(trimws(Sys.getenv("PP_MODE", "")))
LONG_MODE <- MODE_ENV %in% c("long", "full", "big")

sims_arg <- grep("^--sims$", args)
N_SIMS_ARG <- if (length(sims_arg) > 0 && length(args) >= sims_arg + 1L)
  as.numeric(args[sims_arg + 1L]) else NULL

ON_CLUSTER <- nzchar(Sys.getenv("SLURM_JOB_ID")) || FORCE_CLUSTER

resolve_save_dirs <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  script_dir <- if (length(file_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    file.path(getwd(), "inst", "sim_study")
  }
  repo_dir <- if (basename(script_dir) == "sim_study" &&
                  basename(dirname(script_dir)) == "inst") {
    normalizePath(dirname(dirname(script_dir)), winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
  list(
    canonical = PPDisentangle::pp_output_path("sim_study", repo_root = repo_dir),
    script_dir = script_dir,
    repo_dir = repo_dir
  )
}

OMEGA <- c(0, 100, 0, 100)
TREATMENT_TIME <- 10
END_TIME_BASE <- 110
NX <- 10
NY <- 10

local_core_default <- max(1L, parallel::detectCores() - 1L)
if (TEST) {
  n_test <- if (!is.null(N_SIMS_ARG) && is.finite(N_SIMS_ARG)) N_SIMS_ARG else 2L
  if (ON_CLUSTER) {
    slurm_cores_test <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", as.character(n_test))))
    if (!is.finite(slurm_cores_test) || is.na(slurm_cores_test) || slurm_cores_test < 1L) slurm_cores_test <- n_test
    N_CORES <- max(1L, min(as.integer(n_test), slurm_cores_test))
  } else {
    N_CORES <- max(1L, min(as.integer(n_test), parallel::detectCores()))
  }
  SIM_SIZE <- n_test
  N_SIMS <- n_test
  N_TAU_SIMS <- 3
  N_TAU_I <- 3
  N_TAU_I_TRUE <- 5
  N_PROPOSALS <- 5
  SEM_EM_ADAPTIVE_ITER <- 200
  SEM_N_ITER <- 1
  SEM_N_LABELLINGS <- 5
} else if (ON_CLUSTER) {
  N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "100"))
  SIM_SIZE <- 100
  N_SIMS <- 100
  N_TAU_SIMS <- 10
  N_TAU_I <- 10
  N_TAU_I_TRUE <- 100
  N_PROPOSALS <- 10
  # Restore pre-bootstrap baseline SEM settings for cluster stability/accuracy.
  SEM_EM_ADAPTIVE_ITER <- 2000
  SEM_N_ITER <- 10
  SEM_N_LABELLINGS <- 10
} else if (SMALL) {
  N_CORES <- min(8L, local_core_default)
  SIM_SIZE <- N_CORES
  N_SIMS <- 10
  N_TAU_SIMS <- 10
  N_TAU_I <- 10
  N_TAU_I_TRUE <- 10
  N_PROPOSALS <- 5
  SEM_EM_ADAPTIVE_ITER <- 100
  SEM_N_ITER <- 10
  SEM_N_LABELLINGS <- 10
} else {
  N_CORES <- local_core_default
  SIM_SIZE <- N_CORES
  N_SIMS <- 10
  N_TAU_SIMS <- 10
  N_TAU_I <- 10
  N_TAU_I_TRUE <- 10
  N_PROPOSALS <- 100
  SEM_EM_ADAPTIVE_ITER <- 1000
  SEM_N_ITER <- 3
  SEM_N_LABELLINGS <- max(10, N_PROPOSALS %/% 10)
}

save_dirs <- resolve_save_dirs()
SAVE_DIR <- save_dirs$canonical
SCRIPT_DIR <- save_dirs$script_dir

if (!is.null(N_SIMS_ARG) && is.finite(N_SIMS_ARG)) {
  N_SIMS <- N_SIMS_ARG
  SIM_SIZE <- N_SIMS_ARG
}

if (ON_CLUSTER && !TEST) {
  slurm_cores <- suppressWarnings(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", NA_character_)))
  if (is.finite(slurm_cores) && slurm_cores >= 1) {
    if (!isTRUE(all.equal(N_CORES, slurm_cores))) {
      message(sprintf(
        "[sim_study] Using SLURM_CPUS_PER_TASK for worker budget: sims=%s, cores=%s, slurm=%s",
        as.character(SIM_SIZE), as.character(N_CORES), as.character(slurm_cores)
      ))
    }
    N_CORES <- slurm_cores
  }
}

env_int <- function(name, default_value, min_value = 1L) {
  raw <- Sys.getenv(name, "")
  if (!nzchar(raw)) return(as.integer(default_value))
  val <- suppressWarnings(as.integer(raw))
  if (!is.finite(val) || is.na(val)) return(as.integer(default_value))
  as.integer(max(min_value, val))
}

env_num <- function(name, default_value, min_value = NULL) {
  raw <- Sys.getenv(name, "")
  if (!nzchar(raw)) return(as.numeric(default_value))
  val <- suppressWarnings(as.numeric(raw))
  if (!is.finite(val) || is.na(val)) return(as.numeric(default_value))
  if (!is.null(min_value)) val <- max(min_value, val)
  as.numeric(val)
}

env_chr <- function(name, default_value) {
  raw <- trimws(Sys.getenv(name, ""))
  if (!nzchar(raw)) return(as.character(default_value))
  as.character(raw)
}

# Optional runtime overrides for quick/full profiles without code edits.
N_PROPOSALS <- env_int("PP_LABEL_PROPOSALS", N_PROPOSALS, 1L)
SEM_EM_ADAPTIVE_ITER <- env_int("PP_SEM_INNER_ITER", SEM_EM_ADAPTIVE_ITER, 1L)
SEM_N_ITER <- env_int("PP_SEM_OUTER_ITER", SEM_N_ITER, 1L)
SEM_N_LABELLINGS <- env_int("PP_SEM_N_LABELLINGS", SEM_N_LABELLINGS, 0L)
POST_TIME_MULTIPLIER <- env_num("PP_POST_TIME_MULTIPLIER", 1, min_value = 0.1)
END_TIME <- TREATMENT_TIME + POST_TIME_MULTIPLIER * (END_TIME_BASE - TREATMENT_TIME)
if (!is.finite(END_TIME) || END_TIME <= TREATMENT_TIME) {
  stop("Computed END_TIME is invalid; check PP_POST_TIME_MULTIPLIER.")
}
TIME_INT <- END_TIME - TREATMENT_TIME

BASE_SEED <- env_int("PP_BASE_SEED", 123L, 1L)
stage_seed <- function(stage_offset, sim_id = 0L, extra = 0L) {
  as.integer(BASE_SEED + as.integer(stage_offset) * 100000L + as.integer(sim_id) * 1000L + as.integer(extra))
}

# SEM/adaptive configuration
SEM_PARAM_UPDATE_CADENCE <- 10L
SEM_PROPOSAL_UPDATE_CADENCE <- 1L
SEM_PARAM_REFIT_CADENCE <- 1L
SEM_PARAM_UPDATE_CADENCE <- env_int("PP_SEM_PARAM_UPDATE_CADENCE", SEM_PARAM_UPDATE_CADENCE, 1L)
SEM_PROPOSAL_UPDATE_CADENCE <- env_int("PP_SEM_PROPOSAL_UPDATE_CADENCE", SEM_PROPOSAL_UPDATE_CADENCE, 1L)
SEM_PARAM_REFIT_CADENCE <- env_int("PP_SEM_PARAM_REFIT_CADENCE", SEM_PARAM_REFIT_CADENCE, 1L)
SEM_N_PROPS <- 20L
SEM_N_PROPS <- env_int("PP_SEM_N_PROPS", SEM_N_PROPS, 1L)
SEM_CHANGE_FACTOR <- env_num("PP_SEM_CHANGE_FACTOR", 0.005, min_value = 1e-6)
SEM_STALENESS_TRIGGER_EVERY <- 10L
SEM_STALENESS_TRIGGER_EVERY <- env_int("PP_SEM_STALENESS_TRIGGER_EVERY", SEM_STALENESS_TRIGGER_EVERY, 1L)
SEM_INCLUDE_STARTING <- tolower(Sys.getenv("PP_SEM_INCLUDE_STARTING", "true")) %in% c("1", "true", "yes", "y")
SEM_UPDATE_STARTING <- tolower(Sys.getenv("PP_SEM_UPDATE_STARTING", "true")) %in% c("1", "true", "yes", "y")
SEM_UPDATE_CONTROL_PARAMS <- tolower(Sys.getenv("PP_SEM_UPDATE_CONTROL_PARAMS", "false")) %in% c("1", "true", "yes", "y")

TREAT_PROP <- 0.5
TREATMENT_ASSIGNMENT <- env_chr("PP_TREATMENT_ASSIGNMENT", "random")
MAX_TIME   <- 10000 * (END_TIME * OMEGA[2] * OMEGA[4] / 1e6)

SCENARIO_ID <- env_chr("PP_SCENARIO_ID", "baseline_hawkes_exp")
SIM_KERNEL <- normalize_hawkes_kernel(env_chr("PP_SIM_KERNEL", "exponential"))
FIT_KERNEL <- normalize_hawkes_kernel(env_chr("PP_FIT_KERNEL", SIM_KERNEL))
SIM_SPATIAL_KERNEL <- normalize_hawkes_spatial_kernel(env_chr("PP_SIM_SPATIAL_KERNEL", "exponential"))
FIT_SPATIAL_KERNEL <- normalize_hawkes_spatial_kernel(env_chr("PP_FIT_SPATIAL_KERNEL", SIM_SPATIAL_KERNEL))
CONTROL_K <- env_num("PP_CONTROL_K", 0.8, min_value = 0)
TREATED_K <- env_num("PP_TREATED_K", 0.2, min_value = 0)
if (CONTROL_K >= 1 || TREATED_K >= 1) {
  stop("PP_CONTROL_K and PP_TREATED_K must be in [0, 1).")
}
BASE_MU <- env_num("PP_BASE_MU", 8, min_value = 1e-8)
HAWKES_ALPHA <- env_num("PP_HAWKES_ALPHA", 0.01, min_value = 1e-12)
HAWKES_BETA <- env_num("PP_HAWKES_BETA", 10, min_value = 1e-12)
HAWKES_POWER_C <- env_num("PP_HAWKES_POWER_C", 0.1, min_value = 1e-12)
HAWKES_POWER_P <- env_num("PP_HAWKES_POWER_P", 2, min_value = 1.000001)
HAWKES_SPATIAL_POWER_Q <- env_num("PP_HAWKES_SPATIAL_POWER_Q", 2, min_value = 1.500001)
compute_hawkes_temporal_trunc <- function(c_param, p_param, rel_level = 0.05) {
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
hawkes_t_trunc_raw <- Sys.getenv("PP_HAWKES_T_TRUNC", "")
HAWKES_T_TRUNC_USER <- if (nzchar(hawkes_t_trunc_raw)) {
  suppressWarnings(as.numeric(hawkes_t_trunc_raw))
} else {
  NA_real_
}
HAWKES_T_TRUNC_REL <- suppressWarnings(as.numeric(Sys.getenv("PP_HAWKES_T_TRUNC_REL", "0.05")))
if (!is.finite(HAWKES_T_TRUNC_REL) || is.na(HAWKES_T_TRUNC_REL) ||
    HAWKES_T_TRUNC_REL <= 0 || HAWKES_T_TRUNC_REL >= 1) {
  HAWKES_T_TRUNC_REL <- 0.05
}
resolve_hawkes_t_trunc <- function(fit_kernel) {
  if (!identical(normalize_hawkes_kernel(fit_kernel), "power_law")) {
    return(list(value = NULL, source = "disabled"))
  }
  if (is.finite(HAWKES_T_TRUNC_USER) && !is.na(HAWKES_T_TRUNC_USER) && HAWKES_T_TRUNC_USER > 0) {
    return(list(value = as.numeric(HAWKES_T_TRUNC_USER), source = "env"))
  }
  trunc <- compute_hawkes_temporal_trunc(HAWKES_POWER_C, HAWKES_POWER_P, HAWKES_T_TRUNC_REL)
  if (is.null(trunc)) {
    return(list(value = NULL, source = "none"))
  }
  list(
    value = trunc,
    source = sprintf(
      "auto(c=%.4g,p=%.4g,rel=%.3f)",
      HAWKES_POWER_C, HAWKES_POWER_P, HAWKES_T_TRUNC_REL
    )
  )
}
hawkes_t_trunc_cfg <- resolve_hawkes_t_trunc(FIT_KERNEL)
HAWKES_T_TRUNC <- hawkes_t_trunc_cfg$value
HAWKES_T_TRUNC_SOURCE <- hawkes_t_trunc_cfg$source
MU_SCALE <- env_num("PP_MU_SCALE", 1, min_value = 1e-8)
TARGET_POINTS <- env_num("PP_TARGET_POINTS", NA_real_, min_value = 1)
RUN_DECAY_VALIDATION <- tolower(Sys.getenv("PP_DECAY_VALIDATION", "true")) %in% c("1", "true", "yes", "y")
DECAY_VALIDATION_REPS <- env_int("PP_DECAY_REPS", if (TEST) 20L else 2000L, 1L)
DECAY_ANNULUS_WIDTH <- env_num("PP_DECAY_ANNULUS_WIDTH", 1, min_value = 1e-8)
DECAY_FLIP_CELL <- env_int("PP_DECAY_FLIP_CELL", NA_integer_, 1L)
make_hawkes_params_for_kernel <- function(mu, K, kernel, spatial_kernel = "exponential") {
  if (identical(normalize_hawkes_kernel(kernel), "power_law")) {
    out <- list(mu = mu, alpha = HAWKES_ALPHA, c = HAWKES_POWER_C, p = HAWKES_POWER_P,
                beta = HAWKES_POWER_C, K = K, kernel = "power_law")
  } else {
    out <- list(mu = mu, alpha = HAWKES_ALPHA, beta = HAWKES_BETA, K = K,
                kernel = "exponential")
  }
  out$spatial_kernel <- normalize_hawkes_spatial_kernel(spatial_kernel)
  if (identical(out$spatial_kernel, "power_law")) {
    out$spatial_q <- HAWKES_SPATIAL_POWER_Q
  }
  as_hawkes_params(out, out$kernel, out$spatial_kernel)
}
expected_points_per_mu_analytic <- function(k_control, k_treated) {
  # Expected total event count per unit mu under the sim-study design:
  # pre-treatment control Hawkes on Omega, then post-treatment mix on the partition.
  TREATMENT_TIME / (1 - k_control) +
    TIME_INT * ((1 - TREAT_PROP) / (1 - k_control) + TREAT_PROP / (1 - k_treated))
}
calibrate_mu_for_target <- function(target_points, base_mu, k_control, k_treated, mu_scale = 1) {
  if (!is.finite(target_points) || target_points <= 0) {
    return(list(mu = base_mu * mu_scale, expected_points_per_mu = NA_real_))
  }
  expected_per_mu <- expected_points_per_mu_analytic(k_control, k_treated)
  if (!is.finite(expected_per_mu) || expected_per_mu <= 0) {
    return(list(mu = base_mu * mu_scale, expected_points_per_mu = NA_real_))
  }
  list(
    mu = (target_points / expected_per_mu) * mu_scale,
    expected_points_per_mu = expected_per_mu
  )
}
mu_calibration <- calibrate_mu_for_target(TARGET_POINTS, BASE_MU, CONTROL_K, TREATED_K, MU_SCALE)
TRUE_MU <- mu_calibration$mu
EXPECTED_POINTS_PER_MU <- mu_calibration$expected_points_per_mu

# ------------------------------------------------------------------
# Logging
# ------------------------------------------------------------------
JOB_ID <- trimws(Sys.getenv("PP_OUTPUT_BASENAME", ""))
if (JOB_ID == "") JOB_ID <- Sys.getenv("SLURM_JOB_ID", "")
if (JOB_ID == "") JOB_ID <- paste0("local_", format(Sys.time(), "%Y%m%d_%H%M%S"))
OUTPUT_TAG <- trimws(Sys.getenv("PP_OUTPUT_TAG", ""))
RUN_ID <- if (nzchar(OUTPUT_TAG)) paste0(JOB_ID, "_", OUTPUT_TAG) else JOB_ID

for (d in unique(c(SAVE_DIR))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
LOG_FILE <- file.path(SAVE_DIR, paste0(RUN_ID, ".log"))
log_con <- file(LOG_FILE, open = "wt")
on.exit(tryCatch(close(log_con), error = function(e) NULL), add = TRUE)

log_msg <- function(...) {
  line <- sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(line)
  if (isOpen(log_con)) cat(line, file = log_con)
}

log_elapsed <- function(phase, elapsed_sec, n_done = NULL, n_total = NULL) {
  if (is.null(n_done) || is.null(n_total) || n_total <= 0) {
    log_msg(sprintf("%s: %.1f s", phase, elapsed_sec))
  } else {
    eta_sec <- if (n_done > 0) (elapsed_sec / n_done) * (n_total - n_done) else NA
    log_msg(sprintf("%s: %.1f s (%d/%d, ETA %.0f s)", phase, elapsed_sec, n_done, n_total, eta_sec))
  }
}

LOG_MEMORY <- FALSE
SAVE_LIGHT <- tolower(Sys.getenv("PP_SAVE_LIGHT", "true")) %in% c("1", "true", "yes", "y")
log_memory <- function(phase = "") {
  if (!LOG_MEMORY) return(invisible(NULL))
  g <- gc(verbose = FALSE)
  log_msg(sprintf("[MEM %s] used=%.0f Mb  max=%.0f Mb", phase, sum(g[, "used"]), sum(g[, "max used"])))
}

params_are_crazy <- function(control_pp, treated_pp, K_max = 0.98, mu_max = 1e5) {
  check_one <- function(pp) {
    if (is.null(pp)) return(FALSE)
    K <- if (is.list(pp)) pp$K else pp[4]
    mu <- if (is.list(pp)) pp$mu else pp[1]
    (K >= K_max) || (mu > mu_max) || (K < 0) || (K >= 1)
  }
  check_one(control_pp) || check_one(treated_pp)
}

MODE <- if (TEST) "TEST" else if (ON_CLUSTER) "CLUSTER" else if (SMALL) "SMALL" else "LOCAL"
ATE_WORKERS <- if (ON_CLUSTER) {
  # Keep ATE parallelism conservative on cluster to avoid memory spikes.
  max(1L, min(as.integer(N_CORES), as.integer(floor(N_CORES / 8))))
} else {
  N_CORES
}
ATE_WORKERS <- env_int("PP_ATE_WORKERS", ATE_WORKERS, 1L)
ATE_WORKERS <- min(as.integer(N_CORES), as.integer(ATE_WORKERS))
ATE_BATCH_SIZE <- max(ATE_WORKERS, 2L * ATE_WORKERS)
ATE_N_SIMS <- if (TEST) 1L else as.integer(N_SIMS)
ATE_N_TAU_SIMS <- if (TEST) 1L else as.integer(N_TAU_SIMS)
ATE_N_TAU_I <- if (TEST) 1L else as.integer(N_TAU_I)
ATE_MAXIT <- if (TEST) 300L else 1000L
ATE_N_SIMS <- env_int("PP_ATE_N_SIMS", ATE_N_SIMS, 1L)
ATE_N_TAU_SIMS <- env_int("PP_ATE_N_TAU_SIMS", ATE_N_TAU_SIMS, 1L)
ATE_N_TAU_I <- env_int("PP_ATE_N_TAU_I", ATE_N_TAU_I, 1L)
ATE_MAXIT <- env_int("PP_ATE_MAXIT", ATE_MAXIT, 1L)
ATE_COMPUTE_TAU <- tolower(Sys.getenv("PP_ATE_COMPUTE_TAU", "true")) %in% c("1", "true", "yes", "y")
SIM_FILTRATION_AWARE <- tolower(Sys.getenv("PP_SIM_FILTRATION_AWARE", "true")) %in% c("1", "true", "yes", "y")
ATE_CONTROL_FILTRATION_AWARE <- SIM_FILTRATION_AWARE &&
  (tolower(Sys.getenv("PP_ATE_CONTROL_FILTRATION_AWARE", "true")) %in% c("1", "true", "yes", "y"))
FILTER_HIGH_MU_FITS <- tolower(Sys.getenv("PP_FILTER_HIGH_MU_FITS", "true")) %in% c("1", "true", "yes", "y")
MU_FAIL_MULTIPLIER <- env_num("PP_MU_FAIL_MULTIPLIER", 2, min_value = 1)
RUN_SEM_PILOT <- tolower(Sys.getenv("PP_RUN_SEM_PILOT", "false")) %in% c("1", "true", "yes", "y")
if (RUN_SEM_PILOT && ON_CLUSTER) {
  log_msg("[SEM PILOT] Disabled on cluster (local tuning only).")
  RUN_SEM_PILOT <- FALSE
}
PILOT_ONLY <- tolower(Sys.getenv("PP_PILOT_ONLY", "false")) %in% c("1", "true", "yes", "y")
SEM_PILOT_ITERS <- suppressWarnings(as.integer(Sys.getenv("PP_SEM_PILOT_ITERS", "200")))
if (!is.finite(SEM_PILOT_ITERS) || is.na(SEM_PILOT_ITERS) || SEM_PILOT_ITERS < 10L) SEM_PILOT_ITERS <- 200L
SEM_PILOT_SIMS <- suppressWarnings(as.integer(Sys.getenv("PP_SEM_PILOT_SIMS", "1")))
if (!is.finite(SEM_PILOT_SIMS) || is.na(SEM_PILOT_SIMS) || SEM_PILOT_SIMS < 1L) SEM_PILOT_SIMS <- min(SIM_SIZE, 12L)
SEM_PILOT_CORES <- suppressWarnings(as.integer(Sys.getenv("PP_SEM_PILOT_CORES", as.character(max(1L, floor(0.8 * N_CORES))))))
if (!is.finite(SEM_PILOT_CORES) || is.na(SEM_PILOT_CORES) || SEM_PILOT_CORES < 1L) SEM_PILOT_CORES <- max(1L, floor(0.8 * N_CORES))
SEM_WORKERS_DEFAULT <- if (ON_CLUSTER && TEST) {
  min(8L, as.integer(N_CORES), as.integer(SIM_SIZE))
} else {
  min(as.integer(N_CORES), as.integer(SIM_SIZE))
}
SEM_WORKERS <- env_int("PP_SEM_WORKERS", SEM_WORKERS_DEFAULT, 1L)
SEM_WORKERS <- min(as.integer(N_CORES), as.integer(SEM_WORKERS))
log_msg("=== ", JOB_ID, " | ", MODE, " | ", N_CORES, " cores x ", SIM_SIZE, " sims ===")
log_msg("SEM workers=", SEM_WORKERS)
log_msg("SEM adaptive inner=", SEM_EM_ADAPTIVE_ITER, " | outer=", SEM_N_ITER, " | labellings=", SEM_N_LABELLINGS)
log_msg("SEM spec: n_props=", SEM_N_PROPS,
        " | param_cadence=", SEM_PARAM_UPDATE_CADENCE,
        " | proposal_cadence=", SEM_PROPOSAL_UPDATE_CADENCE,
        " | param_refit_cadence=", SEM_PARAM_REFIT_CADENCE,
        " | staleness_trigger_every=", SEM_STALENESS_TRIGGER_EVERY,
        " | change_factor=", SEM_CHANGE_FACTOR,
        " | include_starting=", SEM_INCLUDE_STARTING,
        " | update_starting=", SEM_UPDATE_STARTING,
        " | update_control_params=", SEM_UPDATE_CONTROL_PARAMS)
log_msg("ATE workers=", ATE_WORKERS)
log_msg("ATE batch size=", ATE_BATCH_SIZE)
log_msg("ATE config: n_sims=", ATE_N_SIMS,
        " | n_tau_sims=", ATE_N_TAU_SIMS,
        " | n_tau_i=", ATE_N_TAU_I,
        " | compute_tau=", ATE_COMPUTE_TAU,
        " | maxit=", ATE_MAXIT,
        " | sim_filtration_aware=", SIM_FILTRATION_AWARE,
        " | control_filtration_aware=", ATE_CONTROL_FILTRATION_AWARE)
log_msg("High-mu failure filter: enabled=", FILTER_HIGH_MU_FITS,
        " | mu_fail_multiplier=", MU_FAIL_MULTIPLIER)
log_msg("Save mode: ", ifelse(SAVE_LIGHT, "LIGHT (trimmed traces)", "FULL (all objects)"))
if (RUN_SEM_PILOT) {
  log_msg("SEM pilot: iters=", SEM_PILOT_ITERS,
          " | sims=", SEM_PILOT_SIMS,
          " | cores=", SEM_PILOT_CORES,
          " | pilot_only=", PILOT_ONLY)
}
log_msg("Base seed=", BASE_SEED)
log_msg("Post-treatment time multiplier=", POST_TIME_MULTIPLIER, " | END_TIME=", END_TIME)
log_msg("Scenario=", SCENARIO_ID,
        " | sim_kernel=", SIM_KERNEL,
        " | fit_kernel=", FIT_KERNEL,
        " | sim_spatial_kernel=", SIM_SPATIAL_KERNEL,
        " | fit_spatial_kernel=", FIT_SPATIAL_KERNEL,
        " | spatial_q=", HAWKES_SPATIAL_POWER_Q,
        " | treatment_assignment=", TREATMENT_ASSIGNMENT,
        " | control_K=", CONTROL_K,
        " | treated_K=", TREATED_K,
        " | base_mu=", BASE_MU,
        " | mu_scale=", MU_SCALE,
        " | target_points=", ifelse(is.finite(TARGET_POINTS), TARGET_POINTS, NA),
        " | expected_points_per_mu=", signif(EXPECTED_POINTS_PER_MU, 6),
        " | true_mu=", signif(TRUE_MU, 6))
log_msg(
  "Hawkes t_trunc: ",
  if (is.null(HAWKES_T_TRUNC)) {
    "disabled"
  } else {
    sprintf("%.4g (%s)", HAWKES_T_TRUNC, HAWKES_T_TRUNC_SOURCE)
  }
)
log_msg("Output (canonical): ", SAVE_DIR)

# ------------------------------------------------------------------
# Helper: create a parallel cluster with PPDisentangle loaded
# ------------------------------------------------------------------
make_cluster <- function(n_cores, rng_seed = NULL) {
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, {
    prepend_user_lib_paths <- function() {
      user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
      if (!nzchar(user_lib)) return(invisible(NULL))
      extra_libs <- strsplit(user_lib, .Platform$path.sep, fixed = TRUE)[[1]]
      extra_libs <- extra_libs[nzchar(extra_libs)]
      if (length(extra_libs) > 0L) {
        .libPaths(c(extra_libs, .libPaths()))
      }
    }
    wait_for_namespace <- function(pkg, timeout_s = 180L, sleep_s = 5L) {
      prepend_user_lib_paths()
      lock_dir <- file.path(Sys.getenv("R_LIBS_USER", unset = ""), paste0("00LOCK-", pkg))
      start_time <- Sys.time()
      repeat {
        prepend_user_lib_paths()
        if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
        waited_s <- as.integer(difftime(Sys.time(), start_time, units = "secs"))
        if (waited_s >= timeout_s) {
          stop(sprintf(
            "Worker could not load '%s' after %ds. .libPaths()=%s",
            pkg, timeout_s, paste(.libPaths(), collapse = " | ")
          ))
        }
        if (nzchar(lock_dir) && dir.exists(lock_dir)) {
          message(sprintf("Worker waiting for '%s' lock (%s); waited %ds...", pkg, lock_dir, waited_s))
        } else {
          message(sprintf("Worker waiting for '%s' visibility; waited %ds...", pkg, waited_s))
        }
        Sys.sleep(sleep_s)
      }
    }
    user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
    if (nzchar(user_lib)) {
      extra_libs <- strsplit(user_lib, .Platform$path.sep, fixed = TRUE)[[1]]
      extra_libs <- extra_libs[nzchar(extra_libs)]
      if (length(extra_libs) > 0L) {
        .libPaths(c(extra_libs, .libPaths()))
      }
    }
    wait_for_namespace("PPDisentangle")
    library(PPDisentangle)
    library(R.utils)
  })
  if (!is.null(rng_seed) && is.finite(rng_seed)) {
    parallel::clusterSetRNGStream(cl, iseed = as.integer(rng_seed))
  }
  return(cl)
}

run_maybe_parallel <- function(cl, jobs, fun, use_parallel) {
  if (use_parallel) {
    parLapply(cl = cl, X = jobs, fun = fun)
  } else {
    lapply(jobs, fun)
  }
}

format_sem_error_message <- function(err_list) {
  paste(
    vapply(err_list, function(x) sprintf("sim %d: %s", x$sim_id, x$message), character(1)),
    collapse = " | "
  )
}

export_globals <- function(cl) {
  clusterExport(cl, c(
    "OMEGA", "partition", "N_SIMS", "MAX_TIME",
    "N_TAU_I", "N_TAU_SIMS", "N_TAU_I_TRUE",
    "TREATMENT_TIME", "END_TIME",
    "TIME_INT", "state_spaces", "partition_processes",
    "treated_partitions", "treated_state_space",
    "control_state_space", "hawkes_par_1", "hawkes_par_2", "hawkes_fit_par_1", "hawkes_fit_par_2",
    "N_PROPOSALS", "SEM_EM_ADAPTIVE_ITER", "SEM_N_ITER", "SEM_N_LABELLINGS",
    "SEM_PARAM_UPDATE_CADENCE", "SEM_PROPOSAL_UPDATE_CADENCE",
    "SEM_PARAM_REFIT_CADENCE", "SEM_N_PROPS", "SEM_CHANGE_FACTOR", "SEM_STALENESS_TRIGGER_EVERY",
    "SEM_INCLUDE_STARTING", "SEM_UPDATE_STARTING", "SEM_UPDATE_CONTROL_PARAMS",
    "SIM_FILTRATION_AWARE", "FIT_KERNEL", "FIT_SPATIAL_KERNEL", "HAWKES_SPATIAL_POWER_Q",
    "HAWKES_T_TRUNC"
  ), envir = .GlobalEnv)
}

export_ate_globals <- function(cl) {
  clusterExport(cl, c(
    "OMEGA", "partition", "N_SIMS", "MAX_TIME",
    "N_TAU_SIMS", "N_TAU_I", "TREATMENT_TIME", "END_TIME",
    "treated_partitions", "params_are_crazy"
  ), envir = .GlobalEnv)
}

# ------------------------------------------------------------------
# 1. Define true Hawkes parameters
# ------------------------------------------------------------------
hawkes_par_1 <- make_hawkes_params_for_kernel(TRUE_MU, CONTROL_K, SIM_KERNEL, SIM_SPATIAL_KERNEL)
hawkes_par_2 <- make_hawkes_params_for_kernel(TRUE_MU, TREATED_K, SIM_KERNEL, SIM_SPATIAL_KERNEL)
hawkes_fit_par_1 <- make_hawkes_params_for_kernel(TRUE_MU, CONTROL_K, FIT_KERNEL, FIT_SPATIAL_KERNEL)
hawkes_fit_par_2 <- make_hawkes_params_for_kernel(TRUE_MU, TREATED_K, FIT_KERNEL, FIT_SPATIAL_KERNEL)

partition <- quadrats(X = OMEGA, nx = NX, ny = NY)

all_nothing_ATE <- (hawkes_par_2$mu * TIME_INT * (1 / (1 - hawkes_par_2$K)) -
                      hawkes_par_1$mu * TIME_INT * (1 / (1 - hawkes_par_1$K))) / partition$n

count_points_by_partition <- function(df, partition) {
  counts <- numeric(partition$n)
  if (is.null(df) || nrow(df) < 1L) return(counts)
  for (ii in seq_len(partition$n)) {
    wi <- as.owin(partition[ii])
    counts[[ii]] <- sum(spatstat.geom::inside.owin(df$x, df$y, wi), na.rm = TRUE)
  }
  counts
}

simulate_assignment_reference_counts <- function() {
  set.seed(stage_seed(0L, 0L, 84L))
  reference <- sim_hawkes(
    params = hawkes_par_1,
    windowT = c(0, TREATMENT_TIME), windowS = OMEGA,
    optimized = TRUE
  )
  count_points_by_partition(data.frame(x = reference$x, y = reference$y), partition)
}

simulate_assignment_reference_points <- function() {
  set.seed(stage_seed(0L, 0L, 84L))
  reference <- sim_hawkes(
    params = hawkes_par_1,
    windowT = c(0, TREATMENT_TIME), windowS = OMEGA,
    optimized = TRUE
  )
  data.frame(x = reference$x, y = reference$y)
}

# Pretreatment Voronoi geometry with random tile-level treatment, projected
# onto the analysis grid via cell centroids. Uses the first n_seeds events
# from the reference catalogue as Dirichlet generators (~100 Voronoi areas).
assign_voronoi_random <- function(n_treated, n_seeds = 100L) {
  n_cells <- partition$n
  n_seeds <- as.integer(n_seeds)
  if (!is.finite(n_seeds) || n_seeds < 2L) n_seeds <- 100L
  out <- rep("control", n_cells)
  ref_df <- simulate_assignment_reference_points()
  counts <- count_points_by_partition(ref_df, partition)

  # Fallback if the reference catalogue is empty.
  if (is.null(ref_df) || nrow(ref_df) < 1L) {
    set.seed(stage_seed(0L, 0L, 42L))
    out[sample(seq_len(n_cells), n_treated)] <- "treated"
    return(list(processes = out, counts = counts))
  }

  # First n_seeds pretreatment events (time order from the simulator).
  n_take <- min(n_seeds, nrow(ref_df))
  seed_df <- ref_df[seq_len(n_take), , drop = FALSE]
  # spatstat::ppp requires a true owin; OMEGA is the numeric box c(xmin,xmax,ymin,ymax).
  omega_win <- if (inherits(OMEGA, "owin")) OMEGA else spatstat.geom::as.owin(OMEGA)
  seeds <- spatstat.geom::ppp(
    x = seed_df$x,
    y = seed_df$y,
    window = omega_win,
    check = FALSE
  )
  seeds <- spatstat.geom::unique.ppp(seeds)
  if (seeds$n < 2L) {
    set.seed(stage_seed(0L, 0L, 42L))
    out[sample(seq_len(n_cells), n_treated)] <- "treated"
    return(list(processes = out, counts = counts))
  }

  vor <- spatstat.geom::dirichlet(seeds)
  n_tiles <- vor$n
  n_tile_treated <- as.integer(round(n_tiles * TREAT_PROP))
  n_tile_treated <- max(1L, min(n_tiles - 1L, n_tile_treated))
  set.seed(stage_seed(0L, 0L, 87L))
  treated_tiles <- sample.int(n_tiles, n_tile_treated)

  tile_windows <- lapply(seq_len(n_tiles), function(j) as.owin(vor[j]))
  for (i in seq_len(n_cells)) {
    wi <- as.owin(partition[i])
    cx <- mean(wi$xrange)
    cy <- mean(wi$yrange)
    for (j in treated_tiles) {
      if (isTRUE(spatstat.geom::inside.owin(cx, cy, tile_windows[[j]]))) {
        out[[i]] <- "treated"
        break
      }
    }
  }

  # Match the usual exact treated-cell count used by other assignment rules.
  treated_now <- which(out == "treated")
  if (length(treated_now) > n_treated) {
    set.seed(stage_seed(0L, 0L, 88L))
    drop <- sample(treated_now, length(treated_now) - n_treated)
    out[drop] <- "control"
  } else if (length(treated_now) < n_treated) {
    rem <- which(out == "control")
    set.seed(stage_seed(0L, 0L, 88L))
    add <- sample(rem, n_treated - length(treated_now))
    out[add] <- "treated"
  }

  list(processes = out, counts = counts)
}

# Rook adjacency on the regular cell grid via tile centroids.
partition_rook_neighbors <- function(partition) {
  n <- partition$n
  cents <- matrix(NA_real_, n, 2L)
  for (i in seq_len(n)) {
    wi <- as.owin(partition[i])
    cents[i, ] <- c(mean(wi$xrange), mean(wi$yrange))
  }
  xs <- sort(unique(round(cents[, 1L], 8L)))
  ys <- sort(unique(round(cents[, 2L], 8L)))
  dx <- if (length(xs) > 1L) min(diff(xs)) else Inf
  dy <- if (length(ys) > 1L) min(diff(ys)) else Inf
  tol <- 0.15 * min(dx, dy)
  lapply(seq_len(n), function(i) {
    ddx <- abs(cents[, 1L] - cents[i, 1L])
    ddy <- abs(cents[, 2L] - cents[i, 2L])
    which(
      (ddx < tol & abs(ddy - dy) < tol) |
        (ddy < tol & abs(ddx - dx) < tol)
    )
  })
}

# Grow a contiguous treated region from the hottest seed cell.
grow_contiguous_high_count <- function(counts, n_treated, neighbors) {
  n <- length(counts)
  set.seed(stage_seed(0L, 0L, 85L))
  seed <- which.max(counts + stats::runif(n) * 1e-9)
  selected <- seed
  while (length(selected) < n_treated) {
    frontier <- setdiff(unique(unlist(neighbors[selected], use.names = FALSE)), selected)
    if (length(frontier) < 1L) {
      rem <- setdiff(seq_len(n), selected)
      pick <- rem[which.max(counts[rem] + stats::runif(length(rem)) * 1e-9)]
    } else {
      pick <- frontier[which.max(counts[frontier] + stats::runif(length(frontier)) * 1e-9)]
    }
    selected <- c(selected, pick)
  }
  selected
}

# Logistic propensity on log(1 + pretreatment counts), calibrated to ~n_treated,
# then sample exactly n_treated cells without replacement.
sample_count_propensity <- function(counts, n_treated, slope = 2) {
  n <- length(counts)
  x <- log1p(pmax(counts, 0))
  sx <- stats::sd(x)
  if (!is.finite(sx) || sx < 1e-12) {
    set.seed(stage_seed(0L, 0L, 85L))
    return(sample(seq_len(n), n_treated))
  }
  x <- as.numeric(scale(x))
  target <- as.numeric(n_treated)
  f <- function(a) sum(stats::plogis(a + slope * x)) - target
  lo <- f(-40)
  hi <- f(40)
  a <- if (is.finite(lo) && is.finite(hi) && lo * hi < 0) {
    stats::uniroot(f, c(-40, 40))$root
  } else {
    0
  }
  pi <- pmax(stats::plogis(a + slope * x), 1e-12)
  set.seed(stage_seed(0L, 0L, 85L))
  sample(seq_len(n), n_treated, prob = pi)
}

make_partition_processes <- function(mode) {
  mode <- tolower(trimws(mode))
  n_treated <- as.integer(round(partition$n * TREAT_PROP))
  n_treated <- max(1L, min(partition$n - 1L, n_treated))
  out <- rep("control", partition$n)

  if (mode %in% c("random", "random_50pct", "random_50")) {
    set.seed(stage_seed(0L, 0L, 42L))
    out[sample(seq_len(partition$n), n_treated)] <- "treated"
    return(list(processes = out, counts = rep(NA_real_, partition$n)))
  }

  if (mode %in% c("lowest_count", "lowest_counts", "lowest_points", "lowest_count_50pct")) {
    counts <- simulate_assignment_reference_counts()
    set.seed(stage_seed(0L, 0L, 85L))
    tie_break <- stats::runif(length(counts))
    treated_cells <- order(counts, tie_break)[seq_len(n_treated)]
    out[treated_cells] <- "treated"
    return(list(processes = out, counts = counts))
  }

  if (mode %in% c("highest_count", "highest_counts", "highest_points", "highest_count_50pct")) {
    counts <- simulate_assignment_reference_counts()
    set.seed(stage_seed(0L, 0L, 85L))
    tie_break <- stats::runif(length(counts))
    treated_cells <- order(-counts, tie_break)[seq_len(n_treated)]
    out[treated_cells] <- "treated"
    return(list(processes = out, counts = counts))
  }

  if (mode %in% c(
    "count_propensity", "count_propensity_50pct", "propensity", "propensity_50pct"
  )) {
    counts <- simulate_assignment_reference_counts()
    treated_cells <- sample_count_propensity(counts, n_treated)
    out[treated_cells] <- "treated"
    return(list(processes = out, counts = counts))
  }

  if (mode %in% c(
    "contiguous_aoi", "contiguous_aoi_50pct",
    "contiguous_high_count", "contiguous_high_count_50pct"
  )) {
    counts <- simulate_assignment_reference_counts()
    treated_cells <- grow_contiguous_high_count(
      counts, n_treated, partition_rook_neighbors(partition)
    )
    out[treated_cells] <- "treated"
    return(list(processes = out, counts = counts))
  }

  if (mode %in% c(
    "voronoi_random", "voronoi_random_50pct",
    "voronoi", "pretreatment_voronoi_random"
  )) {
    return(assign_voronoi_random(n_treated))
  }

  stop(sprintf("Unknown PP_TREATMENT_ASSIGNMENT value: %s", mode))
}

partition_processes <- rep("control", partition$n)
assignment_info <- make_partition_processes(TREATMENT_ASSIGNMENT)
partition_processes <- assignment_info$processes
partition_assignment_counts <- assignment_info$counts
log_msg(
  "Treatment assignment: mode=", TREATMENT_ASSIGNMENT,
  " | treated_cells=", sum(partition_processes == "treated"),
  " | control_cells=", sum(partition_processes == "control")
)
if (any(is.finite(partition_assignment_counts))) {
  treated_counts <- partition_assignment_counts[partition_processes == "treated"]
  control_counts <- partition_assignment_counts[partition_processes == "control"]
  log_msg(
    "Assignment reference counts: treated median=", signif(stats::median(treated_counts), 4),
    " | control median=", signif(stats::median(control_counts), 4)
  )
}
treated_idx <- partition_processes == "treated"
control_state_space <- as.owin(partition[!treated_idx])
treated_state_space <- as.owin(partition[treated_idx])
treated_partitions  <- tilenames(partition[which(partition_processes == "treated")])
if (partition_processes[1] == "treated") {
  state_spaces <- list(treated_state_space, control_state_space)
} else {
  state_spaces <- list(control_state_space, treated_state_space)
}

log_msg("True all/nothing ATE:", round(all_nothing_ATE, 4))
pilot_only_mode <- RUN_SEM_PILOT && PILOT_ONLY

# Create cluster now that all globals are defined
cl <- NULL
if (!pilot_only_mode) {
  cl <- make_cluster(SEM_WORKERS, rng_seed = stage_seed(11L))
  export_globals(cl)
}

# ------------------------------------------------------------------
# 2. Compute true average one-flip ATE
# ------------------------------------------------------------------
if (!pilot_only_mode) {
  log_msg("Computing true tau_i ...")
  t0 <- proc.time()[3]
  tau_i_estim <- parSapply(cl = cl, X = 1:partition$n, FUN = function(i) {
    tau_i(i,
      partition = partition, treated_partitions = treated_partitions,
      statespace = OMEGA, windowT = c(TREATMENT_TIME, END_TIME),
      control_pp = hawkes_par_1, treated_pp = hawkes_par_2,
      n_sim = N_TAU_I_TRUE
    )
  })
  true_tau_1 <- mean(tau_i_estim)
  log_elapsed("True tau_i", proc.time()[3] - t0)
  log_msg("True one-flip ATE:", round(true_tau_1, 4))
} else {
  true_tau_1 <- NA_real_
  log_msg("[SEM PILOT] Skipping true tau_i in pilot-only mode.")
}

# ------------------------------------------------------------------
# 3. Generate observed data
# ------------------------------------------------------------------
log_msg("Generating", SIM_SIZE, "observed datasets ...")
t0 <- proc.time()[3]
obs_data <- lapply(1:SIM_SIZE, function(i) {
  set.seed(stage_seed(1L, i))
  pre_treat <- sim_hawkes(
    params = hawkes_par_1,
    windowT = c(0, TREATMENT_TIME), windowS = OMEGA,
    optimized = TRUE
  )
  pre_df <- data.frame(
    x = pre_treat$x, y = pre_treat$y, t = pre_treat$t,
    n = length(pre_treat$t),
    background = pre_treat$background,
    process = rep("control", length(pre_treat$t)),
    location_process = rep("control", length(pre_treat$t))
  )

  combined <- generate_inhomogeneous_hawkes(
    Omega = OMEGA, partition = partition,
    time_window = c(TREATMENT_TIME, END_TIME),
    partition_processes = partition_processes,
    state_spaces = state_spaces, filtration = pre_df,
    hawkes_params = list(control = hawkes_par_1, treated = hawkes_par_2)
  )
  post_df <- as.data.frame(combined)

  shared_cols <- intersect(names(pre_df), names(post_df))
  rbind(pre_df[, shared_cols], post_df[, shared_cols])
})
log_elapsed("Data generation", proc.time()[3] - t0, SIM_SIZE, SIM_SIZE)
log_msg("Points per sim:", paste(sapply(obs_data, nrow), collapse = ", "))
log_memory("post_data_gen")

# ------------------------------------------------------------------
# 4. Baseline labellings (oracle + naive)
# ------------------------------------------------------------------
if (!pilot_only_mode) {
  log_msg("Computing oracle and naive labellings ...")
  pp_labeled_oracle <- lapply(obs_data, function(s) {
    pre  <- as.data.frame(s) %>% filter(.data$t < TREATMENT_TIME) %>%
      mutate(inferred_process = "control", location_process = "control")
    post <- as.data.frame(s) %>% filter(.data$t >= TREATMENT_TIME)
    rbind(pre, oracle_labeling(post))
  })

  pp_labeled_naive <- lapply(obs_data, function(s) {
    pre  <- as.data.frame(s) %>% filter(.data$t < TREATMENT_TIME) %>%
      mutate(inferred_process = "control", location_process = "control")
    post <- as.data.frame(s) %>% filter(.data$t >= TREATMENT_TIME)
    rbind(pre, naive_labeling(post))
  })
}

# ------------------------------------------------------------------
# 5. Labelling proposals
# ------------------------------------------------------------------
if (!pilot_only_mode) {
  log_msg("Generating", N_PROPOSALS, "labelling proposals per sim ...")
  t0 <- proc.time()[3]
  proposal_jobs <- lapply(seq_along(obs_data), function(i) {
    list(i = i, seed = stage_seed(2L, i), data = obs_data[[i]])
  })
  gen_proposals <- function(job) {
    set.seed(job$seed)
    s <- job$data
    df <- as.data.frame(s)
    pre  <- df[df$t < TREATMENT_TIME, , drop = FALSE]
    post <- df[df$t >= TREATMENT_TIME, , drop = FALSE]
    pre$location_process <- "control"
    pre$inferred_process <- NULL
    Filter(Negate(is.null), lapply(1:N_PROPOSALS, function(i) {
      tryCatch({
        tmp <- simulation_labeling_hawkes_hawkes_fast(
          post, partition = partition,
          partition_process = partition_processes,
          statespace = OMEGA, state_spaces = state_spaces,
          windowT = c(TREATMENT_TIME, END_TIME),
          hawkes_params_control = hawkes_fit_par_1,
          hawkes_params_treated = NULL,
          change_factor = 1, filtration = pre, proximity_weight = 0,
          kernel = FIT_KERNEL,
          spatial_kernel = FIT_SPATIAL_KERNEL,
          spatial_q = HAWKES_SPATIAL_POWER_Q
        )
        pre$inferred_process <- "control"
        rbind(pre, tmp)
      }, error = function(e) NULL)
    }))
  }
  if (SEM_WORKERS > 1) {
    labelling_proposals <- run_maybe_parallel(cl, proposal_jobs, gen_proposals, TRUE)
  } else {
    labelling_proposals <- run_maybe_parallel(cl, proposal_jobs, gen_proposals, FALSE)
  }
  log_elapsed("Labelling proposals", proc.time()[3] - t0, SIM_SIZE, SIM_SIZE)

  pp_labeled_best_proposal <- lapply(seq_along(labelling_proposals), function(i) {
    props <- labelling_proposals[[i]]
    if (length(props) == 0) return(pp_labeled_naive[[i]])
    accs <- sapply(props, function(y) {
      keep <- which(y$t > TREATMENT_TIME)
      if (length(keep) < 2) return(0)
      mean(y$inferred_process[keep] == y$process[keep])
    })
    props[[which.max(accs)]]
  })
}

# ------------------------------------------------------------------
# 7. Adaptive SEM
# ------------------------------------------------------------------
log_msg("Running adaptive SEM ...")
t0 <- proc.time()[3]
sem_jobs <- lapply(seq_along(obs_data), function(i) {
  list(i = i, seed = stage_seed(3L, i), data = obs_data[[i]])
})
run_sem_core <- function(job, tuning = NULL, sem_inner_iter_override = NULL) {
  set.seed(job$seed)
  dat <- job$data
  # Start treated SEM parameters at true control Hawkes parameters.
  params_init  <- as.list(hawkes_fit_par_1)
  local_tuning <- if (is.null(tuning)) {
    list(
      param_update_cadence = SEM_PARAM_UPDATE_CADENCE,
      proposal_update_cadence = SEM_PROPOSAL_UPDATE_CADENCE,
      state_spaces = NULL,
      update_control_params = SEM_UPDATE_CONTROL_PARAMS,
      param_refit_cadence = SEM_PARAM_REFIT_CADENCE,
      iter = if (is.null(sem_inner_iter_override)) SEM_EM_ADAPTIVE_ITER else sem_inner_iter_override,
      n_props = SEM_N_PROPS,
      change_factor = SEM_CHANGE_FACTOR,
      stagnation_trigger_every = SEM_STALENESS_TRIGGER_EVERY,
      include_starting_data = SEM_INCLUDE_STARTING,
      update_starting_data = SEM_UPDATE_STARTING,
      verbose = FALSE
    )
  } else {
    tuning
  }

  adaptive_SEM(
    pp_data = dat, partition = partition,
    partition_processes = partition_processes,
    statespace = OMEGA,
    time_window = c(TREATMENT_TIME, END_TIME),
    treatment_time = TREATMENT_TIME,
    hawkes_params_control = hawkes_fit_par_1,
    hawkes_params_treated = params_init,
    N_labellings = SEM_N_LABELLINGS,
    N_iter = SEM_N_ITER, verbose = FALSE,
    hawkes_use_filtration_history = SIM_FILTRATION_AWARE,
    kernel = FIT_KERNEL,
    spatial_kernel = FIT_SPATIAL_KERNEL,
    spatial_q = HAWKES_SPATIAL_POWER_Q,
    t_trunc = HAWKES_T_TRUNC,
    adaptive_control = local_tuning
  )
}
run_sem <- function(job) {
  tryCatch(
    run_sem_core(job),
    error = function(e) {
      structure(
        list(sim_id = job$i, message = conditionMessage(e)),
        class = "ppdis_sem_error"
      )
    }
  )
}
is_sem_error <- function(x) inherits(x, "ppdis_sem_error")

sem_pilot_summary <- NULL
if (RUN_SEM_PILOT) {
  budget_constant <- as.integer(SEM_PILOT_ITERS * 10L) # keep n_props * iter fixed
  mk_cfg <- function(name, change_factor, n_props, include_starting, staleness = 10L,
                     proposal_cadence = 1L, param_cadence = 10L, update_starting = TRUE) {
    n_props <- as.integer(n_props)
    iters <- max(20L, as.integer(floor(budget_constant / max(1L, n_props))))
    list(
      name = name,
      change_factor = change_factor,
      n_props = n_props,
      iter = iters,
      include_starting = include_starting,
      staleness = as.integer(staleness),
      proposal_cadence = as.integer(proposal_cadence),
      param_cadence = as.integer(param_cadence),
      update_starting = isTRUE(update_starting)
    )
  }
  pilot_cfgs <- list(
    mk_cfg("cf005_np5_start0",  0.005,  5L, FALSE),
    mk_cfg("cf005_np10_start0", 0.005, 10L, FALSE),
    mk_cfg("cf005_np20_start0", 0.005, 20L, FALSE),
    mk_cfg("cf010_np5_start0",  0.010,  5L, FALSE),
    mk_cfg("cf010_np10_start0", 0.010, 10L, FALSE),
    mk_cfg("cf010_np20_start0", 0.010, 20L, FALSE),
    mk_cfg("cf020_np5_start0",  0.020,  5L, FALSE),
    mk_cfg("cf020_np10_start0", 0.020, 10L, FALSE),
    mk_cfg("cf020_np20_start0", 0.020, 20L, FALSE),
    mk_cfg("cf005_np10_start1", 0.005, 10L, TRUE),
    mk_cfg("cf010_np10_start1", 0.010, 10L, TRUE),
    mk_cfg("cf010_np10_start0_st5",  0.010, 10L, FALSE, staleness = 5L),
    mk_cfg("cf010_np10_start0_st20", 0.010, 10L, FALSE, staleness = 20L)
  )
  pilot_ids <- seq_len(min(SEM_PILOT_SIMS, length(obs_data)))
  pilot_jobs <- list()
  for (cfg_id in seq_along(pilot_cfgs)) {
    for (sid in pilot_ids) {
      pilot_jobs[[length(pilot_jobs) + 1L]] <- list(
        cfg_id = cfg_id,
        cfg = pilot_cfgs[[cfg_id]],
        job = list(i = sid, seed = stage_seed(31L, sid, cfg_id), data = obs_data[[sid]])
      )
    }
  }
  run_sem_pilot <- function(x) {
    cfg <- x$cfg
    tuning <- list(
      param_update_cadence = cfg$param_cadence,
      proposal_update_cadence = cfg$proposal_cadence,
      state_spaces = NULL,
      update_control_params = FALSE, # no control-parameter adaptation in this pilot
      iter = cfg$iter,
      n_props = cfg$n_props,
      change_factor = cfg$change_factor,
      stagnation_trigger_every = cfg$staleness,
      include_starting_data = cfg$include_starting,
      update_starting_data = cfg$update_starting,
      verbose = FALSE
    )
    out <- tryCatch(run_sem_core(x$job, tuning = tuning), error = function(e) NULL)
    if (is.null(out) || is.null(out$adaptive)) {
      return(data.frame(config = cfg$name, sim_id = x$job$i, ok = FALSE,
                        n_props = cfg$n_props, iter = cfg$iter,
                        end_metric = NA_real_, metric_gain = NA_real_,
                        mean_flips = NA_real_, stringsAsFactors = FALSE))
    }
    ad <- out$adaptive
    metric_vec <- if (!is.null(ad$metrics)) as.numeric(ad$metrics) else numeric(0)
    metric_vec <- metric_vec[is.finite(metric_vec)]
    m0 <- if (length(metric_vec) > 0) metric_vec[[1]] else NA_real_
    mT <- if (length(metric_vec) > 0) metric_vec[[length(metric_vec)]] else NA_real_
    data.frame(
      config = cfg$name,
      sim_id = x$job$i,
      ok = TRUE,
      n_props = cfg$n_props,
      iter = cfg$iter,
      end_metric = mT,
      metric_gain = mT - m0,
      mean_flips = if (!is.null(ad$average_flips) && length(ad$average_flips) > 0) mean(ad$average_flips, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  log_msg("Running SEM pilot grid ...")
  pilot_results <- NULL
  if (SEM_PILOT_CORES > 1L) {
    cl_pilot <- make_cluster(SEM_PILOT_CORES, rng_seed = stage_seed(12L))
    export_globals(cl_pilot)
    clusterExport(cl_pilot, c("run_sem_core", "run_sem_pilot", "stage_seed", "SEM_PILOT_ITERS"), envir = environment())
    pilot_results <- parLapply(cl_pilot, pilot_jobs, run_sem_pilot)
    stopCluster(cl_pilot)
  } else {
    pilot_results <- lapply(pilot_jobs, run_sem_pilot)
  }
  pilot_df <- do.call(rbind, pilot_results)
  sem_pilot_summary <- pilot_df %>%
    group_by(.data$config) %>%
    summarize(
      n = n(),
      n_ok = sum(.data$ok, na.rm = TRUE),
      n_props = mean(.data$n_props, na.rm = TRUE),
      iter = mean(.data$iter, na.rm = TRUE),
      mean_end_metric = mean(.data$end_metric, na.rm = TRUE),
      mean_metric_gain = mean(.data$metric_gain, na.rm = TRUE),
      mean_flips = mean(.data$mean_flips, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(.data$mean_metric_gain), desc(.data$mean_end_metric), .data$mean_flips)
  log_msg("SEM pilot summary:")
  print(sem_pilot_summary)
  if (PILOT_ONLY) {
    log_msg("[SEM PILOT] pilot_only=TRUE; skipping full SEM/ATE/plot/save.")
    if (exists("cl") && !is.null(cl)) {
      try(stopCluster(cl), silent = TRUE)
    }
    close(log_con)
    quit(save = "no", status = 0)
  }
}

if (SEM_WORKERS > 1) {
  # PSOCK workers need these helpers explicitly exported because `run_sem`
  # resolves `run_sem_core` in the worker global environment.
  clusterExport(cl, c("run_sem_core", "run_sem", "is_sem_error"), envir = .GlobalEnv)
  EM_results <- run_maybe_parallel(cl, sem_jobs, run_sem, TRUE)
  failed <- which(vapply(EM_results, is_sem_error, logical(1)))
  if (length(failed) > 0L) {
    fail_msg <- format_sem_error_message(EM_results[failed])
    log_msg("[SEM PARALLEL ERROR] ", fail_msg)
    log_msg("Retrying failed SEM jobs sequentially: ", paste(failed, collapse = ","))
    retry <- run_maybe_parallel(cl, sem_jobs[failed], run_sem, FALSE)
    retry_failed <- which(vapply(retry, is_sem_error, logical(1)))
    if (length(retry_failed) > 0L) {
      retry_msg <- format_sem_error_message(retry[retry_failed])
      stop("Adaptive SEM failed after sequential retry: ", retry_msg)
    }
    EM_results[failed] <- retry
  }
} else {
  EM_results <- run_maybe_parallel(cl, sem_jobs, run_sem, FALSE)
  failed <- which(vapply(EM_results, is_sem_error, logical(1)))
  if (length(failed) > 0L) {
    fail_msg <- format_sem_error_message(EM_results[failed])
    stop("Adaptive SEM failed: ", fail_msg)
  }
}
log_elapsed("Adaptive SEM", proc.time()[3] - t0, SIM_SIZE, SIM_SIZE)
log_memory("post_SEM")

# Free heavy intermediates that are no longer needed before ATE.
# Keep pp_labeled_best_proposal because it is still used in downstream summaries.
rm(labelling_proposals)
gc(verbose = FALSE)

# ------------------------------------------------------------------
# 7b. Extract adaptive SEM diagnostics (flips + accuracy per iteration)
# ------------------------------------------------------------------
sem_diagnostics <- lapply(seq_along(EM_results), function(i) {
  a <- EM_results[[i]]$adaptive
  if (is.null(a)) return(NULL)
  n_iter <- length(a$accuracies)
  if (n_iter < 1L) return(NULL)
  metric_vec <- rep(NA_real_, n_iter)
  if (!is.null(a$metrics) && length(a$metrics) > 0) {
    use_n <- min(n_iter, length(a$metrics))
    metric_vec[seq_len(use_n)] <- as.numeric(a$metrics[seq_len(use_n)])
  }
  data.frame(
    sim_id = i,
    iteration = seq_len(n_iter),
    accuracy = a$accuracies,
    metric = metric_vec,
    average_flips = a$average_flips,
    max_metric_flips = a$max_metric_flips,
    stringsAsFactors = FALSE
  )
})
sem_diagnostics_df <- do.call(rbind, Filter(Negate(is.null), sem_diagnostics))
if (!is.null(sem_diagnostics_df) && nrow(sem_diagnostics_df) > 0) {
  log_msg("SEM diagnostics: ", nrow(sem_diagnostics_df), " rows across ", length(unique(sem_diagnostics_df$sim_id)), " sims")
}

# ------------------------------------------------------------------
# 8. Assemble all labellings
# ------------------------------------------------------------------
labellings <- list(
  oracle   = pp_labeled_oracle,
  naive    = pp_labeled_naive,
  best     = pp_labeled_best_proposal,
  SEM_adaptive = lapply(EM_results, function(x) x$adaptive$adaptive_labelling)
)

# ------------------------------------------------------------------
# 9. Classification accuracy summary
# ------------------------------------------------------------------
label_recovery_metrics <- function(labellings_obj, kept_idx = NULL) {
  lapply(names(labellings_obj), function(nm) {
    idx <- if (is.null(kept_idx)) seq_along(labellings_obj[[nm]]) else kept_idx(nm)
    rows <- lapply(idx, function(i) {
      y <- labellings_obj[[nm]][[i]]
      keep <- which(y$t > TREATMENT_TIME)
      if (length(keep) < 1L) return(NULL)
      truth <- factor(y$process[keep], levels = c("control", "treated"))
      pred <- factor(y$inferred_process[keep], levels = c("control", "treated"))
      tab <- table(truth, pred)
      tn <- unname(tab["control", "control"])
      fp <- unname(tab["control", "treated"])
      fn <- unname(tab["treated", "control"])
      tp <- unname(tab["treated", "treated"])
      denom <- tn + fp + fn + tp
      precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
      recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
      specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
      f1 <- if (is.finite(precision + recall) && (precision + recall) > 0) {
        2 * precision * recall / (precision + recall)
      } else {
        NA_real_
      }
      data.frame(
        method = nm,
        sim_id = i,
        n_post = denom,
        accuracy = if (denom > 0) (tp + tn) / denom else NA_real_,
        balanced_accuracy = mean(c(recall, specificity), na.rm = TRUE),
        precision_treated = precision,
        recall_treated = recall,
        specificity_control = specificity,
        f1_treated = f1,
        true_control_as_control = tn,
        true_control_as_treated = fp,
        true_treated_as_control = fn,
        true_treated_as_treated = tp,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, Filter(Negate(is.null), rows))
  }) %>%
    Filter(Negate(is.null), .) %>%
    do.call(rbind, .)
}

normalize_reported_labelling <- function(x) {
  x <- as.character(x)
  x[x %in% c("SEM_full", "SEM_adaptive")] <- "SEM"
  x
}

summarize_label_recovery <- function(label_recovery_df) {
  if (is.null(label_recovery_df) || nrow(label_recovery_df) < 1L) return(NULL)
  label_recovery_df %>%
    mutate(method = normalize_reported_labelling(.data$method)) %>%
    group_by(.data$method) %>%
    summarize(
      n_kept = n(),
      mean_accuracy = mean(.data$accuracy, na.rm = TRUE),
      mean_balanced_accuracy = mean(.data$balanced_accuracy, na.rm = TRUE),
      mean_precision_treated = mean(.data$precision_treated, na.rm = TRUE),
      mean_recall_treated = mean(.data$recall_treated, na.rm = TRUE),
      mean_specificity_control = mean(.data$specificity_control, na.rm = TRUE),
      mean_f1_treated = mean(.data$f1_treated, na.rm = TRUE),
      .groups = "drop"
    )
}

log_msg("")
log_msg("=== Classification Accuracy (pre-fit filter; final reported after high-mu filter) ===")
label_recovery_all_prefilter <- label_recovery_metrics(labellings)
class_metrics <- summarize_label_recovery(label_recovery_all_prefilter)
print(class_metrics)

# ------------------------------------------------------------------
# 10. ATE estimation
# ------------------------------------------------------------------
log_memory("pre_ATE")
log_msg("Estimating ATEs ...")
t0 <- proc.time()[3]

slim_for_ate <- function(df) {
  keep <- intersect(c("x", "y", "t", "inferred_process", "location_process"), names(df))
  out <- as.data.frame(df[, keep, drop = FALSE])
  if (!"inferred_process" %in% names(out) && "location_process" %in% names(out)) {
    out$inferred_process <- out$location_process
  }
  out
}

estimate_control_params_from_labeling <- function(post_data, filt_data) {
  fit_hawkes_with_filtration_fn <- getFromNamespace("fit_hawkes_with_filtration", "PPDisentangle")
  fit_hawkes_fn <- getFromNamespace("fit_hawkes", "PPDisentangle")
  ctrl_realiz <- post_data[post_data$inferred_process == "control", , drop = FALSE]
  dt_fit <- END_TIME - TREATMENT_TIME
  if (!is.finite(dt_fit) || dt_fit <= 0) dt_fit <- 1
  if (nrow(ctrl_realiz) < 2) {
    return(make_hawkes_params_for_kernel(max(1e-8, nrow(ctrl_realiz) / dt_fit), 0.01, FIT_KERNEL, FIT_SPATIAL_KERNEL))
  }
  ctrl_init <- make_hawkes_params_for_kernel(max(1e-8, nrow(ctrl_realiz) / dt_fit), 0.01, FIT_KERNEL, FIT_SPATIAL_KERNEL)
  if (isTRUE(ATE_CONTROL_FILTRATION_AWARE)) {
    fit_ctrl <- fit_hawkes_with_filtration_fn(
      params_init = ctrl_init,
      realiz = ctrl_realiz,
      filtration = filt_data,
      windowT = c(TREATMENT_TIME, END_TIME),
      windowS = OMEGA,
      maxit = ATE_MAXIT,
      poisson_flag = FALSE,
      zero_background_region = treated_state_space,
      kernel = FIT_KERNEL,
      spatial_kernel = FIT_SPATIAL_KERNEL,
      spatial_q = HAWKES_SPATIAL_POWER_Q,
      t_trunc = HAWKES_T_TRUNC
    )$par
    empirical_rate <- if (dt_fit > 0) nrow(ctrl_realiz) / dt_fit else Inf
    fitted_rate <- as.numeric(fit_ctrl$mu) / max(1e-6, 1 - as.numeric(fit_ctrl$K))
    degenerate_fit <- (!is.finite(fitted_rate) || !is.finite(empirical_rate) ||
                       (as.numeric(fit_ctrl$K) >= 0.98) ||
                       (fitted_rate < 0.2 * empirical_rate))
    if (isTRUE(degenerate_fit)) {
      legacy <- fit_hawkes_fn(
        ctrl_init,
        realiz = ctrl_realiz,
        zero_background_region = treated_state_space,
        windowT = c(TREATMENT_TIME, END_TIME),
        windowS = OMEGA,
        trace = 0,
        maxit = ATE_MAXIT,
        density_approx = FALSE,
        numeric_integral = FALSE,
        poisson_flag = FALSE,
        t_trunc = HAWKES_T_TRUNC,
        kernel = FIT_KERNEL,
        spatial_kernel = FIT_SPATIAL_KERNEL,
        spatial_q = HAWKES_SPATIAL_POWER_Q
      )$par
      fit_ctrl <- as_hawkes_params(as.list(legacy), FIT_KERNEL, FIT_SPATIAL_KERNEL, HAWKES_SPATIAL_POWER_Q)
    }
  } else {
    fit_ctrl <- fit_hawkes_fn(
      ctrl_init,
      realiz = ctrl_realiz,
      zero_background_region = treated_state_space,
      windowT = c(TREATMENT_TIME, END_TIME),
      windowS = OMEGA,
      trace = 0,
      maxit = ATE_MAXIT,
      density_approx = FALSE,
      numeric_integral = FALSE,
      poisson_flag = FALSE,
      t_trunc = HAWKES_T_TRUNC,
      kernel = FIT_KERNEL,
      spatial_kernel = FIT_SPATIAL_KERNEL,
      spatial_q = HAWKES_SPATIAL_POWER_Q
    )$par
    fit_ctrl <- as_hawkes_params(as.list(fit_ctrl), FIT_KERNEL, FIT_SPATIAL_KERNEL, HAWKES_SPATIAL_POWER_Q)
  }
  fit_ctrl
}

tasks <- list()
pre_histories <- lapply(obs_data, function(df) {
  pre <- as.data.frame(df) %>% filter(.data$t < TREATMENT_TIME)
  if (nrow(pre) > 0) {
    pre$location_process <- "control"
    pre$inferred_process <- "control"
  }
  slim_for_ate(pre)
})
for (nm in names(labellings)) {
  for (i in seq_along(labellings[[nm]])) {
    post_x <- labellings[[nm]][[i]] %>% filter(.data$t >= TREATMENT_TIME)
    filt_x <- if (SIM_FILTRATION_AWARE) pre_histories[[i]] else pre_histories[[i]][0, , drop = FALSE]
    tasks[[length(tasks) + 1]] <- list(
      x = slim_for_ate(post_x),
      filtration_data = filt_x,
      labelling_name = nm, hawkes_params = NULL,
      seed = stage_seed(4L, i, length(tasks) + 1L)
    )
  }
}
for (i in seq_along(EM_results)) {
  post_tmp <- EM_results[[i]]$adaptive$adaptive_labelling %>% filter(.data$t >= TREATMENT_TIME)
  filt_x <- if (SIM_FILTRATION_AWARE) pre_histories[[i]] else pre_histories[[i]][0, , drop = FALSE]
  post_slim <- slim_for_ate(post_tmp)
  ctrl_sem_full <- tryCatch(
    estimate_control_params_from_labeling(post_slim, filt_x),
    error = function(e) {
      post_fallback <- post_slim
      post_fallback$inferred_process <- post_fallback$location_process
      estimate_control_params_from_labeling(post_fallback, filt_x)
    }
  )
  tasks[[length(tasks) + 1]] <- list(
    x = post_slim,
    filtration_data = filt_x,
    labelling_name = "SEM_full",
    seed = stage_seed(4L, i, length(tasks) + 1L),
    hawkes_params = list(
      control = ctrl_sem_full,
      treated = EM_results[[i]]$hawkes_params_treated
    )
  )
}

SKIP_CRAZY <- ON_CLUSTER
log_msg(sprintf("[CRAZY PARAMS] skip mode: %s (mode default)", ifelse(SKIP_CRAZY, "ON", "OFF")))
crazy_idx <- integer(0)
for (k in seq_along(tasks)) {
  tsk <- tasks[[k]]
  if (tsk$labelling_name == "SEM_full" && !is.null(tsk$hawkes_params)) {
    if (params_are_crazy(tsk$hawkes_params$control, tsk$hawkes_params$treated)) {
      crazy_idx <- c(crazy_idx, k)
      ctrl <- tsk$hawkes_params$control
      treat <- tsk$hawkes_params$treated
      log_msg(sprintf("[CRAZY PARAMS] task %d (SEM_full sim %d): control mu=%.2g K=%.3f | treated mu=%.2g K=%.3f",
        k, ((k - 1) %% SIM_SIZE) + 1,
        ctrl$mu, ctrl$K, treat$mu, treat$K))
    }
  }
}
if (length(crazy_idx) > 0) {
  log_msg(sprintf("[CRAZY PARAMS] %d tasks with explosive params (K>0.98 or mu>1e5)", length(crazy_idx)))
  log_msg("[CRAZY PARAMS] Retaining all tasks for inspection; explosive ones are excluded from ATE execution queue")
}
skipped_explosive_tasks <- if (length(crazy_idx) > 0) {
  do.call(rbind, lapply(crazy_idx, function(k) {
    tsk <- tasks[[k]]
    ctrl <- tsk$hawkes_params$control
    treat <- tsk$hawkes_params$treated
    data.frame(
      task_idx = k,
      sim_id = ((k - 1) %% SIM_SIZE) + 1,
      labelling = tsk$labelling_name,
      control_mu = ctrl$mu, control_K = ctrl$K,
      treated_mu = treat$mu, treated_K = treat$K,
      stringsAsFactors = FALSE
    )
  }))
} else {
  data.frame(
    task_idx = integer(0), sim_id = integer(0), labelling = character(0),
    control_mu = numeric(0), control_K = numeric(0),
    treated_mu = numeric(0), treated_K = numeric(0),
    stringsAsFactors = FALSE
  )
}
ate_run_idx <- setdiff(seq_along(tasks), crazy_idx)

ATE_env <- new.env(parent = baseenv())
ATE_env$ATE_estim_hawkes <- getFromNamespace("ATE_estim_hawkes", "PPDisentangle")
ATE_env$withTimeout <- getFromNamespace("withTimeout", "R.utils")
ATE_env$params_are_crazy <- params_are_crazy
ATE_env$OMEGA <- OMEGA
ATE_env$partition <- partition
ATE_env$treated_partitions <- treated_partitions
ATE_env$N_SIMS <- N_SIMS
ATE_env$N_TAU_SIMS <- N_TAU_SIMS
ATE_env$N_TAU_I <- N_TAU_I
ATE_env$ATE_N_SIMS <- ATE_N_SIMS
ATE_env$ATE_N_TAU_SIMS <- ATE_N_TAU_SIMS
ATE_env$ATE_N_TAU_I <- ATE_N_TAU_I
ATE_env$ATE_MAXIT <- ATE_MAXIT
ATE_env$ATE_COMPUTE_TAU <- ATE_COMPUTE_TAU
ATE_env$ATE_CONTROL_FILTRATION_AWARE <- ATE_CONTROL_FILTRATION_AWARE
ATE_env$TRUE_CONTROL_HAWKES_INIT <- hawkes_fit_par_1
ATE_env$FIT_KERNEL <- FIT_KERNEL
ATE_env$FIT_SPATIAL_KERNEL <- FIT_SPATIAL_KERNEL
ATE_env$HAWKES_SPATIAL_POWER_Q <- HAWKES_SPATIAL_POWER_Q
ATE_env$HAWKES_T_TRUNC <- HAWKES_T_TRUNC
ATE_env$TREATMENT_TIME <- TREATMENT_TIME
ATE_env$END_TIME <- END_TIME
ATE_env$MAX_TIME <- MAX_TIME

task_function <- function(task) {
  if (!is.null(task$seed) && is.finite(task$seed)) {
    set.seed(as.integer(task$seed))
  }
  if (!is.null(task$hawkes_params) && params_are_crazy(task$hawkes_params$control, task$hawkes_params$treated)) {
    return(NULL)
  }
  r <- tryCatch(
    withTimeout(
      ATE_estim_hawkes(
        statespace = OMEGA, partition = partition,
        observed_data = task$x,
        filtration_data = task$filtration_data,
        control_filtration_aware = ATE_CONTROL_FILTRATION_AWARE,
        treated_params_init = TRUE_CONTROL_HAWKES_INIT,
        treated_partitions = treated_partitions,
        hawkes_params = task$hawkes_params,
        n_sims = ATE_N_SIMS, n_tau_sims = ATE_N_TAU_SIMS, n_tau_i = ATE_N_TAU_I,
        windowT = c(TREATMENT_TIME, END_TIME), windowS = OMEGA,
        maxit = ATE_MAXIT,
        explosive_K_threshold = 0.98,
        poisson_flags = list(control = FALSE, treated = FALSE),
        keep_all_nothing_sim = FALSE,
        compute_tau = ATE_COMPUTE_TAU,
        kernel = FIT_KERNEL,
        spatial_kernel = FIT_SPATIAL_KERNEL,
        spatial_q = HAWKES_SPATIAL_POWER_Q,
        t_trunc = HAWKES_T_TRUNC
      ),
      timeout = MAX_TIME, onTimeout = "error"
    ),
    error = function(e) NULL
  )
  if (!is.null(r) && isTRUE(r$skipped_explosive)) return(NULL)
  r
}
environment(task_function) <- ATE_env

if (N_CORES > 1) gc(verbose = FALSE)
ATE_SEQUENTIAL <- TEST
results_flat <- vector("list", length(tasks))
if (ATE_SEQUENTIAL) {
  log_msg("ATE estimation: sequential (TEST mode)")
  stopCluster(cl)
  if (length(ate_run_idx) > 0L) {
    results_flat[ate_run_idx] <- lapply(tasks[ate_run_idx], task_function)
  }
} else {
  # Recreate a fresh, smaller cluster for ATE to avoid stale worker state and
  # reduce memory pressure from too many concurrent heavy ATE jobs.
  stopCluster(cl)
  if (length(ate_run_idx) > 0L) {
    gc(verbose = FALSE)
    cl_ate <- make_cluster(min(ATE_WORKERS, length(ate_run_idx)), rng_seed = stage_seed(13L))
    export_ate_globals(cl_ate)
    task_ids <- ate_run_idx
    task_batches <- split(task_ids, ceiling(seq_along(task_ids) / ATE_BATCH_SIZE))
    batch_failed <- FALSE
    done_total <- 0L
    ok_total <- 0L
    null_total <- 0L
    for (b in seq_along(task_batches)) {
      idx <- task_batches[[b]]
      log_msg(sprintf("[ATE BATCH %d/%d] start (tasks %d..%d, n=%d)",
                      b, length(task_batches), min(idx), max(idx), length(idx)))
      tb0 <- proc.time()[3]
      batch_res <- tryCatch(
        parLapply(cl_ate, tasks[idx], fun = task_function),
        error = function(e) {
          log_msg("[ATE PARALLEL ERROR] batch ", b, "/", length(task_batches), ": ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(batch_res)) {
        batch_failed <- TRUE
        break
      }
      results_flat[idx] <- batch_res
      n_ok <- sum(vapply(batch_res, function(z) !is.null(z), logical(1)))
      n_null <- length(batch_res) - n_ok
      done_total <- done_total + length(batch_res)
      ok_total <- ok_total + n_ok
      null_total <- null_total + n_null
      batch_elapsed <- proc.time()[3] - tb0
      log_msg(sprintf("[ATE BATCH %d/%d] done in %.1f s | ok=%d null=%d | cumulative %d/%d (ok=%d null=%d)",
                      b, length(task_batches), batch_elapsed, n_ok, n_null,
                      done_total, length(ate_run_idx), ok_total, null_total))
      rm(batch_res); gc(verbose = FALSE)
      log_elapsed("ATE parallel batch", proc.time()[3] - t0, done_total, length(ate_run_idx))
    }
    stopCluster(cl_ate)
    if (batch_failed) {
      log_msg("[ATE PARALLEL ERROR] Falling back to sequential for robustness.")
      results_flat[ate_run_idx] <- lapply(tasks[ate_run_idx], task_function)
    }
  } else {
    log_msg("[ATE] No runnable tasks after explosive-parameter screening.")
  }
}
log_elapsed("ATE estimation", proc.time()[3] - t0, length(ate_run_idx), length(tasks))
log_memory("post_ATE")

# ------------------------------------------------------------------
# 11. Collect and save results
# ------------------------------------------------------------------
# Screen fitted Hawkes mu values and mark failed fits when they exceed
# a configurable multiple of the true generating mu values.
extract_mu_pair <- function(r, task_obj = NULL) {
  safe_mu <- function(pp) {
    if (is.null(pp)) return(NA_real_)
    mu <- suppressWarnings(as.numeric(pp$mu))
    if (length(mu) < 1L) return(NA_real_)
    mu[[1]]
  }
  ctrl_pp <- NULL
  treat_pp <- NULL
  if (!is.null(r)) {
    ctrl_pp <- r$control_pp
    treat_pp <- r$treated_pp
  }
  if (is.null(ctrl_pp) && !is.null(task_obj) && !is.null(task_obj$hawkes_params)) {
    ctrl_pp <- task_obj$hawkes_params$control
  }
  if (is.null(treat_pp) && !is.null(task_obj) && !is.null(task_obj$hawkes_params)) {
    treat_pp <- task_obj$hawkes_params$treated
  }
  c(
    control_mu = safe_mu(ctrl_pp),
    treated_mu = safe_mu(treat_pp)
  )
}
is_failed_mu_fit <- function(r, task_obj, mu_mult) {
  mu_pair <- extract_mu_pair(r, task_obj)
  ctrl_mu <- mu_pair[["control_mu"]]
  treat_mu <- mu_pair[["treated_mu"]]
  if (!is.finite(ctrl_mu) || !is.finite(treat_mu)) return(TRUE)
  (ctrl_mu > (mu_mult * hawkes_par_1$mu)) || (treat_mu > (mu_mult * hawkes_par_2$mu))
}

high_mu_failed_idx <- integer(0)
if (isTRUE(FILTER_HIGH_MU_FITS)) {
  high_mu_failed_idx <- which(vapply(seq_along(results_flat), function(k) {
    is_failed_mu_fit(results_flat[[k]], tasks[[k]], mu_mult = MU_FAIL_MULTIPLIER)
  }, logical(1)))
}
high_mu_failed_tasks <- if (length(high_mu_failed_idx) > 0) {
  do.call(rbind, lapply(high_mu_failed_idx, function(k) {
    mu_pair <- extract_mu_pair(results_flat[[k]], tasks[[k]])
    data.frame(
      task_idx = k,
      sim_id = ((k - 1) %% SIM_SIZE) + 1,
      labelling = tasks[[k]]$labelling_name,
      control_mu = mu_pair[["control_mu"]],
      treated_mu = mu_pair[["treated_mu"]],
      control_mu_threshold = MU_FAIL_MULTIPLIER * hawkes_par_1$mu,
      treated_mu_threshold = MU_FAIL_MULTIPLIER * hawkes_par_2$mu,
      stringsAsFactors = FALSE
    )
  }))
} else {
  data.frame(
    task_idx = integer(0), sim_id = integer(0), labelling = character(0),
    control_mu = numeric(0), treated_mu = numeric(0),
    control_mu_threshold = numeric(0), treated_mu_threshold = numeric(0),
    stringsAsFactors = FALSE
  )
}
include_results <- rep(TRUE, length(results_flat))
if (length(crazy_idx) > 0L) include_results[crazy_idx] <- FALSE
if (length(high_mu_failed_idx) > 0L) include_results[high_mu_failed_idx] <- FALSE
if (length(crazy_idx) > 0L) {
  log_msg(sprintf(
    "[FIT FILTER] %d/%d fits excluded due to explosive/crazy parameters (pre-ATE screen).",
    length(crazy_idx), length(results_flat)
  ))
}
if (length(high_mu_failed_idx) > 0L) {
  log_msg(sprintf(
    "[HIGH MU FILTER] %d/%d fits flagged as failed (mu > %.2fx true); excluded from outputs/estimates but retained in saved raw results.",
    length(high_mu_failed_idx), length(results_flat), MU_FAIL_MULTIPLIER
  ))
} else if (isTRUE(FILTER_HIGH_MU_FITS)) {
  log_msg(sprintf("[HIGH MU FILTER] No fits exceeded %.2fx true mu.", MU_FAIL_MULTIPLIER))
}

# Central fit-status table used by all downstream summaries and plots.
fit_status_df <- do.call(rbind, lapply(seq_along(tasks), function(k) {
  mu_pair <- extract_mu_pair(results_flat[[k]], tasks[[k]])
  data.frame(
    task_idx = k,
    sim_id = ((k - 1L) %% SIM_SIZE) + 1L,
    labelling = tasks[[k]]$labelling_name,
    control_mu = mu_pair[["control_mu"]],
    treated_mu = mu_pair[["treated_mu"]],
    crazy_fit = k %in% crazy_idx,
    high_mu_fit = k %in% high_mu_failed_idx,
    include_result = isTRUE(include_results[[k]]),
    stringsAsFactors = FALSE
  )
}))
if (!is.null(fit_status_df) && nrow(fit_status_df) > 0) {
  exclusion_summary <- fit_status_df %>%
    group_by(.data$labelling) %>%
    summarize(
      n_total = n(),
      n_excluded = sum(!.data$include_result),
      n_high_mu = sum(.data$high_mu_fit),
      n_crazy = sum(.data$crazy_fit),
      .groups = "drop"
    )
  log_msg("=== Fit exclusion summary by method ===")
  print(exclusion_summary)
}
method_kept_idx <- function(labelling_name) {
  if (is.null(fit_status_df) || nrow(fit_status_df) < 1L) return(integer(0))
  fit_status_df %>%
    filter(.data$labelling == labelling_name, .data$include_result) %>%
    pull(.data$sim_id)
}

label_recovery_kept_idx <- function(method_name) {
  fit_label <- if (identical(method_name, "SEM_adaptive")) "SEM_full" else method_name
  method_kept_idx(fit_label)
}

# Recompute classification metrics using only kept fits so paper summaries
# match the high-mu filtering used for ATE and parameter outputs.
label_recovery_all <- label_recovery_metrics(labellings, kept_idx = label_recovery_kept_idx)
class_metrics <- summarize_label_recovery(label_recovery_all)
log_msg("")
log_msg("=== Classification Accuracy (after high-mu filter) ===")
print(class_metrics)

results_df_all <- do.call(rbind, lapply(seq_along(results_flat), function(k) {
  r <- results_flat[[k]]
  if (is.null(r)) return(NULL)
  data.frame(
    task_idx = k,
    sim_id = ((k - 1L) %% SIM_SIZE) + 1L,
    labelling = tasks[[k]]$labelling_name,
    control_mu = extract_mu_pair(r, tasks[[k]])[["control_mu"]],
    treated_mu = extract_mu_pair(r, tasks[[k]])[["treated_mu"]],
    all_nothing_theory = mean(r$all_nothing_theory$ATE),
    tau_1_estim = r$tau_1_estim,
    ATE_total = r$ATE_total,
    ATE_treatment = r$ATE_treatment,
    ATE_spillover = r$ATE_spillover,
    ATE_naive = r$ATE_naive,
    stringsAsFactors = FALSE
  )
}))
results_df <- if (!is.null(results_df_all) && nrow(results_df_all) > 0) {
  results_df_all %>%
    left_join(fit_status_df %>% select(.data$task_idx, .data$include_result), by = "task_idx") %>%
    filter(!is.na(.data$include_result) & .data$include_result) %>%
    select(-.data$include_result)
} else {
  results_df_all
}

log_msg("")
log_msg("=== ATE Results ===")
if (!is.null(results_df) && nrow(results_df) > 0) {
  summary_df <- results_df %>%
    group_by(.data$labelling) %>%
    summarize(
      mean_all_nothing = mean(.data$all_nothing_theory, na.rm = TRUE),
      mean_tau1 = mean(.data$tau_1_estim, na.rm = TRUE),
      mean_ATE_total = mean(.data$ATE_total, na.rm = TRUE),
      mean_ATE_naive = mean(.data$ATE_naive, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_df)
} else {
  log_msg("No valid ATE results (all tasks timed out or errored).")
}

all_nothing_dtaite_per_unit <- function(control_pp, treated_pp, n_units = partition$n) {
  (as.numeric(treated_pp$mu) * TIME_INT * (1 / max(1e-8, 1 - as.numeric(treated_pp$K))) -
     as.numeric(control_pp$mu) * TIME_INT * (1 / max(1e-8, 1 - as.numeric(control_pp$K)))) / n_units
}

TRUE_CTRL_TAU_I <- suppressWarnings(as.integer(Sys.getenv("PP_TRUECTRL_TAU_I", "3")))
if (!is.finite(TRUE_CTRL_TAU_I) || is.na(TRUE_CTRL_TAU_I) || TRUE_CTRL_TAU_I < 1L) TRUE_CTRL_TAU_I <- 3L
TRUE_CTRL_TAU_SIMS <- suppressWarnings(as.integer(Sys.getenv("PP_TRUECTRL_TAU_SIMS", "5")))
if (!is.finite(TRUE_CTRL_TAU_SIMS) || is.na(TRUE_CTRL_TAU_SIMS) || TRUE_CTRL_TAU_SIMS < 1L) TRUE_CTRL_TAU_SIMS <- 5L

results_df_true_control_all <- do.call(rbind, lapply(seq_along(results_flat), function(k) {
  r <- results_flat[[k]]
  if (is.null(r) || is.null(r$treated_pp)) return(NULL)
  treated_pp <- r$treated_pp
  ctrl_true <- hawkes_par_1
  tau_i_true <- if (isTRUE(ATE_COMPUTE_TAU)) {
    vapply(seq_len(TRUE_CTRL_TAU_I), function(j) {
      tau_i(
        sample(partition$n, 1),
        partition = partition, treated_partitions = treated_partitions,
        statespace = OMEGA, windowT = c(TREATMENT_TIME, END_TIME),
        control_pp = ctrl_true, treated_pp = treated_pp, n_sim = TRUE_CTRL_TAU_SIMS
      )
    }, numeric(1))
  } else {
    NA_real_
  }
  data.frame(
    task_idx = k,
    sim_id = ((k - 1L) %% SIM_SIZE) + 1L,
    labelling = tasks[[k]]$labelling_name,
    control_mu = extract_mu_pair(r, tasks[[k]])[["control_mu"]],
    treated_mu = extract_mu_pair(r, tasks[[k]])[["treated_mu"]],
    all_nothing_true_control = all_nothing_dtaite_per_unit(ctrl_true, treated_pp),
    tau_1_true_control = if (all(is.na(tau_i_true))) NA_real_ else mean(tau_i_true, na.rm = TRUE),
    ATE_total = r$ATE_total,
    ATE_treatment = r$ATE_treatment,
    ATE_spillover = r$ATE_spillover,
    ATE_naive = r$ATE_naive,
    stringsAsFactors = FALSE
  )
}))
results_df_true_control <- if (!is.null(results_df_true_control_all) && nrow(results_df_true_control_all) > 0) {
  results_df_true_control_all %>%
    left_join(fit_status_df %>% select(.data$task_idx, .data$include_result), by = "task_idx") %>%
    filter(!is.na(.data$include_result) & .data$include_result) %>%
    select(-.data$include_result)
} else {
  results_df_true_control_all
}
summary_df_true_control <- if (!is.null(results_df_true_control) && nrow(results_df_true_control) > 0) {
  results_df_true_control %>%
    group_by(.data$labelling) %>%
    summarize(
      mean_all_nothing_true_control = mean(.data$all_nothing_true_control, na.rm = TRUE),
      mean_tau1_true_control = mean(.data$tau_1_true_control, na.rm = TRUE),
      mean_ATE_total = mean(.data$ATE_total, na.rm = TRUE),
      mean_ATE_naive = mean(.data$ATE_naive, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  NULL
}
if (!is.null(summary_df_true_control)) {
  log_msg("=== ATE Results (true control parameters fixed) ===")
  print(summary_df_true_control)
}

# ------------------------------------------------------------------
# 11b. Off-support allocation contrasts
# ------------------------------------------------------------------
expected_count_per_tile_hawkes <- function(z, control_pp, treated_pp) {
  z <- as.logical(z)
  n_treated <- sum(z, na.rm = TRUE)
  n_control <- length(z) - n_treated
  control_rate <- as.numeric(control_pp$mu) / max(1e-8, 1 - as.numeric(control_pp$K))
  treated_rate <- as.numeric(treated_pp$mu) / max(1e-8, 1 - as.numeric(treated_pp$K))
  TIME_INT * (n_control * control_rate + n_treated * treated_rate) / partition$n
}

make_support_contrasts <- function() {
  n_cells <- partition$n
  n_treated_50 <- as.integer(round(n_cells * TREAT_PROP))
  n_treated_50 <- max(1L, min(n_cells - 1L, n_treated_50))
  set.seed(stage_seed(5L, 0L, 1L))
  z_random_50 <- rep(FALSE, n_cells)
  z_random_50[sample(seq_len(n_cells), n_treated_50)] <- TRUE
  global_all_treated <- rep(TRUE, n_cells)
  global_all_control <- rep(FALSE, n_cells)

  single_cell_rows <- lapply(seq_len(n_cells), function(j) {
    z_single <- rep(FALSE, n_cells)
    z_single[j] <- TRUE
    data.frame(
      contrast_id = paste0("single_cell_vs_all_control_", j),
      contrast_family = "single_cell_flip",
      cell_id = j,
      hamming_distance_from_zobs = NA_integer_,
      psi_truth = expected_count_per_tile_hawkes(z_single, hawkes_par_1, hawkes_par_2) -
        expected_count_per_tile_hawkes(global_all_control, hawkes_par_1, hawkes_par_2),
      stringsAsFactors = FALSE
    )
  })

  rbind(
    do.call(rbind, single_cell_rows),
    data.frame(
      contrast_id = "random_50pct_vs_all_control",
      contrast_family = "random_50pct_flip",
      cell_id = NA_integer_,
      hamming_distance_from_zobs = NA_integer_,
      psi_truth = expected_count_per_tile_hawkes(z_random_50, hawkes_par_1, hawkes_par_2) -
        expected_count_per_tile_hawkes(global_all_control, hawkes_par_1, hawkes_par_2),
      stringsAsFactors = FALSE
    ),
    data.frame(
      contrast_id = "global_all_treated_vs_all_control",
      contrast_family = "global_1_0",
      cell_id = NA_integer_,
      hamming_distance_from_zobs = NA_integer_,
      psi_truth = expected_count_per_tile_hawkes(global_all_treated, hawkes_par_1, hawkes_par_2) -
        expected_count_per_tile_hawkes(global_all_control, hawkes_par_1, hawkes_par_2),
      stringsAsFactors = FALSE
    )
  )
}

support_contrast_specs <- make_support_contrasts()
estimate_support_contrast_for_task <- function(k, spec_row) {
  r <- results_flat[[k]]
  task_k <- tasks[[k]]
  treat_fit <- if (!is.null(r) && !is.null(r$treated_pp)) r$treated_pp else NULL
  if (is.null(treat_fit)) return(NULL)
  ctrl_true <- hawkes_par_1
  if (identical(spec_row$contrast_family, "global_1_0")) {
    return(data.frame(
      task_idx = k,
      sim_id = ((k - 1L) %% SIM_SIZE) + 1L,
      labelling = task_k$labelling_name,
      contrast_id = spec_row$contrast_id,
      contrast_family = spec_row$contrast_family,
      cell_id = spec_row$cell_id,
      hamming_distance_from_zobs = spec_row$hamming_distance_from_zobs,
      psi_truth = spec_row$psi_truth,
      psi_estimate = all_nothing_dtaite_per_unit(ctrl_true, treat_fit),
      stringsAsFactors = FALSE
    ))
  }
  n_cells <- length(partition_processes)
  z_a <- switch(
    spec_row$contrast_family,
    single_cell_flip = {
      z <- rep(FALSE, n_cells)
      z[as.integer(spec_row$cell_id)] <- TRUE
      z
    },
    random_50pct_flip = {
      n_treated_50 <- as.integer(round(n_cells * TREAT_PROP))
      n_treated_50 <- max(1L, min(n_cells - 1L, n_treated_50))
      set.seed(stage_seed(5L, 0L, 1L))
      z <- rep(FALSE, n_cells)
      z[sample(seq_len(n_cells), n_treated_50)] <- TRUE
      z
    },
    global_1_0 = rep(TRUE, n_cells)
  )
  z_b <- rep(FALSE, n_cells)
  data.frame(
    task_idx = k,
    sim_id = ((k - 1L) %% SIM_SIZE) + 1L,
    labelling = task_k$labelling_name,
    contrast_id = spec_row$contrast_id,
    contrast_family = spec_row$contrast_family,
    cell_id = spec_row$cell_id,
    hamming_distance_from_zobs = spec_row$hamming_distance_from_zobs,
    psi_truth = spec_row$psi_truth,
    psi_estimate = expected_count_per_tile_hawkes(z_a, ctrl_true, treat_fit) -
      expected_count_per_tile_hawkes(z_b, ctrl_true, treat_fit),
    stringsAsFactors = FALSE
  )
}

support_contrast_df_all <- do.call(rbind, lapply(seq_along(tasks), function(k) {
  do.call(rbind, lapply(seq_len(nrow(support_contrast_specs)), function(j) {
    estimate_support_contrast_for_task(k, support_contrast_specs[j, , drop = FALSE])
  }))
}))
support_contrast_df <- if (!is.null(support_contrast_df_all) && nrow(support_contrast_df_all) > 0) {
  global_ate_truth_ref <- {
    hit <- support_contrast_specs$psi_truth[support_contrast_specs$contrast_family == "global_1_0"]
    if (length(hit) < 1L) NA_real_ else hit[[1L]]
  }
  support_contrast_df_all %>%
    left_join(fit_status_df %>% select(.data$task_idx, .data$include_result), by = "task_idx") %>%
    filter(!is.na(.data$include_result) & .data$include_result) %>%
    mutate(
      abs_error = abs(.data$psi_estimate - .data$psi_truth),
      pct_error = ifelse(
        is.finite(.data$psi_truth) & abs(.data$psi_truth) > 1e-8,
        100 * .data$abs_error / abs(.data$psi_truth),
        ifelse(
          .data$contrast_family == "single_cell_flip" &
            is.finite(global_ate_truth_ref) & abs(global_ate_truth_ref) > 1e-8,
          100 * .data$abs_error / abs(global_ate_truth_ref),
          NA_real_
        )
      )
    ) %>%
    select(-.data$include_result)
} else {
  support_contrast_df_all
}
support_contrast_summary <- if (!is.null(support_contrast_df) && nrow(support_contrast_df) > 0) {
  support_contrast_df %>%
    mutate(labelling = normalize_reported_labelling(.data$labelling)) %>%
    group_by(.data$labelling, .data$contrast_family) %>%
    summarize(
      n = n(),
      mean_psi_estimate = mean(.data$psi_estimate, na.rm = TRUE),
      mean_psi_truth = mean(.data$psi_truth, na.rm = TRUE),
      mean_abs_error = mean(.data$abs_error, na.rm = TRUE),
      mean_pct_error = mean(.data$pct_error, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  NULL
}
if (!is.null(support_contrast_summary)) {
  log_msg("=== Off-support allocation contrasts ===")
  print(support_contrast_summary)
}

# ------------------------------------------------------------------
# 11c. Decay-validation diagnostic (forward simulation only)
# ------------------------------------------------------------------
source(file.path(SCRIPT_DIR, "decay_validation_utils.R"), local = FALSE)

run_decay_validation <- function() {
  if (!isTRUE(RUN_DECAY_VALIDATION)) return(NULL)
  out <- run_decay_validation_scenario(
    omega = OMEGA,
    nx = NX,
    ny = NY,
    treatment_time = TREATMENT_TIME,
    end_time = END_TIME,
    partition_processes = partition_processes,
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    sim_kernel = SIM_KERNEL,
    decay_reps = DECAY_VALIDATION_REPS,
    annulus_width = DECAY_ANNULUS_WIDTH,
    flip_cell = DECAY_FLIP_CELL,
    stage_seed_fn = stage_seed,
    log_fn = log_msg
  )
  if (is.null(out)) return(NULL)
  decay_summary <- out$summary
  rate_lines <- out$rate_lines
  eps <- min(decay_summary$mean_abs_delta[decay_summary$mean_abs_delta > 0], na.rm = TRUE)
  if (!is.finite(eps)) eps <- 1e-6
  eps <- eps / 2
  decay_plot <- ggplot(decay_summary, aes(x = .data$d_mid, y = .data$mean_abs_delta_plot)) +
    geom_ribbon(aes(ymin = pmax(.data$q10_abs_delta, eps),
                    ymax = pmax(.data$q90_abs_delta, eps)),
                alpha = 0.18, fill = "#0072B2") +
    geom_point(size = 1.5, color = "#0072B2") +
    geom_line(linewidth = 0.7, color = "#0072B2") +
    {if (!is.null(rate_lines) && nrow(rate_lines) > 0)
      geom_line(data = rate_lines,
                aes(x = .data$d_mid, y = .data$rate_value_plot,
                    color = .data$reference),
                linetype = "dashed", linewidth = 0.8, inherit.aes = FALSE)} +
    scale_y_log10() +
    labs(
      title = "Allocation-influence decay from a single-cell flip",
      subtitle = paste0(
        "Forward simulation only: one treated cell flipped to control under CRN; ",
        "mean |Delta N| per distance annulus."
      ),
      x = "Distance from flipped cell (unit-width annuli)",
      y = "Mean |Delta N| per annulus",
      color = "Reference slope"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  list(
    specs = out$specs,
    df = out$df,
    summary = decay_summary,
    rate_lines = rate_lines,
    plot = decay_plot
  )
}

decay_validation <- run_decay_validation()
if (!is.null(decay_validation)) {
  log_msg("=== Decay-validation annulus summary ===")
  print(utils::head(decay_validation$summary, 12))
}

log_msg("")
log_msg("True all/nothing ATE:", round(all_nothing_ATE, 4))
log_msg("True one-flip ATE:  ", round(true_tau_1, 4))

# ------------------------------------------------------------------
# 12. Create ggplots (as in old sim study) and save as elements of results
# ------------------------------------------------------------------
log_msg("")
log_msg("Building plots ...")
sim_study_plots <- list()
if (!is.null(decay_validation) && inherits(decay_validation$plot, "ggplot")) {
  sim_study_plots$plot_decay_validation <- decay_validation$plot
}
boxplot_method_levels <- c("oracle", "naive", "best", "SEM_adaptive", "SEM_full")

subset_core_methods <- function(df, label_col = "labelling") {
  if (is.null(df) || nrow(df) < 1L || !(label_col %in% names(df))) return(df[0, , drop = FALSE])
  out <- df %>% filter(.data[[label_col]] %in% c("oracle", "naive", "SEM_full"))
  out[[label_col]] <- as.character(out[[label_col]])
  out[[label_col]][out[[label_col]] == "SEM_full"] <- "SEM"
  out[[label_col]] <- factor(out[[label_col]], levels = c("oracle", "naive", "SEM"))
  out
}

strip_titles_from_plots <- function(x) {
  if (inherits(x, "ggplot")) {
    return(
      x +
        labs(title = NULL, subtitle = NULL) +
        theme(plot.title = element_blank(), plot.subtitle = element_blank())
    )
  }
  if (is.list(x)) return(lapply(x, strip_titles_from_plots))
  x
}

# Point pattern plots (first realization)
if (length(obs_data) > 0) {
  post_first <- as.data.frame(obs_data[[1]]) %>% filter(.data$t > TREATMENT_TIME)
  pre_first  <- as.data.frame(obs_data[[1]]) %>% filter(.data$t < TREATMENT_TIME)
  if (nrow(post_first) > 0) {
    sim_study_plots$plot_pp_post <- plot_pp(post_first, partition = partition,
      title = "")
  }
  if (nrow(pre_first) > 0) {
    sim_study_plots$plot_pp_pre <- plot_pp(pre_first, partition = partition,
      title = "")
  }
}

# Control and treated parameter estimates from results_flat
compute_post_loglik <- function(task_obj, par_obj, process_name) {
  if (is.null(task_obj) || is.null(par_obj) || is.null(task_obj$x)) return(NA_real_)
  post_df <- as.data.frame(task_obj$x)
  if (!all(c("x", "y", "t", "inferred_process") %in% names(post_df))) return(NA_real_)
  proc_df <- post_df[post_df$inferred_process == process_name, , drop = FALSE]
  if (nrow(proc_df) < 1L) return(NA_real_)
  zero_bg_region <- if (identical(process_name, "control")) treated_state_space else control_state_space
  val <- tryCatch(
    loglik_hawk_fast(
      params = c(par_obj$mu, par_obj$alpha, par_obj$beta, par_obj$K),
      realiz = proc_df,
      windowT = c(TREATMENT_TIME, END_TIME),
      windowS = OMEGA,
      zero_background_region = zero_bg_region,
      poisson_flag = FALSE,
      optimized = TRUE
    ),
    error = function(e) NA_real_
  )
  as.numeric(val)
}

extract_param_rows <- function(results_flat, tasks, sim_size, field_name) {
  lapply(seq_along(results_flat), function(k) {
    r <- results_flat[[k]]
    task_k <- tasks[[k]]
    par_obj <- if (!is.null(r)) r[[field_name]] else NULL
    # If SEM_full ATE was skipped (e.g., explosive guard), still expose
    # fitted SEM parameters in parameter plots for transparency.
    if (is.null(par_obj) &&
        !is.null(task_k) &&
        identical(task_k$labelling_name, "SEM_full") &&
        !is.null(task_k$hawkes_params)) {
      par_obj <- if (identical(field_name, "control_pp")) task_k$hawkes_params$control else task_k$hawkes_params$treated
    }
    if (is.null(par_obj)) return(NULL)
    process_name <- if (identical(field_name, "control_pp")) "control" else "treated"
    data.frame(
      task_idx = k,
      labelling = task_k$labelling_name,
      sim_id = ((k - 1) %% sim_size) + 1,
      mu = par_obj$mu, alpha = par_obj$alpha,
      beta = par_obj$beta, K = par_obj$K,
      loglik = compute_post_loglik(task_k, par_obj, process_name),
      stringsAsFactors = FALSE
    )
  })
}
extract_oracle_param_means <- function(params_df) {
  if (is.null(params_df) || nrow(params_df) < 1) return(NULL)
  oracle_df <- params_df %>% filter(.data$labelling == "oracle")
  if (nrow(oracle_df) < 1) return(NULL)
  out <- c(
    mu = mean(oracle_df$mu, na.rm = TRUE),
    alpha = mean(oracle_df$alpha, na.rm = TRUE),
    beta = mean(oracle_df$beta, na.rm = TRUE),
    K = mean(oracle_df$K, na.rm = TRUE)
  )
  if (all(!is.finite(out))) return(NULL)
  out
}
build_param_boxplots <- function(params_df, truth_params, oracle_param_means = NULL) {
  if (is.null(params_df) || nrow(params_df) < 1) return(NULL)
  params_df$labelling <- factor(
    params_df$labelling,
    levels = unique(c(boxplot_method_levels, as.character(params_df$labelling)))
  )
  params_long <- reshape2::melt(
    params_df,
    id.vars = c("labelling", "sim_id"),
    variable.name = "param",
    value.name = "value"
  )
  plots <- lapply(c("mu", "alpha", "beta", "K"), function(p) {
    tmp <- params_long %>% filter(.data$param == p)
    if (nrow(tmp) == 0) return(NULL)
    oracle_mean <- NA_real_
    if (!is.null(oracle_param_means) && p %in% names(oracle_param_means)) {
      oracle_mean <- as.numeric(oracle_param_means[[p]])
    }
    if (!is.finite(oracle_mean)) {
      oracle_vals <- tmp %>% filter(.data$labelling == "oracle")
      if (nrow(oracle_vals) > 0) oracle_mean <- mean(oracle_vals$value, na.rm = TRUE)
    }
    ggplot(tmp) +
      geom_boxplot(aes(x = .data$labelling, y = .data$value)) +
      # Show sample means so the plot aligns with the mean/sd summary tables.
      stat_summary(aes(x = .data$labelling, y = .data$value), fun = mean,
                   geom = "point", shape = 23, size = 2, fill = "black", color = "black") +
      geom_hline(yintercept = truth_params[[p]], linetype = "dashed", color = "red") +
      {if (is.finite(oracle_mean)) geom_hline(yintercept = oracle_mean, linetype = "dotted", color = "blue") } +
      labs(x = "Method", y = paste0(p, " estimate")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  Filter(Negate(is.null), plots)
}

build_mu_k_ridge_plot <- function(params_df, truth_params, process_name) {
  if (is.null(params_df) || nrow(params_df) < 2) return(NULL)
  dat <- params_df %>%
    filter(is.finite(.data$mu), is.finite(.data$K))
  if (nrow(dat) < 2) return(NULL)
  dat$labelling <- factor(
    dat$labelling,
    levels = unique(c(boxplot_method_levels, as.character(dat$labelling)))
  )
  p <- ggplot(dat, aes(x = .data$K, y = .data$mu))
  if ("loglik" %in% names(dat) && any(is.finite(dat$loglik))) {
    p <- p +
      geom_point(aes(color = .data$loglik), alpha = 0.65, size = 1.9) +
      scale_color_viridis_c(option = "C", na.value = "gray75", name = "Log-likelihood")
  } else {
    p <- p + geom_point(alpha = 0.5, size = 1.8, color = "#1f77b4")
  }
  p +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "#d62728") +
    geom_vline(xintercept = truth_params$K, linetype = "dashed", color = "gray35") +
    geom_hline(yintercept = truth_params$mu, linetype = "dashed", color = "gray35") +
    facet_wrap(~ labelling, scales = "free") +
    labs(
      x = "K estimate",
      y = "mu estimate",
      title = paste0(process_name, ": mu vs K by fit method"),
      subtitle = "Point color shows fitted post-treatment log-likelihood"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 9),
      plot.title = element_text(size = 11, face = "bold")
    )
}

control_param_rows <- extract_param_rows(results_flat, tasks, SIM_SIZE, "control_pp")
treated_param_rows <- extract_param_rows(results_flat, tasks, SIM_SIZE, "treated_pp")
control_params_df_all <- do.call(rbind, control_param_rows)
treated_params_df_all <- do.call(rbind, treated_param_rows)
control_params_df <- if (!is.null(control_params_df_all) && nrow(control_params_df_all) > 0) {
  control_params_df_all %>%
    left_join(fit_status_df %>% select(.data$task_idx, .data$include_result), by = "task_idx") %>%
    filter(!is.na(.data$include_result) & .data$include_result) %>%
    select(-.data$include_result)
} else {
  control_params_df_all
}
treated_params_df <- if (!is.null(treated_params_df_all) && nrow(treated_params_df_all) > 0) {
  treated_params_df_all %>%
    left_join(fit_status_df %>% select(.data$task_idx, .data$include_result), by = "task_idx") %>%
    filter(!is.na(.data$include_result) & .data$include_result) %>%
    select(-.data$include_result)
} else {
  treated_params_df_all
}
if (!is.null(treated_params_df) && nrow(treated_params_df) > 0) {
  present_treated <- unique(as.character(treated_params_df$labelling))
  if (!("oracle" %in% present_treated)) {
    log_msg("[WARN] Treated parameter rows contain no oracle fits.")
  }
}

summarize_param_table <- function(df, process_name) {
  if (is.null(df) || nrow(df) < 1) return(NULL)
  out <- df %>%
    group_by(.data$labelling) %>%
    summarize(
      n = n(),
      mu_mean = mean(.data$mu, na.rm = TRUE), mu_sd = sd(.data$mu, na.rm = TRUE),
      alpha_mean = mean(.data$alpha, na.rm = TRUE), alpha_sd = sd(.data$alpha, na.rm = TRUE),
      beta_mean = mean(.data$beta, na.rm = TRUE), beta_sd = sd(.data$beta, na.rm = TRUE),
      K_mean = mean(.data$K, na.rm = TRUE), K_sd = sd(.data$K, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(factor(.data$labelling, levels = boxplot_method_levels))
  out$process <- process_name
  out
}
control_param_summary <- summarize_param_table(control_params_df, "control")
treated_param_summary <- summarize_param_table(treated_params_df, "treated")
control_oracle_means <- extract_oracle_param_means(control_params_df)
treated_oracle_means <- extract_oracle_param_means(treated_params_df)
if (!is.null(control_param_summary)) {
  log_msg("=== Control fitted parameters (mean/sd) ===")
  print(control_param_summary)
}
if (!is.null(treated_param_summary)) {
  log_msg("=== Treated fitted parameters (mean/sd) ===")
  print(treated_param_summary)
}

control_param_plots <- build_param_boxplots(
  control_params_df,
  hawkes_par_1,
  oracle_param_means = control_oracle_means
)
if (!is.null(control_param_plots) && length(control_param_plots) > 0) {
  sim_study_plots$plot_control_params <- control_param_plots
}

treated_param_plots <- build_param_boxplots(
  treated_params_df,
  hawkes_par_2,
  oracle_param_means = treated_oracle_means
)
if (!is.null(treated_param_plots) && length(treated_param_plots) > 0) {
  sim_study_plots$plot_treated_params <- treated_param_plots
}

control_mu_k_plot <- build_mu_k_ridge_plot(control_params_df, hawkes_par_1, "Control")
if (!is.null(control_mu_k_plot)) {
  sim_study_plots$plot_control_mu_k <- control_mu_k_plot
}
treated_mu_k_plot <- build_mu_k_ridge_plot(treated_params_df, hawkes_par_2, "Treated")
if (!is.null(treated_mu_k_plot)) {
  sim_study_plots$plot_treated_mu_k <- treated_mu_k_plot
}

# All-nothing ATE boxplot
if (!is.null(results_df) && nrow(results_df) > 0) {
  results_df$labelling <- factor(results_df$labelling,
    levels = unique(c(boxplot_method_levels, as.character(results_df$labelling))))
  lines_data_ate <- data.frame(all_nothing_ATE = all_nothing_ATE)
  oracle_ate <- results_df %>% filter(.data$labelling == "oracle")
  oracle_mean_ate <- if (nrow(oracle_ate) > 0) mean(oracle_ate$all_nothing_theory, na.rm = TRUE) else NA
  y_lo <- min(all_nothing_ATE * 1.5, min(results_df$all_nothing_theory, na.rm = TRUE) * 0.95)
  y_hi <- max(3, max(results_df$all_nothing_theory, na.rm = TRUE) * 1.05)
  sim_study_plots$plot_all_nothing_ATE <- ggplot(results_df) +
    geom_boxplot(aes(x = .data$labelling, y = .data$all_nothing_theory)) +
    geom_hline(data = lines_data_ate, aes(yintercept = .data$all_nothing_ATE),
               linetype = "solid", color = scales::hue_pal()(3)[1], linewidth = 1) +
    {if (!is.na(oracle_mean_ate)) geom_hline(yintercept = oracle_mean_ate, linetype = "dotted",
               color = "blue", linewidth = 0.8) } +
    labs(x = "Method", y = "All-Nothing ATE Estimate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(y_lo, y_hi))
  results_df_core <- subset_core_methods(results_df, "labelling")
  if (!is.null(results_df_core) && nrow(results_df_core) > 0) {
    oracle_ate_core <- results_df_core %>% filter(.data$labelling == "oracle")
    oracle_mean_ate_core <- if (nrow(oracle_ate_core) > 0) mean(oracle_ate_core$all_nothing_theory, na.rm = TRUE) else NA_real_
    sim_study_plots$plot_all_nothing_ATE_core <- ggplot(results_df_core) +
      geom_boxplot(aes(x = .data$labelling, y = .data$all_nothing_theory)) +
      geom_hline(data = lines_data_ate, aes(yintercept = .data$all_nothing_ATE),
        linetype = "solid", color = scales::hue_pal()(3)[1], linewidth = 1) +
      { if (!is.na(oracle_mean_ate_core)) geom_hline(yintercept = oracle_mean_ate_core, linetype = "dotted",
        color = "blue", linewidth = 0.8) } +
      labs(x = "Method", y = "All-Nothing ATE Estimate") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      coord_cartesian(ylim = c(y_lo, y_hi))
  }
}

if (!is.null(results_df_true_control) && nrow(results_df_true_control) > 0) {
  results_df_true_control$labelling <- factor(
    results_df_true_control$labelling,
    levels = unique(c(boxplot_method_levels, as.character(results_df_true_control$labelling)))
  )
  oracle_tc <- results_df_true_control %>% filter(.data$labelling == "oracle")
  oracle_tc_mean <- if (nrow(oracle_tc) > 0) mean(oracle_tc$all_nothing_true_control, na.rm = TRUE) else NA_real_
  y_lo_tc <- min(all_nothing_ATE * 1.5, min(results_df_true_control$all_nothing_true_control, na.rm = TRUE) * 0.95)
  y_hi_tc <- max(3, max(results_df_true_control$all_nothing_true_control, na.rm = TRUE) * 1.05)
  sim_study_plots$plot_all_nothing_ATE_true_control <- ggplot(results_df_true_control) +
    geom_boxplot(aes(x = .data$labelling, y = .data$all_nothing_true_control)) +
    geom_hline(yintercept = all_nothing_ATE, linetype = "solid", color = scales::hue_pal()(3)[1], linewidth = 1) +
    {if (!is.na(oracle_tc_mean)) geom_hline(yintercept = oracle_tc_mean, linetype = "dotted", color = "blue", linewidth = 0.8)} +
    labs(x = "Method", y = "All-Nothing ATE (True Control Fixed)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(y_lo_tc, y_hi_tc))
  results_df_true_control_core <- subset_core_methods(results_df_true_control, "labelling")
  if (!is.null(results_df_true_control_core) && nrow(results_df_true_control_core) > 0) {
    oracle_tc_core <- results_df_true_control_core %>% filter(.data$labelling == "oracle")
    oracle_tc_core_mean <- if (nrow(oracle_tc_core) > 0) mean(oracle_tc_core$all_nothing_true_control, na.rm = TRUE) else NA_real_
    sim_study_plots$plot_all_nothing_ATE_true_control_core <- ggplot(results_df_true_control_core) +
      geom_boxplot(aes(x = .data$labelling, y = .data$all_nothing_true_control)) +
      geom_hline(yintercept = all_nothing_ATE, linetype = "solid", color = scales::hue_pal()(3)[1], linewidth = 1) +
      { if (!is.na(oracle_tc_core_mean)) geom_hline(yintercept = oracle_tc_core_mean, linetype = "dotted",
        color = "blue", linewidth = 0.8) } +
      labs(x = "Method", y = "All-Nothing ATE (True Control Fixed)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      coord_cartesian(ylim = c(y_lo_tc, y_hi_tc))
  }
}

# Points per tile (control vs treated theoretical means)
ate_detail_rows <- lapply(seq_along(results_flat), function(k) {
  r <- results_flat[[k]]
  if (is.null(r) || is.null(r$all_nothing_theory)) return(NULL)
  th <- r$all_nothing_theory
  data.frame(
    task_idx = rep(k, 2),
    labelling = rep(tasks[[k]]$labelling_name, 2),
    method = c("points_per_tile_control_theory", "points_per_tile_treated_theory"),
    ATE_estim = c(th$c_mean[1], th$t_mean[1]),
    stringsAsFactors = FALSE
  )
})
ate_detail_rows <- do.call(rbind, ate_detail_rows)
if (!is.null(ate_detail_rows) && nrow(ate_detail_rows) > 0) {
  ate_detail_rows <- ate_detail_rows %>%
    left_join(fit_status_df %>% select(.data$task_idx, .data$include_result), by = "task_idx") %>%
    filter(!is.na(.data$include_result) & .data$include_result) %>%
    select(-.data$include_result)
}
if (!is.null(ate_detail_rows) && nrow(ate_detail_rows) > 0) {
  ate_detail_rows$labelling <- factor(
    ate_detail_rows$labelling,
    levels = unique(c(boxplot_method_levels, as.character(ate_detail_rows$labelling)))
  )
  true_means_pts <- data.frame(
    method = c("points_per_tile_control_theory", "points_per_tile_treated_theory"),
    mean_ATE = c(
      hawkes_par_1$mu * TIME_INT * (1 / (1 - hawkes_par_1$K)) / partition$n,
      hawkes_par_2$mu * TIME_INT * (1 / (1 - hawkes_par_2$K)) / partition$n
    )
  )
  oracle_means_pts <- ate_detail_rows %>%
    filter(.data$labelling == "oracle") %>%
    group_by(.data$method) %>%
    summarize(mean_ATE = mean(.data$ATE_estim, na.rm = TRUE), .groups = "drop")
  sim_study_plots$plot_points_per_tile <- ggplot(ate_detail_rows, aes(x = .data$labelling, y = .data$ATE_estim, fill = .data$method)) +
    geom_boxplot() +
    facet_wrap(~ method) +
    geom_hline(aes(yintercept = .data$mean_ATE, colour = "Oracle mean"), data = oracle_means_pts, linetype = "dashed", linewidth = 1) +
    geom_hline(aes(yintercept = .data$mean_ATE, colour = "True mean"), data = true_means_pts, linetype = "dashed", linewidth = 1) +
    scale_colour_manual(name = "Reference", values = c("Oracle mean" = "#0072B2", "True mean" = "#D55E00")) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
    labs(x = "Method", y = "Points Per Tile") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
    coord_cartesian(ylim = c(0, max(true_means_pts$mean_ATE) * 1.05))
  ate_detail_rows_core <- subset_core_methods(ate_detail_rows, "labelling")
  if (!is.null(ate_detail_rows_core) && nrow(ate_detail_rows_core) > 0) {
    oracle_means_pts_core <- ate_detail_rows_core %>%
      filter(.data$labelling == "oracle") %>%
      group_by(.data$method) %>%
      summarize(mean_ATE = mean(.data$ATE_estim, na.rm = TRUE), .groups = "drop")
    sim_study_plots$plot_points_per_tile_core <- ggplot(ate_detail_rows_core,
      aes(x = .data$labelling, y = .data$ATE_estim, fill = .data$method)) +
      geom_boxplot() +
      facet_wrap(~ method) +
      geom_hline(aes(yintercept = .data$mean_ATE, colour = "Oracle mean"), data = oracle_means_pts_core, linetype = "dashed", linewidth = 1) +
      geom_hline(aes(yintercept = .data$mean_ATE, colour = "True mean"), data = true_means_pts, linetype = "dashed", linewidth = 1) +
      scale_colour_manual(name = "Reference", values = c("Oracle mean" = "#0072B2", "True mean" = "#D55E00")) +
      guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
      labs(x = "Method", y = "Points Per Tile") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
      coord_cartesian(ylim = c(0, max(true_means_pts$mean_ATE) * 1.05))
  }
}

# One-flip ATE boxplot
if (!is.null(results_df) && nrow(results_df) > 0) {
  lines_data_tau <- data.frame(true_1_flip = true_tau_1)
  oracle_tau <- results_df %>% filter(.data$labelling == "oracle")
  oracle_mean_tau <- if (nrow(oracle_tau) > 0) mean(oracle_tau$tau_1_estim, na.rm = TRUE) else NA
  sim_study_plots$plot_one_flip_ATE <- ggplot(results_df) +
    geom_boxplot(aes(x = .data$labelling, y = .data$tau_1_estim)) +
    geom_hline(data = lines_data_tau, aes(yintercept = .data$true_1_flip),
               linetype = "solid", color = scales::hue_pal()(3)[1], linewidth = 1) +
    {if (!is.na(oracle_mean_tau)) geom_hline(yintercept = oracle_mean_tau, linetype = "dotted",
               color = "blue", linewidth = 0.8) } +
    labs(x = "Method", y = "Single Flip ATE Estimate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  results_df_core_tau <- subset_core_methods(results_df, "labelling")
  if (!is.null(results_df_core_tau) && nrow(results_df_core_tau) > 0) {
    oracle_tau_core <- results_df_core_tau %>% filter(.data$labelling == "oracle")
    oracle_mean_tau_core <- if (nrow(oracle_tau_core) > 0) mean(oracle_tau_core$tau_1_estim, na.rm = TRUE) else NA_real_
    sim_study_plots$plot_one_flip_ATE_core <- ggplot(results_df_core_tau) +
      geom_boxplot(aes(x = .data$labelling, y = .data$tau_1_estim)) +
      geom_hline(data = lines_data_tau, aes(yintercept = .data$true_1_flip),
        linetype = "solid", color = scales::hue_pal()(3)[1], linewidth = 1) +
      { if (!is.na(oracle_mean_tau_core)) geom_hline(yintercept = oracle_mean_tau_core, linetype = "dotted",
        color = "blue", linewidth = 0.8) } +
      labs(x = "Method", y = "Single Flip ATE Estimate") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
}

# Draw helpers (grid.arrange objects don't survive save/load)
sim_study_plots$draw_control_params_combined <- function() {
  if (length(sim_study_plots$plot_control_params) > 0)
    gridExtra::grid.arrange(grobs = sim_study_plots$plot_control_params, ncol = 2)
}
sim_study_plots$draw_treated_params_combined <- function() {
  if (length(sim_study_plots$plot_treated_params) > 0)
    gridExtra::grid.arrange(grobs = sim_study_plots$plot_treated_params, ncol = 2)
}

if (exists("EM_results") && length(EM_results) > 0) {
  sem_keep_idx <- method_kept_idx("SEM_adaptive")
  acc_list <- lapply(sem_keep_idx, function(i) {
    acc <- EM_results[[i]]$adaptive$accuracies
    if (is.null(acc) || length(acc) == 0) return(NULL)
    data.frame(
      iteration = seq_along(acc),
      accuracy = acc,
      sim_id = i
    )
  })
  acc_list <- Filter(Negate(is.null), acc_list)
  acc_df <- if (length(acc_list) > 0) do.call(rbind, acc_list) else NULL

  if (!is.null(acc_df) && nrow(acc_df) > 0) {
    mean_acc_df <- acc_df %>%
      group_by(iteration) %>%
      summarize(mean_accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop")

    oracle_acc <- mean(sapply(pp_labeled_oracle[sem_keep_idx], function(y) {
      keep <- which(y$t > TREATMENT_TIME)
      if (length(keep) < 2) return(NA_real_)
      mean(y$inferred_process[keep] == y$process[keep])
    }), na.rm = TRUE)
    naive_acc <- mean(sapply(pp_labeled_naive[sem_keep_idx], function(y) {
      keep <- which(y$t > TREATMENT_TIME)
      if (length(keep) < 2) return(NA_real_)
      mean(y$inferred_process[keep] == y$process[keep])
    }), na.rm = TRUE)

    sim_study_plots$plot_em_accuracy_iters <- ggplot() +
    geom_line(data = acc_df, aes(x = iteration, y = accuracy, group = sim_id),
              color = "gray80", alpha = 0.5) +
    geom_line(data = mean_acc_df, aes(x = iteration, y = mean_accuracy),
              color = "red", linewidth = 1) +
    geom_hline(yintercept = oracle_acc, linetype = "dashed", color = "blue", linewidth = 0.7) +
    geom_hline(yintercept = naive_acc, linetype = "dotted", color = "orange", linewidth = 0.7) +
    annotate("text", x = max(acc_df$iteration), y = oracle_acc, label = "Oracle",
             hjust = 1, vjust = -0.5, color = "blue", size = 3) +
    annotate("text", x = max(acc_df$iteration), y = naive_acc, label = "Naive",
             hjust = 1, vjust = -0.5, color = "orange", size = 3) +
    labs(x = "Iteration", y = "Accuracy") +
    theme_minimal()
  }
}

sem_diagnostics_df_filtered <- sem_diagnostics_df
if (!is.null(sem_diagnostics_df) && nrow(sem_diagnostics_df) > 0) {
  sem_keep_idx <- method_kept_idx("SEM_adaptive")
  sem_diagnostics_df_filtered <- sem_diagnostics_df %>%
    filter(.data$sim_id %in% sem_keep_idx)
}

# Adaptive SEM objective/likelihood tracking over iterations
if (!is.null(sem_diagnostics_df_filtered) && nrow(sem_diagnostics_df_filtered) > 0 &&
    "metric" %in% names(sem_diagnostics_df_filtered)) {
  metric_df <- sem_diagnostics_df_filtered %>% filter(is.finite(.data$metric))
  if (!is.null(metric_df) && nrow(metric_df) > 0) {
    mean_metric_df <- metric_df %>%
      group_by(iteration) %>%
      summarize(mean_metric = mean(.data$metric, na.rm = TRUE), .groups = "drop")
    sim_study_plots$plot_em_loglik_iters <- ggplot() +
      geom_line(data = metric_df, aes(x = iteration, y = metric, group = sim_id),
                color = "gray80", alpha = 0.45) +
      geom_line(data = mean_metric_df, aes(x = iteration, y = mean_metric),
                color = "darkgreen", linewidth = 1) +
      labs(x = "Iteration", y = "Adaptive SEM objective (higher is better)") +
      theme_minimal()
  }
}

# Flips over iterations plot
if (!is.null(sem_diagnostics_df_filtered) && nrow(sem_diagnostics_df_filtered) > 0) {
  mean_flips_df <- sem_diagnostics_df_filtered %>%
    group_by(iteration) %>%
    summarize(
      mean_avg_flips = mean(average_flips, na.rm = TRUE),
      mean_max_flips = mean(max_metric_flips, na.rm = TRUE),
      .groups = "drop"
    )
  flips_long <- sem_diagnostics_df_filtered %>%
    reshape2::melt(id.vars = c("sim_id", "iteration"),
                   measure.vars = c("average_flips", "max_metric_flips"),
                   variable.name = "flip_type", value.name = "flips")

  sim_study_plots$plot_flips_iters <- ggplot() +
    geom_line(data = flips_long %>% filter(flip_type == "average_flips"),
              aes(x = iteration, y = flips, group = sim_id),
              color = "gray80", alpha = 0.4) +
    geom_line(data = mean_flips_df,
              aes(x = iteration, y = mean_avg_flips, color = "Mean avg flips"),
              linewidth = 1) +
    geom_line(data = mean_flips_df,
              aes(x = iteration, y = mean_max_flips, color = "Mean accepted flips"),
              linewidth = 1) +
    scale_color_manual(name = "", values = c("Mean avg flips" = "blue", "Mean accepted flips" = "red")) +
    labs(x = "Iteration", y = "Number of Flips") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

log_msg("Plots: ", paste(names(sim_study_plots), collapse = ", "))

# Publication-ready default: remove plot titles/subtitles.
sim_study_plots <- strip_titles_from_plots(sim_study_plots)

# ------------------------------------------------------------------
# 13. Timing report and save results (including plots)
# ------------------------------------------------------------------
time_end_global <- proc.time()[3]
elapsed_sec <- time_end_global - time_start_global
timing_report <- list(
  start_iso     = format(Sys.time() - elapsed_sec, "%Y-%m-%d %H:%M:%S %Z"),
  end_iso       = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  elapsed_sec   = elapsed_sec,
  elapsed_min   = round(elapsed_sec / 60, 2),
  run_id        = RUN_ID,
  job_id        = Sys.getenv("SLURM_JOB_ID", NA_character_),
  output_tag    = OUTPUT_TAG,
  n_cores       = N_CORES,
  sim_size      = SIM_SIZE,
  on_cluster    = ON_CLUSTER
)

slim_task_for_save <- function(task_obj) {
  if (is.null(task_obj)) return(NULL)
  list(
    labelling_name = task_obj$labelling_name,
    has_hawkes_params = !is.null(task_obj$hawkes_params),
    n_post = if (!is.null(task_obj$x)) nrow(task_obj$x) else NA_integer_,
    n_filtration = if (!is.null(task_obj$filtration_data)) nrow(task_obj$filtration_data) else NA_integer_
  )
}

slim_result_for_save <- function(res_obj) {
  if (is.null(res_obj)) return(NULL)
  list(
    all_nothing_theory = res_obj$all_nothing_theory,
    tau_1_estim = res_obj$tau_1_estim,
    ATE_total = res_obj$ATE_total,
    ATE_treatment = res_obj$ATE_treatment,
    ATE_spillover = res_obj$ATE_spillover,
    ATE_naive = res_obj$ATE_naive,
    control_pp = res_obj$control_pp,
    treated_pp = res_obj$treated_pp
  )
}

slim_em_result_for_save <- function(em_obj) {
  if (is.null(em_obj)) return(NULL)
  ad <- em_obj$adaptive
  list(
    hawkes_params_control = em_obj$hawkes_params_control,
    hawkes_params_treated = em_obj$hawkes_params_treated,
    time = em_obj$time,
    adaptive_summary = list(
      n_iter = if (!is.null(ad$accuracies)) length(ad$accuracies) else NA_integer_,
      final_accuracy = if (!is.null(ad$accuracies) && length(ad$accuracies) > 0) tail(ad$accuracies, 1) else NA_real_,
      final_metric = if (!is.null(ad$metrics) && length(ad$metrics) > 0) tail(ad$metrics, 1) else NA_real_,
      final_avg_flips = if (!is.null(ad$average_flips) && length(ad$average_flips) > 0) tail(ad$average_flips, 1) else NA_real_,
      final_acc_flips = if (!is.null(ad$max_metric_flips) && length(ad$max_metric_flips) > 0) tail(ad$max_metric_flips, 1) else NA_real_
    )
  )
}

results_flat_save <- if (isTRUE(SAVE_LIGHT)) {
  lapply(results_flat, slim_result_for_save)
} else {
  results_flat
}
tasks_save <- if (isTRUE(SAVE_LIGHT)) {
  lapply(tasks, slim_task_for_save)
} else {
  tasks
}
EM_results_save <- if (isTRUE(SAVE_LIGHT)) {
  lapply(EM_results, slim_em_result_for_save)
} else {
  EM_results
}
sem_diagnostics_all_save <- if (isTRUE(SAVE_LIGHT)) NULL else sem_diagnostics_df

sim_study_results <- list(
  results_df = results_df,
  summary_df = if (exists("summary_df")) summary_df else NULL,
  results_df_true_control = results_df_true_control,
  summary_df_true_control = summary_df_true_control,
  control_param_summary = if (exists("control_param_summary")) control_param_summary else NULL,
  treated_param_summary = if (exists("treated_param_summary")) treated_param_summary else NULL,
  control_params_df = if (exists("control_params_df")) control_params_df else NULL,
  treated_params_df = if (exists("treated_params_df")) treated_params_df else NULL,
  sem_pilot_summary = sem_pilot_summary,
  results_flat = results_flat_save,
  tasks = tasks_save,
  ate_run_idx = ate_run_idx,
  skipped_explosive_tasks = skipped_explosive_tasks,
  high_mu_failed_idx = high_mu_failed_idx,
  high_mu_failed_tasks = high_mu_failed_tasks,
  fit_status = fit_status_df,
  kept_result_idx = fit_status_df$task_idx[fit_status_df$include_result],
  class_metrics = class_metrics,
  label_recovery = label_recovery_all,
  label_recovery_prefilter = label_recovery_all_prefilter,
  support_contrast_specs = support_contrast_specs,
  support_contrast_df = support_contrast_df,
  support_contrast_summary = support_contrast_summary,
  decay_validation = decay_validation,
  all_nothing_ATE = all_nothing_ATE,
  true_tau_1 = true_tau_1,
  timing_report = timing_report,
  EM_results = EM_results_save,
  sem_diagnostics = sem_diagnostics_df_filtered,
  sem_diagnostics_all = sem_diagnostics_all_save,
  config = list(
    SIM_SIZE = SIM_SIZE, N_SIMS = N_SIMS, N_TAU_SIMS = N_TAU_SIMS, N_TAU_I = N_TAU_I,
    ATE_N_SIMS = ATE_N_SIMS, ATE_N_TAU_SIMS = ATE_N_TAU_SIMS, ATE_N_TAU_I = ATE_N_TAU_I, ATE_MAXIT = ATE_MAXIT,
    ATE_COMPUTE_TAU = ATE_COMPUTE_TAU,
    TRUE_CTRL_TAU_I = TRUE_CTRL_TAU_I, TRUE_CTRL_TAU_SIMS = TRUE_CTRL_TAU_SIMS,
    N_PROPOSALS = N_PROPOSALS,
    SEM_EM_ADAPTIVE_ITER = SEM_EM_ADAPTIVE_ITER,
    SEM_N_ITER = SEM_N_ITER, SEM_N_LABELLINGS = SEM_N_LABELLINGS,
    SEM_WORKERS = SEM_WORKERS,
    SEM_N_PROPS = SEM_N_PROPS,
    SEM_PARAM_UPDATE_CADENCE = SEM_PARAM_UPDATE_CADENCE,
    SEM_PROPOSAL_UPDATE_CADENCE = SEM_PROPOSAL_UPDATE_CADENCE,
    SEM_PARAM_REFIT_CADENCE = SEM_PARAM_REFIT_CADENCE,
    SEM_CHANGE_FACTOR = SEM_CHANGE_FACTOR,
    SIM_FILTRATION_AWARE = SIM_FILTRATION_AWARE,
    FILTER_HIGH_MU_FITS = FILTER_HIGH_MU_FITS,
    MU_FAIL_MULTIPLIER = MU_FAIL_MULTIPLIER,
    RUN_SEM_PILOT = RUN_SEM_PILOT,
    PILOT_ONLY = PILOT_ONLY,
    SEM_PILOT_ITERS = SEM_PILOT_ITERS,
    SEM_PILOT_SIMS = SEM_PILOT_SIMS,
    SEM_PILOT_CORES = SEM_PILOT_CORES,
    SEM_STALENESS_TRIGGER_EVERY = SEM_STALENESS_TRIGGER_EVERY,
    POST_TIME_MULTIPLIER = POST_TIME_MULTIPLIER,
    SCENARIO_ID = SCENARIO_ID,
    TREATMENT_ASSIGNMENT = TREATMENT_ASSIGNMENT,
    partition_assignment_counts = partition_assignment_counts,
    SIM_KERNEL = SIM_KERNEL,
    FIT_KERNEL = FIT_KERNEL,
    TARGET_POINTS = TARGET_POINTS,
    EXPECTED_POINTS_PER_MU = EXPECTED_POINTS_PER_MU,
    TRUE_MU = TRUE_MU,
    BASE_MU = BASE_MU,
    BASE_SEED = BASE_SEED,
    HAWKES_ALPHA = HAWKES_ALPHA,
    HAWKES_BETA = HAWKES_BETA,
    HAWKES_POWER_C = HAWKES_POWER_C,
    HAWKES_POWER_P = HAWKES_POWER_P,
    HAWKES_T_TRUNC = HAWKES_T_TRUNC,
    HAWKES_T_TRUNC_SOURCE = HAWKES_T_TRUNC_SOURCE,
    HAWKES_T_TRUNC_REL = HAWKES_T_TRUNC_REL,
    MU_SCALE = MU_SCALE,
    RUN_DECAY_VALIDATION = RUN_DECAY_VALIDATION,
    DECAY_VALIDATION_REPS = DECAY_VALIDATION_REPS,
    DECAY_ANNULUS_WIDTH = DECAY_ANNULUS_WIDTH,
    DECAY_FLIP_CELL = DECAY_FLIP_CELL,
    OUTPUT_TAG = OUTPUT_TAG,
    RUN_ID = RUN_ID,
    SAVE_LIGHT = SAVE_LIGHT,
    OMEGA = OMEGA, END_TIME = END_TIME, TREATMENT_TIME = TREATMENT_TIME,
    NX = NX, NY = NY, hawkes_par_1 = hawkes_par_1, hawkes_par_2 = hawkes_par_2,
    hawkes_fit_par_1 = hawkes_fit_par_1, hawkes_fit_par_2 = hawkes_fit_par_2
  ),
  plots = sim_study_plots
)

outfile <- file.path(SAVE_DIR, paste0(RUN_ID, ".rds"))
saveRDS(sim_study_results, outfile)
obj_mb <- as.numeric(object.size(sim_study_results)) / (1024^2)
log_msg(sprintf("Saved object size: %.1f MB", obj_mb))
log_msg("Results: ", outfile)
log_msg("Log:     ", LOG_FILE)
log_msg("=== DONE ", RUN_ID, " | ", round(elapsed_sec, 1), "s (", round(elapsed_sec / 60, 1), " min) ===")
close(log_con)
