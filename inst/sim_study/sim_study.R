#!/usr/bin/env Rscript
# Simulation study for PPDisentangle.
# Usage:
#   Rscript sim_study.R                  # local (auto-detect cores)
#   Rscript sim_study.R --test --sims 2  # quick smoke test
#   Rscript sim_study.R --small          # local, reduced
#   sbatch run_nesi.sh --sims 100        # NeSI cluster
#
# Output (canonical): output/sim_study/{JOB_ID}.rds  output/sim_study/{JOB_ID}.log

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
    canonical = file.path(repo_dir, "output", "sim_study"),
    legacy = file.path(script_dir, "output")
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
  SEM_N_ITER <- 100
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
LEGACY_SAVE_DIR <- save_dirs$legacy

if (!is.null(N_SIMS_ARG) && is.finite(N_SIMS_ARG)) {
  N_SIMS <- N_SIMS_ARG
  SIM_SIZE <- N_SIMS_ARG
}

if (ON_CLUSTER && !TEST) {
  slurm_cores <- suppressWarnings(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", NA_character_)))
  if (is.finite(slurm_cores) && slurm_cores >= 1) {
    if (!isTRUE(all.equal(SIM_SIZE, slurm_cores)) || !isTRUE(all.equal(N_CORES, slurm_cores))) {
      message(sprintf(
        "[sim_study] Aligning sims/cores to SLURM_CPUS_PER_TASK for efficiency: sims=%s, cores=%s, slurm=%s",
        as.character(SIM_SIZE), as.character(N_CORES), as.character(slurm_cores)
      ))
    }
    SIM_SIZE <- slurm_cores
    N_SIMS <- slurm_cores
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

# Optional runtime overrides for quick/full profiles without code edits.
N_PROPOSALS <- env_int("PP_LABEL_PROPOSALS", N_PROPOSALS, 1L)
SEM_EM_ADAPTIVE_ITER <- env_int("PP_SEM_INNER_ITER", SEM_EM_ADAPTIVE_ITER, 1L)
SEM_N_ITER <- env_int("PP_SEM_OUTER_ITER", SEM_N_ITER, 1L)
SEM_N_LABELLINGS <- env_int("PP_SEM_N_LABELLINGS", SEM_N_LABELLINGS, 1L)
POST_TIME_MULTIPLIER <- env_num("PP_POST_TIME_MULTIPLIER", 1, min_value = 0.1)
END_TIME <- TREATMENT_TIME + POST_TIME_MULTIPLIER * (END_TIME_BASE - TREATMENT_TIME)
if (!is.finite(END_TIME) || END_TIME <= TREATMENT_TIME) {
  stop("Computed END_TIME is invalid; check PP_POST_TIME_MULTIPLIER.")
}
TIME_INT <- END_TIME - TREATMENT_TIME

BASE_SEED <- 123L
stage_seed <- function(stage_offset, sim_id = 0L, extra = 0L) {
  as.integer(BASE_SEED + as.integer(stage_offset) * 100000L + as.integer(sim_id) * 1000L + as.integer(extra))
}

# SEM/adaptive configuration
SEM_PARAM_UPDATE_CADENCE <- 10L
SEM_PROPOSAL_UPDATE_CADENCE <- 1L
SEM_N_PROPS <- 20L
SEM_N_PROPS <- env_int("PP_SEM_N_PROPS", SEM_N_PROPS, 1L)
SEM_CHANGE_FACTOR <- env_num("PP_SEM_CHANGE_FACTOR", 0.005, min_value = 1e-6)
SEM_STALENESS_TRIGGER_EVERY <- 10L
SEM_STALENESS_TRIGGER_EVERY <- env_int("PP_SEM_STALENESS_TRIGGER_EVERY", SEM_STALENESS_TRIGGER_EVERY, 1L)
SEM_INCLUDE_STARTING <- tolower(Sys.getenv("PP_SEM_INCLUDE_STARTING", "true")) %in% c("1", "true", "yes", "y")
SEM_UPDATE_STARTING <- tolower(Sys.getenv("PP_SEM_UPDATE_STARTING", "true")) %in% c("1", "true", "yes", "y")
SEM_UPDATE_CONTROL_PARAMS <- tolower(Sys.getenv("PP_SEM_UPDATE_CONTROL_PARAMS", "false")) %in% c("1", "true", "yes", "y")

TREAT_PROP <- 0.5
MAX_TIME   <- 10000 * (END_TIME * OMEGA[2] * OMEGA[4] / 1e6)

# ------------------------------------------------------------------
# Logging
# ------------------------------------------------------------------
JOB_ID <- Sys.getenv("SLURM_JOB_ID", "")
if (JOB_ID == "") JOB_ID <- paste0("local_", format(Sys.time(), "%Y%m%d_%H%M%S"))

for (d in unique(c(SAVE_DIR, LEGACY_SAVE_DIR))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
save_dir_norm <- normalizePath(SAVE_DIR, winslash = "/", mustWork = FALSE)
legacy_save_dir_norm <- normalizePath(LEGACY_SAVE_DIR, winslash = "/", mustWork = FALSE)
mirror_to_legacy <- !identical(save_dir_norm, legacy_save_dir_norm)
mirror_to_legacy_file <- function(src_path) {
  if (!mirror_to_legacy) return(NA_character_)
  dest <- file.path(LEGACY_SAVE_DIR, basename(src_path))
  ok <- file.copy(src_path, dest, overwrite = TRUE)
  if (!isTRUE(ok)) {
    warning(sprintf("Failed to mirror %s to legacy path %s", basename(src_path), dest))
    return(NA_character_)
  }
  dest
}
LOG_FILE <- file.path(SAVE_DIR, paste0(JOB_ID, ".log"))
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
ATE_BATCH_SIZE <- max(ATE_WORKERS, 2L * ATE_WORKERS)
ATE_N_SIMS <- if (TEST) 1L else as.integer(N_SIMS)
ATE_N_TAU_SIMS <- if (TEST) 1L else as.integer(N_TAU_SIMS)
ATE_N_TAU_I <- if (TEST) 1L else as.integer(N_TAU_I)
ATE_MAXIT <- if (TEST) 300L else 1000L
ATE_N_SIMS <- env_int("PP_ATE_N_SIMS", ATE_N_SIMS, 1L)
ATE_N_TAU_SIMS <- env_int("PP_ATE_N_TAU_SIMS", ATE_N_TAU_SIMS, 1L)
ATE_N_TAU_I <- env_int("PP_ATE_N_TAU_I", ATE_N_TAU_I, 1L)
ATE_MAXIT <- env_int("PP_ATE_MAXIT", ATE_MAXIT, 1L)
SIM_FILTRATION_AWARE <- tolower(Sys.getenv("PP_SIM_FILTRATION_AWARE", "true")) %in% c("1", "true", "yes", "y")
ATE_CONTROL_FILTRATION_AWARE <- SIM_FILTRATION_AWARE &&
  (tolower(Sys.getenv("PP_ATE_CONTROL_FILTRATION_AWARE", "true")) %in% c("1", "true", "yes", "y"))
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
SEM_WORKERS_DEFAULT <- if (ON_CLUSTER && TEST) min(8L, as.integer(N_CORES)) else as.integer(N_CORES)
SEM_WORKERS <- env_int("PP_SEM_WORKERS", SEM_WORKERS_DEFAULT, 1L)
SEM_WORKERS <- min(as.integer(N_CORES), as.integer(SEM_WORKERS))
log_msg("=== ", JOB_ID, " | ", MODE, " | ", N_CORES, " cores x ", SIM_SIZE, " sims ===")
log_msg("SEM workers=", SEM_WORKERS)
log_msg("SEM adaptive inner=", SEM_EM_ADAPTIVE_ITER, " | outer=", SEM_N_ITER, " | labellings=", SEM_N_LABELLINGS)
log_msg("SEM spec: n_props=", SEM_N_PROPS,
        " | param_cadence=", SEM_PARAM_UPDATE_CADENCE,
        " | proposal_cadence=", SEM_PROPOSAL_UPDATE_CADENCE,
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
        " | maxit=", ATE_MAXIT,
        " | sim_filtration_aware=", SIM_FILTRATION_AWARE,
        " | control_filtration_aware=", ATE_CONTROL_FILTRATION_AWARE)
if (RUN_SEM_PILOT) {
  log_msg("SEM pilot: iters=", SEM_PILOT_ITERS,
          " | sims=", SEM_PILOT_SIMS,
          " | cores=", SEM_PILOT_CORES,
          " | pilot_only=", PILOT_ONLY)
}
log_msg("Base seed=", BASE_SEED)
log_msg("Post-treatment time multiplier=", POST_TIME_MULTIPLIER, " | END_TIME=", END_TIME)
log_msg("Output (canonical): ", SAVE_DIR)
if (mirror_to_legacy) {
  log_msg("Output (legacy mirror): ", LEGACY_SAVE_DIR)
}

# ------------------------------------------------------------------
# Helper: create a parallel cluster with PPDisentangle loaded
# ------------------------------------------------------------------
make_cluster <- function(n_cores) {
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
    "control_state_space", "hawkes_par_1", "hawkes_par_2",
    "N_PROPOSALS", "SEM_EM_ADAPTIVE_ITER", "SEM_N_ITER", "SEM_N_LABELLINGS",
    "SEM_PARAM_UPDATE_CADENCE", "SEM_PROPOSAL_UPDATE_CADENCE",
    "SEM_N_PROPS", "SEM_CHANGE_FACTOR", "SEM_STALENESS_TRIGGER_EVERY",
    "SEM_INCLUDE_STARTING", "SEM_UPDATE_STARTING", "SEM_UPDATE_CONTROL_PARAMS",
    "SIM_FILTRATION_AWARE"
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
hawkes_par_1 <- list(mu = 8, alpha = 0.01, beta = 10, K = 0.8)
hawkes_par_2 <- list(mu = 8, alpha = 0.01, beta = 10, K = 0.2)

partition <- quadrats(X = OMEGA, nx = NX, ny = NY)

all_nothing_ATE <- (hawkes_par_2$mu * TIME_INT * (1 / (1 - hawkes_par_2$K)) -
                      hawkes_par_1$mu * TIME_INT * (1 / (1 - hawkes_par_1$K))) / partition$n

partition_processes <- rep("control", partition$n)
set.seed(42)
partition_processes[sample(1:(NX * NY), NX * NY * TREAT_PROP)] <- "treated"
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
  cl <- make_cluster(SEM_WORKERS)
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
    process = "control", location_process = "control"
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
          hawkes_params_control = hawkes_par_1,
          hawkes_params_treated = NULL,
          change_factor = 1, filtration = pre, proximity_weight = 0
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
  total_points <- sum(dat$location_process == "treated" & dat$t >= TREATMENT_TIME)
  mu_start     <- total_points / TIME_INT
  params_init  <- list(mu = mu_start, alpha = 0.1, beta = TIME_INT / 10, K = 0.1)
  local_tuning <- if (is.null(tuning)) {
    list(
      param_update_cadence = SEM_PARAM_UPDATE_CADENCE,
      proposal_update_cadence = SEM_PROPOSAL_UPDATE_CADENCE,
      state_spaces = NULL,
      update_control_params = SEM_UPDATE_CONTROL_PARAMS,
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
    hawkes_params_control = hawkes_par_1,
    hawkes_params_treated = params_init,
    N_labellings = SEM_N_LABELLINGS,
    N_iter = SEM_N_ITER, verbose = FALSE,
    hawkes_use_filtration_history = SIM_FILTRATION_AWARE,
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
    cl_pilot <- make_cluster(SEM_PILOT_CORES)
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
log_msg("")
log_msg("=== Classification Accuracy ===")
class_metrics <- lapply(names(labellings), function(nm) {
  accs <- sapply(labellings[[nm]], function(y) {
    keep <- which(y$t > TREATMENT_TIME)
    if (length(keep) < 2) return(NA_real_)
    mean(y$inferred_process[keep] == y$process[keep])
  })
  data.frame(method = nm, mean_accuracy = mean(accs, na.rm = TRUE))
})
class_metrics <- do.call(rbind, class_metrics)
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
    return(list(mu = max(1e-8, nrow(ctrl_realiz) / dt_fit), alpha = 0, beta = dt_fit / 10, K = 0.01))
  }
  ctrl_init <- list(
    mu = max(1e-8, nrow(ctrl_realiz) / dt_fit),
    alpha = 0,
    beta = dt_fit / 100,
    K = 0.01
  )
  if (isTRUE(ATE_CONTROL_FILTRATION_AWARE)) {
    fit_ctrl <- fit_hawkes_with_filtration_fn(
      params_init = ctrl_init,
      realiz = ctrl_realiz,
      filtration = filt_data,
      windowT = c(TREATMENT_TIME, END_TIME),
      windowS = OMEGA,
      maxit = ATE_MAXIT,
      poisson_flag = FALSE,
      zero_background_region = treated_state_space
    )$par
    empirical_rate <- if (dt_fit > 0) nrow(ctrl_realiz) / dt_fit else Inf
    fitted_rate <- as.numeric(fit_ctrl$mu) / max(1e-6, 1 - as.numeric(fit_ctrl$K))
    degenerate_fit <- (!is.finite(fitted_rate) || !is.finite(empirical_rate) ||
                       (as.numeric(fit_ctrl$K) >= 0.98) ||
                       (fitted_rate < 0.2 * empirical_rate))
    if (isTRUE(degenerate_fit)) {
      legacy <- fit_hawkes_fn(
        unlist(ctrl_init),
        realiz = ctrl_realiz,
        zero_background_region = treated_state_space,
        windowT = c(TREATMENT_TIME, END_TIME),
        windowS = OMEGA,
        trace = 0,
        maxit = ATE_MAXIT,
        density_approx = FALSE,
        numeric_integral = FALSE,
        poisson_flag = FALSE,
        t_trunc = -1
      )$par
      fit_ctrl <- as.list(legacy)
      names(fit_ctrl) <- c("mu", "alpha", "beta", "K")
    }
  } else {
    fit_ctrl <- fit_hawkes_fn(
      unlist(ctrl_init),
      realiz = ctrl_realiz,
      zero_background_region = treated_state_space,
      windowT = c(TREATMENT_TIME, END_TIME),
      windowS = OMEGA,
      trace = 0,
      maxit = ATE_MAXIT,
      density_approx = FALSE,
      numeric_integral = FALSE,
      poisson_flag = FALSE,
      t_trunc = -1
    )$par
    fit_ctrl <- as.list(fit_ctrl)
    names(fit_ctrl) <- c("mu", "alpha", "beta", "K")
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
      labelling_name = nm, hawkes_params = NULL
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
ATE_env$ATE_CONTROL_FILTRATION_AWARE <- ATE_CONTROL_FILTRATION_AWARE
ATE_env$TREATMENT_TIME <- TREATMENT_TIME
ATE_env$END_TIME <- END_TIME
ATE_env$MAX_TIME <- MAX_TIME

task_function <- function(task) {
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
        treated_partitions = treated_partitions,
        hawkes_params = task$hawkes_params,
        n_sims = ATE_N_SIMS, n_tau_sims = ATE_N_TAU_SIMS, n_tau_i = ATE_N_TAU_I,
        windowT = c(TREATMENT_TIME, END_TIME), windowS = OMEGA,
        maxit = ATE_MAXIT,
        explosive_K_threshold = 0.98,
        poisson_flags = list(control = FALSE, treated = FALSE)
      ),
      timeout = MAX_TIME, onTimeout = "error"
    ),
    error = function(e) NULL
  )
  if (!is.null(r) && isTRUE(r$skipped_explosive)) return(NULL)
  if (!is.null(r)) r$all_nothing_sim <- NULL
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
    cl_ate <- make_cluster(min(ATE_WORKERS, length(ate_run_idx)))
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
results_df <- do.call(rbind, lapply(seq_along(results_flat), function(k) {
  r <- results_flat[[k]]
  if (is.null(r)) return(NULL)
  data.frame(
    labelling = tasks[[k]]$labelling_name,
    all_nothing_theory = mean(r$all_nothing_theory$ATE),
    tau_1_estim = r$tau_1_estim,
    ATE_total = r$ATE_total,
    ATE_treatment = r$ATE_treatment,
    ATE_spillover = r$ATE_spillover,
    ATE_naive = r$ATE_naive,
    stringsAsFactors = FALSE
  )
}))

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

TRUE_CTRL_TAU_I <- suppressWarnings(as.integer(Sys.getenv("PP_TRUECTRL_TAU_I", "3")))
if (!is.finite(TRUE_CTRL_TAU_I) || is.na(TRUE_CTRL_TAU_I) || TRUE_CTRL_TAU_I < 1L) TRUE_CTRL_TAU_I <- 3L
TRUE_CTRL_TAU_SIMS <- suppressWarnings(as.integer(Sys.getenv("PP_TRUECTRL_TAU_SIMS", "5")))
if (!is.finite(TRUE_CTRL_TAU_SIMS) || is.na(TRUE_CTRL_TAU_SIMS) || TRUE_CTRL_TAU_SIMS < 1L) TRUE_CTRL_TAU_SIMS <- 5L

results_df_true_control <- do.call(rbind, lapply(seq_along(results_flat), function(k) {
  r <- results_flat[[k]]
  if (is.null(r) || is.null(r$treated_pp)) return(NULL)
  treated_pp <- r$treated_pp
  ctrl_true <- hawkes_par_1
  tau_i_true <- vapply(seq_len(TRUE_CTRL_TAU_I), function(j) {
    tau_i(
      sample(partition$n, 1),
      partition = partition, treated_partitions = treated_partitions,
      statespace = OMEGA, windowT = c(TREATMENT_TIME, END_TIME),
      control_pp = ctrl_true, treated_pp = treated_pp, n_sim = TRUE_CTRL_TAU_SIMS
    )
  }, numeric(1))
  data.frame(
    labelling = tasks[[k]]$labelling_name,
    all_nothing_true_control =
      (treated_pp$mu * TIME_INT * (1 / (1 - treated_pp$K)) -
         ctrl_true$mu * TIME_INT * (1 / (1 - ctrl_true$K))) / partition$n,
    tau_1_true_control = mean(tau_i_true),
    ATE_total = r$ATE_total,
    ATE_treatment = r$ATE_treatment,
    ATE_spillover = r$ATE_spillover,
    ATE_naive = r$ATE_naive,
    stringsAsFactors = FALSE
  )
}))
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

log_msg("")
log_msg("True all/nothing ATE:", round(all_nothing_ATE, 4))
log_msg("True one-flip ATE:  ", round(true_tau_1, 4))

# ------------------------------------------------------------------
# 12. Create ggplots (as in old sim study) and save as elements of results
# ------------------------------------------------------------------
log_msg("")
log_msg("Building plots ...")
sim_study_plots <- list()
boxplot_method_levels <- c("oracle", "naive", "best", "SEM_adaptive", "SEM_full")

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
    data.frame(
      labelling = task_k$labelling_name,
      sim_id = ((k - 1) %% sim_size) + 1,
      mu = par_obj$mu, alpha = par_obj$alpha,
      beta = par_obj$beta, K = par_obj$K,
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

control_param_rows <- extract_param_rows(results_flat, tasks, SIM_SIZE, "control_pp")
treated_param_rows <- extract_param_rows(results_flat, tasks, SIM_SIZE, "treated_pp")
control_params_df <- do.call(rbind, control_param_rows)
treated_params_df <- do.call(rbind, treated_param_rows)
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
}

# Points per tile (control vs treated theoretical means)
ate_detail_rows <- lapply(seq_along(results_flat), function(k) {
  r <- results_flat[[k]]
  if (is.null(r) || is.null(r$all_nothing_theory)) return(NULL)
  th <- r$all_nothing_theory
  data.frame(
    labelling = rep(tasks[[k]]$labelling_name, 2),
    method = c("points_per_tile_control_theory", "points_per_tile_treated_theory"),
    ATE_estim = c(th$c_mean[1], th$t_mean[1]),
    stringsAsFactors = FALSE
  )
})
ate_detail_rows <- do.call(rbind, ate_detail_rows)
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
  acc_list <- lapply(seq_along(EM_results), function(i) {
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

    oracle_acc <- mean(sapply(pp_labeled_oracle, function(y) {
      keep <- which(y$t > TREATMENT_TIME)
      if (length(keep) < 2) return(NA_real_)
      mean(y$inferred_process[keep] == y$process[keep])
    }), na.rm = TRUE)
    naive_acc <- mean(sapply(pp_labeled_naive, function(y) {
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

# Adaptive SEM objective/likelihood tracking over iterations
if (!is.null(sem_diagnostics_df) && nrow(sem_diagnostics_df) > 0 &&
    "metric" %in% names(sem_diagnostics_df)) {
  metric_df <- sem_diagnostics_df %>% filter(is.finite(.data$metric))
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
if (!is.null(sem_diagnostics_df) && nrow(sem_diagnostics_df) > 0) {
  mean_flips_df <- sem_diagnostics_df %>%
    group_by(iteration) %>%
    summarize(
      mean_avg_flips = mean(average_flips, na.rm = TRUE),
      mean_max_flips = mean(max_metric_flips, na.rm = TRUE),
      .groups = "drop"
    )
  flips_long <- sem_diagnostics_df %>%
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
  job_id        = Sys.getenv("SLURM_JOB_ID", NA_character_),
  n_cores       = N_CORES,
  sim_size      = SIM_SIZE,
  on_cluster    = ON_CLUSTER
)

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
  results_flat = results_flat,
  tasks = tasks,
  ate_run_idx = ate_run_idx,
  skipped_explosive_tasks = skipped_explosive_tasks,
  class_metrics = class_metrics,
  all_nothing_ATE = all_nothing_ATE,
  true_tau_1 = true_tau_1,
  timing_report = timing_report,
  EM_results = EM_results,
  sem_diagnostics = sem_diagnostics_df,
  config = list(
    SIM_SIZE = SIM_SIZE, N_SIMS = N_SIMS, N_TAU_SIMS = N_TAU_SIMS, N_TAU_I = N_TAU_I,
    ATE_N_SIMS = ATE_N_SIMS, ATE_N_TAU_SIMS = ATE_N_TAU_SIMS, ATE_N_TAU_I = ATE_N_TAU_I, ATE_MAXIT = ATE_MAXIT,
    TRUE_CTRL_TAU_I = TRUE_CTRL_TAU_I, TRUE_CTRL_TAU_SIMS = TRUE_CTRL_TAU_SIMS,
    N_PROPOSALS = N_PROPOSALS,
    SEM_EM_ADAPTIVE_ITER = SEM_EM_ADAPTIVE_ITER,
    SEM_N_ITER = SEM_N_ITER, SEM_N_LABELLINGS = SEM_N_LABELLINGS,
    SEM_WORKERS = SEM_WORKERS,
    SEM_N_PROPS = SEM_N_PROPS,
    SEM_CHANGE_FACTOR = SEM_CHANGE_FACTOR,
    SIM_FILTRATION_AWARE = SIM_FILTRATION_AWARE,
    RUN_SEM_PILOT = RUN_SEM_PILOT,
    PILOT_ONLY = PILOT_ONLY,
    SEM_PILOT_ITERS = SEM_PILOT_ITERS,
    SEM_PILOT_SIMS = SEM_PILOT_SIMS,
    SEM_PILOT_CORES = SEM_PILOT_CORES,
    SEM_STALENESS_TRIGGER_EVERY = SEM_STALENESS_TRIGGER_EVERY,
    POST_TIME_MULTIPLIER = POST_TIME_MULTIPLIER,
    OMEGA = OMEGA, END_TIME = END_TIME, TREATMENT_TIME = TREATMENT_TIME,
    NX = NX, NY = NY, hawkes_par_1 = hawkes_par_1, hawkes_par_2 = hawkes_par_2
  ),
  plots = sim_study_plots
)

outfile <- file.path(SAVE_DIR, paste0(JOB_ID, ".rds"))
saveRDS(sim_study_results, outfile)
legacy_outfile <- mirror_to_legacy_file(outfile)
if (isTRUE(mirror_to_legacy)) {
  flush(log_con)
  legacy_log_file <- mirror_to_legacy_file(LOG_FILE)
}
log_msg("Results: ", outfile)
log_msg("Log:     ", LOG_FILE)
if (isTRUE(mirror_to_legacy)) {
  if (!is.na(legacy_outfile)) log_msg("Results (legacy mirror): ", legacy_outfile)
  if (!is.na(legacy_log_file)) log_msg("Log (legacy mirror):     ", legacy_log_file)
}
log_msg("=== DONE ", JOB_ID, " | ", round(elapsed_sec, 1), "s (", round(elapsed_sec / 60, 1), " min) ===")
close(log_con)
