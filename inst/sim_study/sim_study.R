#!/usr/bin/env Rscript
# Simulation study for PPDisentangle.
# Usage:
#   Rscript sim_study.R                  # local (auto-detect cores)
#   Rscript sim_study.R --test --sims 2  # quick smoke test
#   Rscript sim_study.R --small          # local, reduced
#   sbatch run_nesi.sh --sims 100        # NeSI cluster
#
# Output: inst/sim_study/output/{JOB_ID}.rds  inst/sim_study/output/{JOB_ID}.log

user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(user_lib)) {
  extra_libs <- strsplit(user_lib, .Platform$path.sep, fixed = TRUE)[[1]]
  extra_libs <- extra_libs[nzchar(extra_libs)]
  if (length(extra_libs) > 0L) {
    .libPaths(c(extra_libs, .libPaths()))
  }
}
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

resolve_save_dir <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  script_dir <- if (length(file_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    file.path(getwd(), "inst", "sim_study")
  }
  file.path(script_dir, "output")
}

OMEGA <- c(0, 100, 0, 100)
END_TIME <- 110
TREATMENT_TIME <- 10
NX <- 10
NY <- 10

local_core_default <- max(1L, parallel::detectCores() - 1L)
if (TEST) {
  n_test <- if (!is.null(N_SIMS_ARG) && is.finite(N_SIMS_ARG)) N_SIMS_ARG else 2L
  N_CORES <- if (ON_CLUSTER) as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", n_test)) else max(1L, min(n_test, parallel::detectCores()))
  SIM_SIZE <- n_test
  N_SIMS <- n_test
  N_TAU_SIMS <- 3
  N_TAU_I <- 3
  N_TAU_I_TRUE <- 5
  N_PROPOSALS <- 5
  SEM_EM_ADAPTIVE_ITER <- 100
  SEM_N_ITER <- 3
  SEM_N_LABELLINGS <- 5
  SAVE_DIR <- resolve_save_dir()
} else if (ON_CLUSTER) {
  N_CORES <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "100"))
  SIM_SIZE <- 100
  N_SIMS <- 100
  N_TAU_SIMS <- 10
  N_TAU_I <- 10
  N_TAU_I_TRUE <- 100
  N_PROPOSALS <- 10
  SEM_EM_ADAPTIVE_ITER <- 1000
  SEM_N_ITER <- if (LONG_MODE) 10 else 3
  SEM_N_LABELLINGS <- if (LONG_MODE) 20 else 10
  SAVE_DIR <- resolve_save_dir()
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
  SAVE_DIR <- resolve_save_dir()
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
  SAVE_DIR <- resolve_save_dir()
}

if (!is.null(N_SIMS_ARG) && is.finite(N_SIMS_ARG)) {
  N_SIMS <- N_SIMS_ARG
  SIM_SIZE <- N_SIMS_ARG
}

if (ON_CLUSTER) {
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

BASE_SEED <- 123L
stage_seed <- function(stage_offset, sim_id = 0L, extra = 0L) {
  as.integer(BASE_SEED + as.integer(stage_offset) * 100000L + as.integer(sim_id) * 1000L + as.integer(extra))
}

# SEM/adaptive configuration
SEM_PARAM_UPDATE_CADENCE <- 10L
SEM_PROPOSAL_UPDATE_CADENCE <- 1L
SEM_N_PROPS <- 10L
SEM_CHANGE_FACTOR <- 0.01
SEM_INCLUDE_STARTING <- TRUE
SEM_UPDATE_STARTING <- TRUE
SEM_UPDATE_CONTROL_PARAMS <- TRUE

TREAT_PROP <- 0.5
TIME_INT   <- END_TIME - TREATMENT_TIME
MAX_TIME   <- 10000 * (END_TIME * OMEGA[2] * OMEGA[4] / 1e6)

# ------------------------------------------------------------------
# Logging
# ------------------------------------------------------------------
JOB_ID <- Sys.getenv("SLURM_JOB_ID", "")
if (JOB_ID == "") JOB_ID <- paste0("local_", format(Sys.time(), "%Y%m%d_%H%M%S"))

dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)
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

params_are_crazy <- function(control_pp, treated_pp, K_max = 0.95, mu_max = 1e5) {
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
  min(12L, max(4L, as.integer(floor(N_CORES / 8))))
} else {
  N_CORES
}
ATE_BATCH_SIZE <- max(ATE_WORKERS, 2L * ATE_WORKERS)
log_msg("=== ", JOB_ID, " | ", MODE, " | ", N_CORES, " cores x ", SIM_SIZE, " sims ===")
log_msg("SEM adaptive inner=", SEM_EM_ADAPTIVE_ITER, " | outer=", SEM_N_ITER, " | labellings=", SEM_N_LABELLINGS)
log_msg("SEM spec: n_props=", SEM_N_PROPS,
        " | param_cadence=", SEM_PARAM_UPDATE_CADENCE,
        " | proposal_cadence=", SEM_PROPOSAL_UPDATE_CADENCE,
        " | change_factor=", SEM_CHANGE_FACTOR,
        " | include_starting=", SEM_INCLUDE_STARTING,
        " | update_starting=", SEM_UPDATE_STARTING,
        " | update_control_params=", SEM_UPDATE_CONTROL_PARAMS)
log_msg("ATE workers=", ATE_WORKERS)
log_msg("ATE batch size=", ATE_BATCH_SIZE)
log_msg("Base seed=", BASE_SEED)
log_msg("Output: ", SAVE_DIR)

# ------------------------------------------------------------------
# Helper: create a parallel cluster with PPDisentangle loaded
# ------------------------------------------------------------------
make_cluster <- function(n_cores) {
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, {
    user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
    if (nzchar(user_lib)) {
      extra_libs <- strsplit(user_lib, .Platform$path.sep, fixed = TRUE)[[1]]
      extra_libs <- extra_libs[nzchar(extra_libs)]
      if (length(extra_libs) > 0L) {
        .libPaths(c(extra_libs, .libPaths()))
      }
    }
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
    "SEM_N_PROPS", "SEM_CHANGE_FACTOR",
    "SEM_INCLUDE_STARTING", "SEM_UPDATE_STARTING", "SEM_UPDATE_CONTROL_PARAMS"
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

# Create cluster now that all globals are defined
cl <- make_cluster(N_CORES)
export_globals(cl)

# ------------------------------------------------------------------
# 2. Compute true average one-flip ATE
# ------------------------------------------------------------------
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

# ------------------------------------------------------------------
# 5. Labelling proposals
# ------------------------------------------------------------------
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
if (N_CORES > 1) {
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

# ------------------------------------------------------------------
# 7. Adaptive SEM
# ------------------------------------------------------------------
log_msg("Running adaptive SEM ...")
t0 <- proc.time()[3]
sem_jobs <- lapply(seq_along(obs_data), function(i) {
  list(i = i, seed = stage_seed(3L, i), data = obs_data[[i]])
})
run_sem_core <- function(job) {
  set.seed(job$seed)
  dat <- job$data
  total_points <- sum(dat$location_process == "treated" & dat$t >= TREATMENT_TIME)
  mu_start     <- total_points / TIME_INT
  params_init  <- list(mu = mu_start, alpha = 0.1, beta = TIME_INT / 10, K = 0.1)

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
    adaptive_control = list(
      param_update_cadence = SEM_PARAM_UPDATE_CADENCE,
      proposal_update_cadence = SEM_PROPOSAL_UPDATE_CADENCE,
      state_spaces = NULL,
      update_control_params = SEM_UPDATE_CONTROL_PARAMS,
      iter = SEM_EM_ADAPTIVE_ITER,
      n_props = SEM_N_PROPS,
      change_factor = SEM_CHANGE_FACTOR,
      include_starting_data = SEM_INCLUDE_STARTING,
      update_starting_data = SEM_UPDATE_STARTING,
      verbose = FALSE
    )
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
if (N_CORES > 1) {
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
  data.frame(
    sim_id = i,
    iteration = seq_len(n_iter),
    accuracy = a$accuracies,
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
    tasks[[length(tasks) + 1]] <- list(
      x = slim_for_ate(post_x),
      filtration_data = pre_histories[[i]],
      labelling_name = nm, hawkes_params = NULL
    )
  }
}
for (i in seq_along(EM_results)) {
  tmp <- obs_data[[i]]
  tmp$inferred_process <- tmp$location_process
  post_tmp <- tmp %>% filter(.data$t >= TREATMENT_TIME)
  tasks[[length(tasks) + 1]] <- list(
    x = slim_for_ate(post_tmp),
    filtration_data = pre_histories[[i]],
    labelling_name = "SEM_full",
    hawkes_params = list(
      control = EM_results[[i]]$hawkes_params_control,
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
  log_msg(sprintf("[CRAZY PARAMS] %d tasks with explosive params (K>=0.95 or mu>1e5)", length(crazy_idx)))
  if (SKIP_CRAZY) {
    tasks <- tasks[-crazy_idx]
    log_msg("[CRAZY PARAMS] Skipping explosive tasks")
  } else {
    log_msg("[CRAZY PARAMS] Retaining tasks in queue, but worker-side guard still returns NULL for explosive params")
  }
}

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
        treated_partitions = treated_partitions,
        hawkes_params = task$hawkes_params,
        n_sims = N_SIMS, n_tau_sims = N_TAU_SIMS, n_tau_i = N_TAU_I,
        windowT = c(TREATMENT_TIME, END_TIME), windowS = OMEGA,
        maxit = 1000,
        poisson_flags = list(control = FALSE, treated = FALSE)
      ),
      timeout = MAX_TIME, onTimeout = "error"
    ),
    error = function(e) NULL
  )
  if (!is.null(r)) r$all_nothing_sim <- NULL
  r
}
environment(task_function) <- ATE_env

if (N_CORES > 1) gc(verbose = FALSE)
ATE_SEQUENTIAL <- TEST
if (ATE_SEQUENTIAL) {
  log_msg("ATE estimation: sequential (TEST mode)")
  stopCluster(cl)
  results_flat <- lapply(tasks, task_function)
} else {
  # Recreate a fresh, smaller cluster for ATE to avoid stale worker state and
  # reduce memory pressure from too many concurrent heavy ATE jobs.
  stopCluster(cl)
  gc(verbose = FALSE)
  cl_ate <- make_cluster(min(ATE_WORKERS, length(tasks)))
  export_ate_globals(cl_ate)
  task_ids <- seq_along(tasks)
  task_batches <- split(task_ids, ceiling(task_ids / ATE_BATCH_SIZE))
  results_flat <- vector("list", length(tasks))
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
                    done_total, length(tasks), ok_total, null_total))
    rm(batch_res); gc(verbose = FALSE)
    log_elapsed("ATE parallel batch", proc.time()[3] - t0, max(idx), length(tasks))
  }
  stopCluster(cl_ate)
  if (batch_failed) {
    log_msg("[ATE PARALLEL ERROR] Falling back to sequential for robustness.")
    results_flat <- lapply(tasks, task_function)
  }
}
log_elapsed("ATE estimation", proc.time()[3] - t0, length(tasks), length(tasks))
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

log_msg("")
log_msg("True all/nothing ATE:", round(all_nothing_ATE, 4))
log_msg("True one-flip ATE:  ", round(true_tau_1, 4))

# ------------------------------------------------------------------
# 12. Create ggplots (as in old sim study) and save as elements of results
# ------------------------------------------------------------------
log_msg("")
log_msg("Building plots ...")
sim_study_plots <- list()

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
    par_obj <- r[[field_name]]
    if (is.null(r) || is.null(par_obj)) return(NULL)
    data.frame(
      labelling = tasks[[k]]$labelling_name,
      sim_id = ((k - 1) %% sim_size) + 1,
      mu = par_obj$mu, alpha = par_obj$alpha,
      beta = par_obj$beta, K = par_obj$K,
      stringsAsFactors = FALSE
    )
  })
}
build_param_boxplots <- function(params_df, truth_params) {
  if (is.null(params_df) || nrow(params_df) < 1) return(NULL)
  params_long <- reshape2::melt(
    params_df,
    id.vars = c("labelling", "sim_id"),
    variable.name = "param",
    value.name = "value"
  )
  plots <- lapply(c("mu", "alpha", "beta", "K"), function(p) {
    tmp <- params_long %>% filter(.data$param == p)
    if (nrow(tmp) == 0) return(NULL)
    oracle_vals <- tmp %>% filter(.data$labelling == "oracle")
    ggplot(tmp) +
      geom_boxplot(aes(x = .data$labelling, y = .data$value)) +
      geom_hline(yintercept = truth_params[[p]], linetype = "dashed", color = "red") +
      {if (nrow(oracle_vals) > 0) geom_hline(yintercept = mean(oracle_vals$value, na.rm = TRUE), linetype = "dotted", color = "blue") } +
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

control_param_plots <- build_param_boxplots(control_params_df, hawkes_par_1)
if (!is.null(control_param_plots) && length(control_param_plots) > 0) {
  sim_study_plots$plot_control_params <- control_param_plots
}

treated_param_plots <- build_param_boxplots(treated_params_df, hawkes_par_2)
if (!is.null(treated_param_plots) && length(treated_param_plots) > 0) {
  sim_study_plots$plot_treated_params <- treated_param_plots
}

# All-nothing ATE boxplot
if (!is.null(results_df) && nrow(results_df) > 0) {
  results_df$labelling <- factor(results_df$labelling,
    levels = unique(c("oracle", "naive", "best", "SEM_adaptive", "SEM_full", results_df$labelling)))
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
  results_flat = results_flat,
  tasks = tasks,
  class_metrics = class_metrics,
  all_nothing_ATE = all_nothing_ATE,
  true_tau_1 = true_tau_1,
  timing_report = timing_report,
  EM_results = EM_results,
  sem_diagnostics = sem_diagnostics_df,
  config = list(
    SIM_SIZE = SIM_SIZE, N_SIMS = N_SIMS, N_TAU_SIMS = N_TAU_SIMS, N_TAU_I = N_TAU_I,
    N_PROPOSALS = N_PROPOSALS,
    SEM_EM_ADAPTIVE_ITER = SEM_EM_ADAPTIVE_ITER,
    SEM_N_ITER = SEM_N_ITER, SEM_N_LABELLINGS = SEM_N_LABELLINGS,
    OMEGA = OMEGA, END_TIME = END_TIME, TREATMENT_TIME = TREATMENT_TIME,
    NX = NX, NY = NY, hawkes_par_1 = hawkes_par_1, hawkes_par_2 = hawkes_par_2
  ),
  plots = sim_study_plots
)

outfile <- file.path(SAVE_DIR, paste0(JOB_ID, ".rds"))
saveRDS(sim_study_results, outfile)
log_msg("Results: ", outfile)
log_msg("Log:     ", LOG_FILE)
log_msg("=== DONE ", JOB_ID, " | ", round(elapsed_sec, 1), "s (", round(elapsed_sec / 60, 1), " min) ===")
close(log_con)
