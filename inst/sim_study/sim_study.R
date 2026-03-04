#!/usr/bin/env Rscript
# Unified simulation study for the PPDisentangle package.
# Detects SLURM environment automatically; falls back to local defaults.
# Usage:
#   Local:   Rscript sim_study.R [--small]
#   Cluster (batch): sbatch run_PPDisentangle_sim.slurm
#   Cluster (interactive): Rscript sim_study.R --cluster   (or run_sim_study.sh --cluster)
#   Quick-test: Rscript sim_study.R --test [--sims N]   (10 iters, 2 sims/cores default; verifies machinery)

library(PPDisentangle)
library(spatstat)
library(ggplot2)
library(dplyr)
library(data.table)
library(parallel)
library(doParallel)
library(R.utils)
library(reshape2)
library(gridExtra)
library(scales)

time_start_global <- proc.time()[3]

args <- commandArgs(trailingOnly = TRUE)
SMALL <- "--small" %in% args
FORCE_CLUSTER <- "--cluster" %in% args
TEST <- "--test" %in% args

# Parse --sims N from command line (e.g. run_sim_study.sh --sims 100)
sims_arg <- grep("^--sims$", args)
N_SIMS_ARG <- if (length(sims_arg) > 0 && length(args) >= sims_arg + 1L)
  as.numeric(args[sims_arg + 1L]) else NULL

# Cluster config: when running under SLURM or when --cluster is passed (e.g. interactive on login node)
ON_CLUSTER <- nzchar(Sys.getenv("SLURM_JOB_ID")) || FORCE_CLUSTER

if (TEST) {
  # Quick-test mode: 10 iters for everything, minimal sims to verify machinery (logging, parallel, etc.)
  n_test <- if (!is.null(N_SIMS_ARG) && is.finite(N_SIMS_ARG)) N_SIMS_ARG else 2L
  N_CORES      <- if (ON_CLUSTER) as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", n_test)) else max(1L, min(n_test, parallel::detectCores()))
  SIM_SIZE     <- n_test
  N_SIMS       <- n_test
  N_TAU_SIMS   <- 10
  N_TAU_I      <- 10
  N_TAU_I_TRUE <- 10
  N_PROPOSALS  <- 10
  EM_ITER      <- 10
  SEM_EM_ADAPTIVE_ITER <- 10
  SEM_N_ITER   <- 10
  SEM_N_LABELLINGS <- 10
  OMEGA        <- c(0, 100, 0, 100)
  END_TIME     <- 110
  TREATMENT_TIME <- 10
  NX <- 10; NY <- 10
  SAVE_DIR     <- file.path(if (nzchar(Sys.getenv("SLURM_JOB_ID"))) Sys.getenv("SLURM_SUBMIT_DIR", getwd()) else getwd(), "cluster_output")
} else if (ON_CLUSTER) {
  N_CORES      <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
  SIM_SIZE     <- N_CORES
  N_SIMS       <- 100
  N_TAU_SIMS   <- 10
  N_TAU_I      <- 10
  N_TAU_I_TRUE <- 100
  N_PROPOSALS  <- 10
  EM_ITER      <- 100
  SEM_EM_ADAPTIVE_ITER <- 100
  SEM_N_ITER   <- 10
  SEM_N_LABELLINGS <- 100
  OMEGA        <- c(0, 100, 0, 100)
  END_TIME     <- 110
  TREATMENT_TIME <- 10
  NX <- 10; NY <- 10
  SAVE_DIR     <- file.path(Sys.getenv("SLURM_SUBMIT_DIR", getwd()), "cluster_output")
} else if (SMALL) {
  N_CORES      <- max(1, parallel::detectCores() - 1)
  SIM_SIZE     <- N_CORES
  N_SIMS       <- 10
  N_TAU_SIMS   <- 10
  N_TAU_I      <- 10
  N_TAU_I_TRUE <- 10
  N_PROPOSALS  <- 5
  EM_ITER      <- 10
  SEM_EM_ADAPTIVE_ITER <- 10
  SEM_N_ITER   <- 10
  SEM_N_LABELLINGS <- 10
  OMEGA        <- c(0, 100, 0, 100)
  END_TIME     <- 110
  TREATMENT_TIME <- 10
  NX <- 10; NY <- 10
  SAVE_DIR     <- file.path(getwd(), "cluster_output")
} else {
  N_CORES      <- max(1, parallel::detectCores() - 1)
  SIM_SIZE     <- N_CORES
  N_SIMS       <- 10
  N_TAU_SIMS   <- 10
  N_TAU_I      <- 10
  N_TAU_I_TRUE <- 10
  N_PROPOSALS  <- 100
  EM_ITER      <- 100
  SEM_EM_ADAPTIVE_ITER <- 100
  SEM_N_ITER   <- 10
  SEM_N_LABELLINGS <- max(10, N_PROPOSALS %/% 10)
  OMEGA        <- c(0, 100, 0, 100)
  END_TIME     <- 110
  TREATMENT_TIME <- 10
  NX <- 10; NY <- 10
  SAVE_DIR     <- file.path(getwd(), "cluster_output")
}

# Environment and CLI overrides (for cluster or shell script)
if (!is.null(N_SIMS_ARG) && is.finite(N_SIMS_ARG)) {
  N_SIMS <- N_SIMS_ARG
  SIM_SIZE <- N_SIMS_ARG
}
if (nzchar(Sys.getenv("CORES_OVERRIDE"))) {
  N_CORES <- as.numeric(Sys.getenv("CORES_OVERRIDE"))
}
if (nzchar(Sys.getenv("N_SIMS_OVERRIDE"))) {
  n_ov <- as.numeric(Sys.getenv("N_SIMS_OVERRIDE"))
  N_SIMS <- n_ov
  SIM_SIZE <- n_ov
}
if (nzchar(Sys.getenv("SIM_SIZE_OVERRIDE"))) {
  SIM_SIZE <- as.numeric(Sys.getenv("SIM_SIZE_OVERRIDE"))
}

TREAT_PROP <- 0.5
TIME_INT   <- END_TIME - TREATMENT_TIME
MAX_TIME   <- 10000 * (END_TIME * OMEGA[2] * OMEGA[4] / 1e6)

# ------------------------------------------------------------------
# Logging: timestamped messages to console and log file
# ------------------------------------------------------------------
# Use RUN_ID from environment (set by shell script) or generate new one
RUN_ID <- Sys.getenv("PP_RUN_ID", "")
if (RUN_ID == "") RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")

LOG_DIR <- file.path(SAVE_DIR, "logs", RUN_ID)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)
LOG_FILE <- file.path(LOG_DIR, "sim_study.log")
log_con <- file(LOG_FILE, open = "wt")
on.exit(tryCatch(close(log_con), error = function(e) NULL), add = TRUE)
log_msg <- function(...) {
  msg <- paste0(...)
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s\n", ts, msg)
  if (!ON_CLUSTER) cat(line)
  if (isOpen(log_con)) cat(line, file = log_con)
}
log_file_only <- function(...) {
  msg <- paste0(...)
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s\n", ts, msg)
  if (isOpen(log_con)) cat(line, file = log_con)
}
log_elapsed <- function(phase, elapsed_sec, n_done = NULL, n_total = NULL) {
  if (is.null(n_done) || is.null(n_total) || n_total <= 0) {
    log_msg(sprintf("%s: %.1f s", phase, elapsed_sec))
  } else {
    eta_sec <- if (n_done > 0) (elapsed_sec / n_done) * (n_total - n_done) else NA
    log_msg(sprintf("%s: %.1f s (done %d/%d, ETA %.0f s)", phase, elapsed_sec, n_done, n_total, eta_sec))
  }
}

log_msg("=== Simulation Study Config ===")
log_msg("Mode:", if (TEST) "TEST" else if (ON_CLUSTER) "CLUSTER" else if (SMALL) "SMALL/LOCAL" else "LOCAL")
log_msg("Cores:", N_CORES, " Sims:", SIM_SIZE)
log_msg("Omega:", paste(OMEGA, collapse = " "), " T:", END_TIME, " t*:", TREATMENT_TIME)
log_msg("Grid:", NX, "x", NY)
log_msg("Save:", SAVE_DIR)
log_msg("Log file:", LOG_FILE)
log_msg("EM-style: iter=", EM_ITER, " n_props=10 | SEM adaptive: iter=", SEM_EM_ADAPTIVE_ITER, " n_props=10 | SEM full: N_iter=", SEM_N_ITER, " N_labellings=", SEM_N_LABELLINGS)
log_msg("===============================")

# ------------------------------------------------------------------
# Helper: create a parallel cluster with PPDisentangle loaded
# ------------------------------------------------------------------
make_cluster <- function(n_cores) {
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, {
    library(PPDisentangle)
    library(R.utils)
  })
  return(cl)
}

export_globals <- function(cl) {
  clusterExport(cl, c(
    "OMEGA", "partition", "N_SIMS", "MAX_TIME",
    "N_TAU_I", "N_TAU_SIMS", "N_TAU_I_TRUE",
    "TREATMENT_TIME", "END_TIME",
    "TIME_INT", "state_spaces", "partition_processes",
    "treated_partitions", "treated_state_space",
    "control_state_space", "hawkes_par_1", "hawkes_par_2",
    "N_PROPOSALS", "EM_ITER", "SEM_EM_ADAPTIVE_ITER", "SEM_N_ITER", "SEM_N_LABELLINGS"
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
set.seed(123)
obs_data <- lapply(1:SIM_SIZE, function(i) {
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
gen_proposals <- function(s) {
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
if (N_CORES > 1 && !SMALL) {
  labelling_proposals <- parLapply(cl = cl, X = obs_data, fun = gen_proposals)
} else {
  labelling_proposals <- lapply(obs_data, gen_proposals)
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
# 6. EM-style labelling
# ------------------------------------------------------------------
log_msg("Running EM-style labelling ...")
t0 <- proc.time()[3]
run_em <- function(x) {
  total_points <- sum(x$location_process == "treated" & x$t >= TREATMENT_TIME)
  mu_start     <- total_points / TIME_INT
  params_init  <- list(mu = mu_start, alpha = 0.1, beta = TIME_INT / 10, K = 0.1)

  em_style_labelling(
    pp_data = x, partition = partition,
    partition_processes = partition_processes,
    statespace = OMEGA, state_spaces = state_spaces,
    time_window = c(TREATMENT_TIME, END_TIME),
    treatment_time = TREATMENT_TIME,
    hawkes_params_control = hawkes_par_1,
    hawkes_params_treated = params_init,
    param_update_cadence = 10, proposal_update_cadence = 1,
    update_starting_data = TRUE, include_starting_data = FALSE,
    metric_name = "post_likelihood", optim_method = "max",
    iter = EM_ITER, n_props = 10, change_factor = 0.02,
    MCMC_style = FALSE, verbose = FALSE
  )
}
if (N_CORES > 1 && !SMALL) {
  EM_max_style <- parLapply(cl = cl, X = obs_data, fun = run_em)
} else {
  EM_max_style <- lapply(obs_data, run_em)
}
log_elapsed("EM-style labelling", proc.time()[3] - t0, SIM_SIZE, SIM_SIZE)

# Log EM-style summary per sim: accuracy, n_iters, timing (log file only)
log_file_only("EM-style per-sim: accuracy | n_iter | sampling_s | likelihood_s | param_update_s")
for (k in seq_along(EM_max_style)) {
  r <- EM_max_style[[k]]
  acc <- if (length(r$accuracies) > 0) r$accuracies[length(r$accuracies)] else NA_real_
  t <- r$timing
  log_file_only(sprintf("  sim %d: %.3f | %d | %.1f | %.1f | %.1f",
    k, acc, t$n_iter, t$sampling_s, t$likelihood_s, t$param_update_s))
}

pp_labelled_em_post <- lapply(EM_max_style, function(x) x$labelling)

# ------------------------------------------------------------------
# 7. Adaptive SEM
# ------------------------------------------------------------------
log_msg("Running adaptive SEM ...")
t0 <- proc.time()[3]
run_sem <- function(dat) {
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
      param_update_cadence = 10,
      proposal_update_cadence = 1,
      state_spaces = NULL,
      update_control_params = FALSE,
      iter = SEM_EM_ADAPTIVE_ITER,
      n_props = 10,
      change_factor = 0.02,
      include_starting_data = FALSE,
      update_starting_data = TRUE,
      verbose = FALSE
    )
  )
}
if (N_CORES > 1 && !SMALL) {
  EM_results <- parLapply(cl = cl, X = obs_data, fun = run_sem)
} else {
  EM_results <- lapply(obs_data, run_sem)
}
log_elapsed("Adaptive SEM", proc.time()[3] - t0, SIM_SIZE, SIM_SIZE)

# ------------------------------------------------------------------
# 8. Assemble all labellings
# ------------------------------------------------------------------
labellings <- list(
  oracle   = pp_labeled_oracle,
  naive    = pp_labeled_naive,
  best     = pp_labeled_best_proposal,
  EM_style = pp_labelled_em_post,
  EM_full  = lapply(EM_results, function(x) x$adaptive$adaptive_labelling)
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
cm_out <- capture.output(print(class_metrics))
for (line in cm_out) log_file_only(line)
if (!ON_CLUSTER) print(class_metrics)

# ------------------------------------------------------------------
# 10. ATE estimation
# ------------------------------------------------------------------
log_msg("Estimating ATEs ...")
t0 <- proc.time()[3]

tasks <- list()
for (nm in names(labellings)) {
  for (i in seq_along(labellings[[nm]])) {
    tasks[[length(tasks) + 1]] <- list(
      x = labellings[[nm]][[i]] %>% filter(.data$t >= TREATMENT_TIME),
      labelling_name = nm, hawkes_params = NULL
    )
  }
}
for (i in seq_along(EM_results)) {
  tmp <- obs_data[[i]]
  tmp$inferred_process <- tmp$location_process
  tasks[[length(tasks) + 1]] <- list(
    x = tmp %>% filter(.data$t >= TREATMENT_TIME),
    labelling_name = "EM_full_params",
    hawkes_params = list(
      control = EM_results[[i]]$hawkes_params_control,
      treated = EM_results[[i]]$hawkes_params_treated
    )
  )
}

# Use minimal env for task_function to avoid serializing .GlobalEnv (obs_data,
# EM_results, labellings, etc.) when sending to workers - that causes OOM.
ATE_env <- new.env(parent = baseenv())
ATE_env$ATE_estim_hawkes <- getFromNamespace("ATE_estim_hawkes", "PPDisentangle")
ATE_env$withTimeout <- getFromNamespace("withTimeout", "R.utils")
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
  r <- tryCatch(
    withTimeout(
      ATE_estim_hawkes(
        statespace = OMEGA, partition = partition,
        observed_data = task$x,
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
  if (!is.null(r)) r$all_nothing_sim <- NULL  # drop sims, keep theory (used in plots)
  r
}
environment(task_function) <- ATE_env

# Free memory before ATE phase; workers get fresh tasks
if (N_CORES > 1 && !SMALL) gc(verbose = FALSE)

# Sequential ATE fallback: set env ATE_SEQUENTIAL=1 to avoid worker OOM (slower but lower memory)
ATE_SEQUENTIAL <- nzchar(Sys.getenv("ATE_SEQUENTIAL")) && tolower(Sys.getenv("ATE_SEQUENTIAL")) %in% c("1", "true", "yes")
if (ATE_SEQUENTIAL) {
  log_msg("ATE estimation: sequential (ATE_SEQUENTIAL=1)")
  stopCluster(cl)
  results_flat <- lapply(tasks, task_function)
} else {
  results_flat <- parLapply(cl, tasks, fun = task_function)
  stopCluster(cl)
}
log_elapsed("ATE estimation", proc.time()[3] - t0, length(tasks), length(tasks))

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
      title = "Simulated Hawkes-Hawkes (post-treatment)")
  }
  if (nrow(pre_first) > 0) {
    sim_study_plots$plot_pp_pre <- plot_pp(pre_first, partition = partition,
      title = "Simulated Hawkes-Hawkes (pre-treatment)")
  }
}

# Control and treated parameter estimates from results_flat
control_param_rows <- lapply(seq_along(results_flat), function(k) {
  r <- results_flat[[k]]
  if (is.null(r) || is.null(r$control_pp)) return(NULL)
  data.frame(
    labelling = tasks[[k]]$labelling_name,
    sim_id = ((k - 1) %% SIM_SIZE) + 1,
    mu = r$control_pp$mu, alpha = r$control_pp$alpha,
    beta = r$control_pp$beta, K = r$control_pp$K,
    stringsAsFactors = FALSE
  )
})
treated_param_rows <- lapply(seq_along(results_flat), function(k) {
  r <- results_flat[[k]]
  if (is.null(r) || is.null(r$treated_pp)) return(NULL)
  data.frame(
    labelling = tasks[[k]]$labelling_name,
    sim_id = ((k - 1) %% SIM_SIZE) + 1,
    mu = r$treated_pp$mu, alpha = r$treated_pp$alpha,
    beta = r$treated_pp$beta, K = r$treated_pp$K,
    stringsAsFactors = FALSE
  )
})
control_params_df <- do.call(rbind, control_param_rows)
treated_params_df <- do.call(rbind, treated_param_rows)

if (!is.null(control_params_df) && nrow(control_params_df) > 0) {
  control_params_long <- reshape2::melt(control_params_df,
    id.vars = c("labelling", "sim_id"),
    variable.name = "param", value.name = "value")
  control_param_plots <- lapply(c("mu", "alpha", "beta", "K"), function(p) {
    tmp <- control_params_long %>% filter(.data$param == p)
    if (nrow(tmp) == 0) return(NULL)
    ggplot(tmp) +
      geom_boxplot(aes(x = .data$labelling, y = .data$value)) +
      geom_hline(yintercept = hawkes_par_1[[p]], linetype = "dashed", color = "red") +
      labs(title = paste0(p, " estimates for Hawkes-Hawkes"), x = "Method", y = "Parameter Estimate") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
  })
  control_param_plots <- Filter(Negate(is.null), control_param_plots)
  if (length(control_param_plots) > 0) {
    sim_study_plots$plot_control_params <- control_param_plots
  }
}

if (!is.null(treated_params_df) && nrow(treated_params_df) > 0) {
  treated_params_long <- reshape2::melt(treated_params_df,
    id.vars = c("labelling", "sim_id"),
    variable.name = "param", value.name = "value")
  treated_param_plots <- lapply(c("mu", "alpha", "beta", "K"), function(p) {
    tmp <- treated_params_long %>% filter(.data$param == p)
    if (nrow(tmp) == 0) return(NULL)
    ggplot(tmp) +
      geom_boxplot(aes(x = .data$labelling, y = .data$value)) +
      geom_hline(yintercept = hawkes_par_2[[p]], linetype = "dashed", color = "red") +
      labs(title = paste0(p, " estimates for Hawkes-Hawkes example"), x = "Method", y = "Parameter Estimate") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
  })
  treated_param_plots <- Filter(Negate(is.null), treated_param_plots)
  if (length(treated_param_plots) > 0) {
    sim_study_plots$plot_treated_params <- treated_param_plots
  }
}

# All-nothing ATE boxplot (match old: title, y label, linetype legend, ylim)
if (!is.null(results_df) && nrow(results_df) > 0) {
  results_df$labelling <- factor(results_df$labelling,
    levels = unique(c("oracle", "naive", "best", "EM_style", "EM_full", "EM_full_params", results_df$labelling)))
  lines_data_ate <- data.frame(all_nothing_ATE = all_nothing_ATE)
  y_lo <- min(all_nothing_ATE * 1.5, min(results_df$all_nothing_theory, na.rm = TRUE) * 0.95)
  y_hi <- max(3, max(results_df$all_nothing_theory, na.rm = TRUE) * 1.05)
  sim_study_plots$plot_all_nothing_ATE <- ggplot(results_df) +
    geom_boxplot(aes(x = .data$labelling, y = .data$all_nothing_theory)) +
    geom_hline(data = lines_data_ate, aes(yintercept = .data$all_nothing_ATE, linetype = "True All/Nothing ATE"),
               color = scales::hue_pal()(3)[1], linewidth = 1) +
    scale_linetype_manual(name = "", values = c("True All/Nothing ATE" = "solid", "Average Single Flip ATE" = "solid")) +
    labs(title = "All-Nothing ATE estimates for Hawkes-Hawkes example", x = "Method", y = "All-Nothing ATE Estimate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
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
    summarize(mean_ATE = median(.data$ATE_estim), .groups = "drop")
  sim_study_plots$plot_points_per_tile <- ggplot(ate_detail_rows, aes(x = .data$labelling, y = .data$ATE_estim, fill = .data$method)) +
    geom_boxplot() +
    facet_wrap(~ method) +
    geom_hline(aes(yintercept = .data$mean_ATE, colour = "Oracle mean"), data = oracle_means_pts, linetype = "dashed", linewidth = 1) +
    geom_hline(aes(yintercept = .data$mean_ATE, colour = "True mean"), data = true_means_pts, linetype = "dashed", linewidth = 1) +
    scale_colour_manual(name = "Reference", values = c("Oracle mean" = "#0072B2", "True mean" = "#D55E00")) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
    labs(title = "Points-Per-Tile Hawkes-Hawkes example", x = "Method", y = "Points Per Tile") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
    coord_cartesian(ylim = c(0, max(true_means_pts$mean_ATE) * 1.05))
}

# One-flip ATE boxplot (match old: title, y label, linetype legend)
if (!is.null(results_df) && nrow(results_df) > 0) {
  lines_data_tau <- data.frame(true_1_flip = true_tau_1)
  sim_study_plots$plot_one_flip_ATE <- ggplot(results_df) +
    geom_boxplot(aes(x = .data$labelling, y = .data$tau_1_estim)) +
    geom_hline(data = lines_data_tau, aes(yintercept = .data$true_1_flip, linetype = "Average Single Flip ATE"),
               color = scales::hue_pal()(3)[1], linewidth = 1) +
    scale_linetype_manual(name = "", values = c("All/Nothing ATE" = "solid", "Average Single Flip ATE" = "solid")) +
    labs(title = "Average Single Flip ATE for Hawkes-Hawkes example", x = "Method", y = "Single Flip ATE Estimate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
}

# Combined param plots: grid.arrange() returns a grid object that does not print after save/load.
# Store draw functions so you can display combined plots after loading the RDS:
sim_study_plots$draw_control_params_combined <- function() {
  if (length(sim_study_plots$plot_control_params) > 0)
    gridExtra::grid.arrange(grobs = sim_study_plots$plot_control_params, ncol = 2)
}
sim_study_plots$draw_treated_params_combined <- function() {
  if (length(sim_study_plots$plot_treated_params) > 0)
    gridExtra::grid.arrange(grobs = sim_study_plots$plot_treated_params, ncol = 2)
}

# EM Accuracy over iterations plot
if (exists("EM_results") && length(EM_results) > 0) {
  # Extract accuracies from each simulation run
  # EM_results is a list of results from run_sem (adaptive_SEM)
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
    # Calculate mean accuracy per iteration
    mean_acc_df <- acc_df %>%
      group_by(iteration) %>%
      summarize(mean_accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop")

    sim_study_plots$plot_em_accuracy_iters <- ggplot() +
    # Individual runs in gray
    geom_line(data = acc_df, aes(x = iteration, y = accuracy, group = sim_id), 
              color = "gray80", alpha = 0.5) +
    # Mean in red
    geom_line(data = mean_acc_df, aes(x = iteration, y = mean_accuracy), 
              color = "red", linewidth = 1) +
    labs(title = "EM Accuracy over Iterations",
         subtitle = "Gray lines: individual runs; Red line: mean accuracy",
         x = "Iteration", y = "Accuracy") +
    theme_minimal()
  }
}

log_msg("Plots done.")
log_msg("Plots stored in sim_study_plots:", paste(names(sim_study_plots), collapse = ", "))
log_msg("  To print combined param plots after load: x <- readRDS(...); x$plots$draw_control_params_combined(); x$plots$draw_treated_params_combined()")

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
  config = list(
    SIM_SIZE = SIM_SIZE, N_SIMS = N_SIMS, N_TAU_SIMS = N_TAU_SIMS, N_TAU_I = N_TAU_I,
    OMEGA = OMEGA, END_TIME = END_TIME, TREATMENT_TIME = TREATMENT_TIME,
    NX = NX, NY = NY, hawkes_par_1 = hawkes_par_1, hawkes_par_2 = hawkes_par_2
  ),
  plots = sim_study_plots
)

dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(LOG_DIR, "results.rds")
saveRDS(sim_study_results, outfile)
log_msg("Results saved to:", outfile)

# Write timing report to a text file (for cluster logs and quick inspection)
timing_txt <- file.path(LOG_DIR, "timing.txt")
writeLines(c(
  "PPDisentangle simulation study timing",
  paste("Start:    ", timing_report$start_iso),
  paste("End:      ", timing_report$end_iso),
  paste("Elapsed:  ", round(timing_report$elapsed_sec, 1), "sec (", timing_report$elapsed_min, "min)"),
  paste("Job ID:   ", timing_report$job_id),
  paste("Cores:    ", timing_report$n_cores),
  paste("Sim size: ", timing_report$sim_size)
), timing_txt)
log_msg("Timing report written to:", timing_txt)
log_msg("=== TOTAL ELAPSED:", round(elapsed_sec, 1), "s (", round(elapsed_sec / 60, 1), "min) ===")
close(log_con)
