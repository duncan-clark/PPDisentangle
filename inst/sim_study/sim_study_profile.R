#!/usr/bin/env Rscript
# Profiled simulation study for PPDisentangle
# Adds better logging, intermediate saving, and runtime estimation.

devtools::load_all("../..")
library(spatstat)
library(ggplot2)
library(dplyr)
library(data.table)
library(parallel)
library(doParallel)
library(R.utils)
library(caret)
library(reshape2)

# --- Configuration ---
args <- commandArgs(trailingOnly = TRUE)
SMALL <- "--small" %in% args
PILOT <- "--pilot" %in% args
FORCE_CLUSTER <- "--cluster" %in% args

# Parse --sims argument if present
sims_arg <- grep("--sims", args)
N_SIMS_OVERRIDE <- if (length(sims_arg) > 0 && length(args) >= sims_arg + 1) as.numeric(args[sims_arg + 1]) else NULL

ON_CLUSTER <- nzchar(Sys.getenv("SLURM_JOB_ID")) || FORCE_CLUSTER

if (ON_CLUSTER) {
  N_CORES      <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "100"))
  N_SIMS       <- if (!is.null(N_SIMS_OVERRIDE)) N_SIMS_OVERRIDE else 100
  # If we are overriding sims, also cap cores to N_SIMS to avoid idle workers
  if (!is.null(N_SIMS_OVERRIDE)) N_CORES <- min(N_CORES, N_SIMS)
  
  SEM_EM_ADAPTIVE_ITER <- 1000
  SEM_N_ITER   <- 10
  SAVE_DIR     <- file.path(Sys.getenv("SLURM_SUBMIT_DIR", getwd()), "cluster_output_profiled")
  LOG_DIR      <- file.path(SAVE_DIR, "logs")
} else if (SMALL) {
  N_CORES      <- 2
  N_SIMS       <- 4
  SEM_EM_ADAPTIVE_ITER <- 50
  SEM_N_ITER   <- 2
  SAVE_DIR     <- getwd()
  LOG_DIR      <- file.path(SAVE_DIR, "logs")
} else {
  N_CORES      <- max(1, parallel::detectCores() - 1)
  N_SIMS       <- 10
  SEM_EM_ADAPTIVE_ITER <- 100
  SEM_N_ITER   <- 5
  SAVE_DIR     <- getwd()
  LOG_DIR      <- file.path(SAVE_DIR, "logs")
}

if (PILOT) {
  cat("!!! PILOT MODE ENABLED !!!\n")
  N_SIMS <- N_CORES # Just one batch
}

# Fixed parameters for the study
OMEGA        <- c(0, 100, 0, 100)
END_TIME     <- 110
TREATMENT_TIME <- 10
NX <- 10; NY <- 10
N_TAU_I      <- 10
N_TAU_SIMS   <- 10
N_TAU_I_TRUE <- 100
N_PROPOSALS  <- 100
EM_ITER      <- 100
SEM_N_LABELLINGS <- 10
TREAT_PROP <- 0.5
TIME_INT   <- END_TIME - TREATMENT_TIME

dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Profiled Simulation Study Config ===\n")
cat("  Mode:", if (ON_CLUSTER) "CLUSTER" else if (SMALL) "SMALL" else "LOCAL", "\n")
cat("  Cores:", N_CORES, " Sims:", N_SIMS, "\n")
cat("  SEM Adaptive Iter:", SEM_EM_ADAPTIVE_ITER, " SEM Outer Iter:", SEM_N_ITER, "\n")
cat("  Save Dir:", SAVE_DIR, "\n")
cat("  Log Dir:", LOG_DIR, "\n")
cat("========================================\n")

# --- Setup ---
hawkes_par_1 <- list(mu = 8, alpha = 0.01, beta = 10, K = 0.8)
hawkes_par_2 <- list(mu = 8, alpha = 0.01, beta = 10, K = 0.2)
partition <- quadrats(X = OMEGA, nx = NX, ny = NY)
partition_processes <- rep("control", partition$n)
set.seed(42)
partition_processes[sample(1:(NX * NY), NX * NY * TREAT_PROP)] <- "treated"
treated_idx <- partition_processes == "treated"
control_state_space <- as.owin(partition[!treated_idx])
treated_state_space <- as.owin(partition[treated_idx])
treated_partitions  <- tilenames(partition[which(partition_processes == "treated")])
state_spaces <- if (partition_processes[1] == "treated") list(treated_state_space, control_state_space) else list(control_state_space, treated_state_space)

# --- Data Generation ---
cat("Generating datasets...\n")
set.seed(123)
obs_data <- lapply(1:N_SIMS, function(i) {
  pre_treat <- sim_hawkes(params = hawkes_par_1, windowT = c(0, TREATMENT_TIME), windowS = OMEGA, optimized = TRUE)
  pre_df <- data.frame(x = pre_treat$x, y = pre_treat$y, t = pre_treat$t, background = pre_treat$background, process = "control", location_process = "control")
  combined <- generate_inhomogeneous_hawkes(Omega = OMEGA, partition = partition, time_window = c(TREATMENT_TIME, END_TIME), partition_processes = partition_processes, state_spaces = state_spaces, filtration = pre_df, hawkes_params = list(control = hawkes_par_1, treated = hawkes_par_2))
  post_df <- as.data.frame(combined)
  rbind(pre_df, post_df[, names(pre_df)])
})

# --- Parallel Loop with Logging ---
cl <- makeCluster(N_CORES, outfile = file.path(LOG_DIR, "worker_log.txt"))
registerDoParallel(cl)
clusterExport(cl, ls())
clusterEvalQ(cl, devtools::load_all("../.."))

cat("Starting parallel simulations. Check logs/worker_log.txt for detailed progress.\n")

# We'll use a wrapper that logs start/end and saves intermediate results
run_single_sim <- function(i, dat) {
  log_file <- file.path(LOG_DIR, paste0("sim_", i, "_progress.txt"))
  cat(paste0("Sim ", i, " started at ", Sys.time(), "\n"), file = log_file)
  
  t0 <- proc.time()[3]
  
  # 1. EM-style
  cat(paste0("Sim ", i, ": Running EM...\n"), file = log_file, append = TRUE)
  total_points <- sum(dat$location_process == "treated" & dat$t >= TREATMENT_TIME)
  params_init  <- list(mu = total_points / TIME_INT, alpha = 0.1, beta = TIME_INT / 10, K = 0.1)
  
  em_res <- em_style_labelling(
    pp_data = dat, partition = partition, partition_processes = partition_processes,
    statespace = OMEGA, state_spaces = state_spaces, time_window = c(TREATMENT_TIME, END_TIME),
    treatment_time = TREATMENT_TIME, hawkes_params_control = hawkes_par_1, hawkes_params_treated = params_init,
    iter = EM_ITER, n_props = 10, verbose = FALSE
  )
  
  # 2. Adaptive SEM
  cat(paste0("Sim ", i, ": Running SEM...\n"), file = log_file, append = TRUE)
  sem_res <- adaptive_SEM(
    pp_data = dat, partition = partition, partition_processes = partition_processes,
    statespace = OMEGA, time_window = c(TREATMENT_TIME, END_TIME), treatment_time = TREATMENT_TIME,
    hawkes_params_control = hawkes_par_1, hawkes_params_treated = params_init,
    N_labellings = SEM_N_LABELLINGS, N_iter = SEM_N_ITER, verbose = TRUE,
    adaptive_control = list(iter = SEM_EM_ADAPTIVE_ITER, n_props = 10, verbose = TRUE)
  )
  
  # 3. ATE Estimation
  cat(paste0("Sim ", i, ": Estimating ATE...\n"), file = log_file, append = TRUE)
  # (Simplified for profiling)
  ate_res <- ATE_estim_hawkes(
    statespace = OMEGA, partition = partition, observed_data = sem_res$adaptive$adaptive_labelling,
    treated_partitions = treated_partitions, hawkes_params = list(control = sem_res$hawkes_params_control, treated = sem_res$hawkes_params_treated),
    n_sims = 10, n_tau_sims = 5, n_tau_i = 5, windowT = c(TREATMENT_TIME, END_TIME), windowS = OMEGA
  )
  
  t1 <- proc.time()[3]
  elapsed <- t1 - t0
  
  cat(paste0("Sim ", i, " finished in ", round(elapsed, 1), "s at ", Sys.time(), "\n"), file = log_file, append = TRUE)
  
  res <- list(em = em_res, sem = sem_res, ate = ate_res, time = elapsed)
  saveRDS(res, file.path(SAVE_DIR, paste0("res_sim_", i, ".rds")))
  return(res)
}

t_all_start <- proc.time()[3]
results <- foreach(i = 1:N_SIMS) %dopar% {
  run_single_sim(i, obs_data[[i]])
}
t_all_end <- proc.time()[3]
stopCluster(cl)

total_elapsed <- t_all_end - t_all_start
avg_time <- mean(sapply(results, function(x) x$time))

cat("\n=== Timing Summary ===\n")
cat("  Total time for", N_SIMS, "sims on", N_CORES, "cores:", round(total_elapsed / 60, 2), "min\n")
cat("  Average time per sim:", round(avg_time, 2), "s\n")
cat("  Estimated time for 100 sims on", N_CORES, "cores with full iterations:\n")

# Extrapolation logic
# Current avg_time is for SEM_EM_ADAPTIVE_ITER iterations.
# We want to estimate for 5000 iterations.
# ATE and EM also take time, but SEM is the main part.
# Let's assume linear scaling for SEM iterations.

if (SMALL) {
  actual_iter <- 50
} else {
  actual_iter <- 1000
}
target_iter <- 5000
scale_factor <- target_iter / actual_iter

est_sim_time <- avg_time * scale_factor
est_total_time_100_sims_cores <- (est_sim_time * ceiling(100/N_CORES)) / 3600 # hours

cat("  Scale factor:", scale_factor, "\n")
cat("  Estimated time per full sim:", round(est_sim_time / 3600, 2), "hours\n")
cat("  Estimated total time (100 sims,", N_CORES, "cores):", round(est_total_time_100_sims_cores, 2), "hours (", round(est_total_time_100_sims_cores/24, 2), "days )\n")
cat("=======================\n")

if (est_total_time_100_sims_cores > 72) {
  cat("WARNING: Estimated time exceeds 3 days! Consider reducing iterations or increasing cores.\n")
} else {
  cat("SUCCESS: Run seems feasible within 3 days.\n")
}
