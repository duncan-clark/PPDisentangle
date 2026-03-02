#!/usr/bin/env Rscript
# Profiled simulation study for PPDisentangle

# --- Find package root ---
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  return(getwd())
}
SCRIPT_DIR <- get_script_dir()
PKG_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))

devtools::load_all(PKG_ROOT)
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
FORCE_CLUSTER <- "--cluster" %in% args

sims_arg <- grep("--sims", args)
N_SIMS_OVERRIDE <- if (length(sims_arg) > 0 && length(args) >= sims_arg + 1) as.numeric(args[sims_arg + 1]) else NULL

ON_CLUSTER <- nzchar(Sys.getenv("SLURM_JOB_ID")) || FORCE_CLUSTER

if (ON_CLUSTER) {
  N_CORES      <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "100"))
  N_SIMS       <- if (!is.null(N_SIMS_OVERRIDE)) N_SIMS_OVERRIDE else 100
  if (!is.null(N_SIMS_OVERRIDE)) N_CORES <- min(N_CORES, N_SIMS)
  SEM_EM_ADAPTIVE_ITER <- 1000
  SEM_N_ITER   <- 10
  SAVE_DIR     <- file.path(Sys.getenv("SLURM_SUBMIT_DIR", getwd()), "cluster_output_profiled")
} else if (SMALL) {
  N_CORES <- 2; N_SIMS <- 4; SEM_EM_ADAPTIVE_ITER <- 50; SEM_N_ITER <- 2
  SAVE_DIR <- getwd()
} else {
  N_CORES <- max(1, parallel::detectCores() - 1); N_SIMS <- 10
  SEM_EM_ADAPTIVE_ITER <- 100; SEM_N_ITER <- 5
  SAVE_DIR <- getwd()
}

OMEGA <- c(0, 100, 0, 100); END_TIME <- 110; TREATMENT_TIME <- 10
NX <- 10; NY <- 10; EM_ITER <- 100; SEM_N_LABELLINGS <- 10; TREAT_PROP <- 0.5
TIME_INT <- END_TIME - TREATMENT_TIME

dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Config: %s mode | %d sims | %d cores | %d SEM iter | save: %s\n",
    if (ON_CLUSTER) "CLUSTER" else if (SMALL) "SMALL" else "LOCAL",
    N_SIMS, N_CORES, SEM_EM_ADAPTIVE_ITER, SAVE_DIR))

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

# --- Parallel Loop ---
cl <- makeCluster(N_CORES, outfile = file.path(SAVE_DIR, "worker_log.txt"))
registerDoParallel(cl)
clusterExport(cl, ls())
clusterEvalQ(cl, devtools::load_all(PKG_ROOT))

run_single_sim <- function(i, dat) {
  t0 <- proc.time()[3]
  cat(sprintf("[Sim %d] started %s\n", i, Sys.time()))

  total_points <- sum(dat$location_process == "treated" & dat$t >= TREATMENT_TIME)
  params_init  <- list(mu = total_points / TIME_INT, alpha = 0.1, beta = TIME_INT / 10, K = 0.1)

  em_res <- em_style_labelling(
    pp_data = dat, partition = partition, partition_processes = partition_processes,
    statespace = OMEGA, state_spaces = state_spaces, time_window = c(TREATMENT_TIME, END_TIME),
    treatment_time = TREATMENT_TIME, hawkes_params_control = hawkes_par_1, hawkes_params_treated = params_init,
    iter = EM_ITER, n_props = 10, verbose = FALSE
  )

  cat(sprintf("[Sim %d] EM done, starting SEM\n", i))

  sem_res <- adaptive_SEM(
    pp_data = dat, partition = partition, partition_processes = partition_processes,
    statespace = OMEGA, time_window = c(TREATMENT_TIME, END_TIME), treatment_time = TREATMENT_TIME,
    hawkes_params_control = hawkes_par_1, hawkes_params_treated = params_init,
    N_labellings = SEM_N_LABELLINGS, N_iter = SEM_N_ITER, verbose = FALSE,
    adaptive_control = list(iter = SEM_EM_ADAPTIVE_ITER, n_props = 10, verbose = FALSE)
  )

  cat(sprintf("[Sim %d] SEM done, starting ATE\n", i))

  ate_res <- ATE_estim_hawkes(
    statespace = OMEGA, partition = partition, observed_data = sem_res$adaptive$adaptive_labelling,
    treated_partitions = treated_partitions, hawkes_params = list(control = sem_res$hawkes_params_control, treated = sem_res$hawkes_params_treated),
    n_sims = 10, n_tau_sims = 5, n_tau_i = 5, windowT = c(TREATMENT_TIME, END_TIME), windowS = OMEGA
  )

  elapsed <- proc.time()[3] - t0
  cat(sprintf("[Sim %d] finished in %.1fs\n", i, elapsed))

  res <- list(em = em_res, sem = sem_res, ate = ate_res, time = elapsed)
  saveRDS(res, file.path(SAVE_DIR, paste0("res_sim_", i, ".rds")))
  return(res)
}

cat(sprintf("Running %d sims on %d cores...\n", N_SIMS, N_CORES))
t_all_start <- proc.time()[3]
results <- foreach(i = 1:N_SIMS) %dopar% {
  run_single_sim(i, obs_data[[i]])
}
t_all_end <- proc.time()[3]
stopCluster(cl)

total_elapsed <- t_all_end - t_all_start
avg_time <- mean(sapply(results, function(x) x$time))
actual_iter <- if (SMALL) 50 else 1000
scale_factor <- 5000 / actual_iter
est_sim_time <- avg_time * scale_factor
est_total <- (est_sim_time * ceiling(100 / N_CORES)) / 3600

cat(sprintf("\n=== Timing ===\n"))
cat(sprintf("  Total: %.1f min (%d sims, %d cores)\n", total_elapsed / 60, N_SIMS, N_CORES))
cat(sprintf("  Per sim: %.1fs (at %d iter)\n", avg_time, SEM_EM_ADAPTIVE_ITER))
cat(sprintf("  Est full run (100 sims, 5000 iter, %d cores): %.1f hours (%.1f days)\n", N_CORES, est_total, est_total / 24))
if (est_total > 72) cat("  WARNING: >3 days!\n") else cat("  Looks feasible.\n")
