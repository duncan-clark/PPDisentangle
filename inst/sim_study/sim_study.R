#!/usr/bin/env Rscript
# Unified simulation study for the PPDisentangle package.
# Detects SLURM environment automatically; falls back to local defaults.
# Usage:
#   Local:   Rscript sim_study.R [--small]
#   Cluster: sbatch sim_study.slurm

library(PPDisentangle)
library(spatstat)
library(ggplot2)
library(dplyr)
library(data.table)
library(parallel)
library(doParallel)
library(R.utils)
library(caret)

args <- commandArgs(trailingOnly = TRUE)
SMALL <- "--small" %in% args

ON_CLUSTER <- nzchar(Sys.getenv("SLURM_JOB_ID"))

if (ON_CLUSTER) {
  N_CORES      <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
  SIM_SIZE     <- 100
  N_SIMS       <- 10
  N_TAU_SIMS   <- 10
  N_TAU_I      <- 10
  N_TAU_I_TRUE <- 100
  N_PROPOSALS  <- 1000
  EM_ITER      <- 1000
  SEM_EM_ADAPTIVE_ITER <- 1000
  OMEGA        <- c(0, 100, 0, 100)
  END_TIME     <- 110
  TREATMENT_TIME <- 10
  NX <- 10; NY <- 10
  SAVE_DIR <- Sys.getenv("SLURM_SUBMIT_DIR", getwd())
} else if (SMALL) {
  N_CORES      <- 2
  SIM_SIZE     <- 2
  N_SIMS       <- 3
  N_TAU_SIMS   <- 3
  N_TAU_I      <- 3
  N_TAU_I_TRUE <- 5
  N_PROPOSALS  <- 5
  EM_ITER      <- 10
  SEM_EM_ADAPTIVE_ITER <- 10
  OMEGA        <- c(0, 50, 0, 50)
  END_TIME     <- 30
  TREATMENT_TIME <- 5
  NX <- 4; NY <- 4
  SAVE_DIR <- getwd()
} else {
  N_CORES      <- max(1, parallel::detectCores() - 1)
  SIM_SIZE     <- N_CORES
  N_SIMS       <- 10
  N_TAU_SIMS   <- 10
  N_TAU_I      <- 10
  N_TAU_I_TRUE <- 100
  N_PROPOSALS  <- 1000
  EM_ITER      <- 1000
  SEM_EM_ADAPTIVE_ITER <- 1000
  OMEGA        <- c(0, 100, 0, 100)
  END_TIME     <- 110
  TREATMENT_TIME <- 10
  NX <- 10; NY <- 10
  SAVE_DIR <- getwd()
}

TREAT_PROP <- 0.5
TIME_INT   <- END_TIME - TREATMENT_TIME
MAX_TIME   <- 10000 * (END_TIME * OMEGA[2] * OMEGA[4] / 1e6)

cat("=== Simulation Study Config ===\n")
cat("  Mode:", if (ON_CLUSTER) "CLUSTER" else if (SMALL) "SMALL/LOCAL" else "LOCAL", "\n")
cat("  Cores:", N_CORES, " Sims:", SIM_SIZE, "\n")
cat("  Omega:", OMEGA, " T:", END_TIME, " t*:", TREATMENT_TIME, "\n")
cat("  Grid:", NX, "x", NY, "\n")
cat("  Save:", SAVE_DIR, "\n")
cat("===============================\n")

# ------------------------------------------------------------------
# Helper: create a parallel cluster with PPDisentangle loaded
# ------------------------------------------------------------------
make_cluster <- function(n_cores) {
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, {
    library(PPDisentangle)
    library(spatstat)
    library(data.table)
    library(dplyr)
    library(R.utils)
    library(caret)
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
    "N_PROPOSALS", "EM_ITER", "SEM_EM_ADAPTIVE_ITER"
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

cat("True all/nothing ATE:", round(all_nothing_ATE, 4), "\n")

# Create cluster now that all globals are defined
cl <- make_cluster(N_CORES)
export_globals(cl)

# ------------------------------------------------------------------
# 2. Compute true average one-flip ATE
# ------------------------------------------------------------------
cat("Computing true tau_i ... ")
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
cat(round(proc.time()[3] - t0, 1), "s\n")
cat("True one-flip ATE:", round(true_tau_1, 4), "\n")

# ------------------------------------------------------------------
# 3. Generate observed data
# ------------------------------------------------------------------
cat("Generating", SIM_SIZE, "observed datasets ... ")
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
cat(round(proc.time()[3] - t0, 1), "s\n")
cat("Points per sim:", paste(sapply(obs_data, nrow), collapse = ", "), "\n")

# ------------------------------------------------------------------
# 4. Baseline labellings (oracle + naive)
# ------------------------------------------------------------------
cat("Computing oracle and naive labellings ...\n")
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
cat("Generating", N_PROPOSALS, "labelling proposals per sim ... ")
t0 <- proc.time()[3]
labelling_proposals <- lapply(obs_data, function(s) {
  pre  <- as.data.frame(s) %>% filter(.data$t < TREATMENT_TIME)
  post <- as.data.frame(s) %>% filter(.data$t >= TREATMENT_TIME)
  pre$location_process <- "control"
  pre$inferred_process <- NULL
  Filter(Negate(is.null), lapply(1:N_PROPOSALS, function(i) {
    tryCatch({
      tmp <- simulation_labeling_hawkes_hawkes(
        post, partition = partition,
        partition_process = partition_processes,
        statespace = OMEGA, state_spaces = state_spaces,
        windowT = c(TREATMENT_TIME, END_TIME),
        hawkes_params_control = hawkes_par_1,
        hawkes_params_treated = NULL,
        change_factor = 1, filtration = pre
      )
      pre$inferred_process <- "control"
      rbind(pre, tmp)
    }, error = function(e) NULL)
  }))
})
cat(round(proc.time()[3] - t0, 1), "s\n")

pp_labeled_best_proposal <- lapply(seq_along(labelling_proposals), function(i) {
  props <- labelling_proposals[[i]]
  if (length(props) == 0) return(pp_labeled_naive[[i]])
  accs <- sapply(props, function(y) {
    keep <- which(y$t > TREATMENT_TIME)
    if (length(keep) < 2) return(0)
    tryCatch(
      caret::confusionMatrix(
        as.factor(y$inferred_process[keep]),
        as.factor(y$process[keep])
      )$overall[["Accuracy"]],
      error = function(e) 0
    )
  })
  props[[which.max(accs)]]
})

# ------------------------------------------------------------------
# 6. EM-style labelling
# ------------------------------------------------------------------
cat("Running EM-style labelling ... ")
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
cat(round(proc.time()[3] - t0, 1), "s\n")

pp_labelled_em_post <- lapply(EM_max_style, function(x) x$labelling)

# ------------------------------------------------------------------
# 7. Adaptive SEM
# ------------------------------------------------------------------
cat("Running adaptive SEM ... ")
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
    N_labellings = max(10, N_PROPOSALS / 10),
    N_iter = 10, verbose = FALSE,
    adaptive_control = list(
      param_update_cadence = 10,
      proposal_update_cadence = 1,
      state_spaces = NULL,
      update_control_params = FALSE,
      iter = SEM_EM_ADAPTIVE_ITER,
      n_props = 10, change_factor = 0.02,
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
cat(round(proc.time()[3] - t0, 1), "s\n")

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
cat("\n=== Classification Accuracy ===\n")
class_metrics <- lapply(names(labellings), function(nm) {
  accs <- sapply(labellings[[nm]], function(y) {
    keep <- which(y$t > TREATMENT_TIME)
    if (length(keep) < 2) return(NA_real_)
    caret::confusionMatrix(
      as.factor(y$inferred_process[keep]),
      as.factor(y$process[keep])
    )$overall[["Accuracy"]]
  })
  data.frame(method = nm, mean_accuracy = mean(accs, na.rm = TRUE))
})
class_metrics <- do.call(rbind, class_metrics)
print(class_metrics)

# ------------------------------------------------------------------
# 10. ATE estimation
# ------------------------------------------------------------------
cat("Estimating ATEs ... ")
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

task_function <- function(task) {
  tryCatch(
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
}

results_flat <- parLapply(cl, tasks, fun = task_function)
stopCluster(cl)
cat(round(proc.time()[3] - t0, 1), "s\n")

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

cat("\n=== ATE Results ===\n")
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
  cat("No valid ATE results (all tasks timed out or errored).\n")
}

cat("\nTrue all/nothing ATE:", round(all_nothing_ATE, 4), "\n")
cat("True one-flip ATE:  ", round(true_tau_1, 4), "\n")

outfile <- file.path(SAVE_DIR, paste0("sim_study_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData"))
save.image(outfile)
cat("Results saved to:", outfile, "\n")
