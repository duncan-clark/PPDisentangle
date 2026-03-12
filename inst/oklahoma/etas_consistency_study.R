#!/usr/bin/env Rscript
# ============================================================================
# ETAS Consistency Study
#
# Verifies that fit_etas recovers known ETAS parameters from simulated data.
# Runs N_REPS independent replicates, each simulating from the true model
# and fitting via MLE with structural parameters optionally held fixed.
#
# Supports --cluster mode (saves to cluster_output/) and --test mode
# (reduced reps for CI).
#
# Usage:
#   Rscript etas_consistency_study.R
#   Rscript etas_consistency_study.R --cluster --test
#   Rscript etas_consistency_study.R --reps 50
# ============================================================================

library(PPDisentangle)
library(spatstat)
library(parallel)

# ---- CLI arguments ----
args <- commandArgs(trailingOnly = TRUE)
CLUSTER_MODE <- "--cluster" %in% args
TEST_MODE    <- "--test" %in% args

get_arg <- function(flag, default) {
  idx <- match(flag, args)
  if (!is.na(idx) && length(args) > idx) return(args[idx + 1])
  return(default)
}
N_REPS <- as.integer(get_arg("--reps", if (TEST_MODE) "3" else "100"))

# ---- Cluster output paths ----
JOB_ID <- Sys.getenv("SLURM_JOB_ID", format(Sys.time(), "local_%Y%m%d_%H%M%S"))
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
N_CORES <- max(1L, min(N_CORES, N_REPS, detectCores()))

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
PKG_ROOT <- tryCatch(
  normalizePath(file.path(SCRIPT_DIR, "../..")),
  error = function(e) getwd()
)

if (CLUSTER_MODE) {
  OUT_DIR <- file.path(PKG_ROOT, "cluster_output", "results")
  LOG_DIR <- file.path(PKG_ROOT, "cluster_output", "logs", JOB_ID)
  if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
  if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
  OUT_FILE <- file.path(OUT_DIR, paste0("etas_consistency_", JOB_ID, ".rds"))
} else {
  OUT_DIR <- file.path(SCRIPT_DIR, "output")
  if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
  OUT_FILE <- file.path(OUT_DIR, "etas_consistency_results.rds")
}

# ---- True parameters ----
# Moderate ETAS with subcritical branching ratio
true_params <- list(mu = 3.0, A = 0.25, alpha_m = 0.8,
                    c = 0.05, p = 1.2, D = 5.0, gamma = 0.5, q = 1.5)
m0 <- 2.5
beta_gr <- 2.3  # b * ln(10) with b ~ 1.0
win <- owin(c(0, 50), c(0, 50))
T_window <- c(0, 100)

# Fix structural parameters (c, p, D, gamma, q) and estimate (mu, A, alpha_m)
fixed_structural <- list(c = 0.05, p = 1.2, D = 5.0, gamma = 0.5, q = 1.5)

MAXIT <- if (TEST_MODE) 1000 else 5000

cat("=== ETAS Consistency Study ===\n")
cat(sprintf("Reps: %d | Cores: %d | Mode: %s\n", N_REPS, N_CORES,
            if (CLUSTER_MODE) "cluster" else "local"))
cat(sprintf("True: mu=%.2f A=%.2f alpha_m=%.2f c=%.3f p=%.1f D=%.1f gamma=%.1f q=%.1f\n",
            true_params$mu, true_params$A, true_params$alpha_m,
            true_params$c, true_params$p, true_params$D,
            true_params$gamma, true_params$q))
cat(sprintf("Fixed: %s\n", paste(names(fixed_structural), fixed_structural, sep="=", collapse=", ")))
cat(sprintf("Output: %s\n\n", OUT_FILE))

# ---- Single replicate function ----
run_one_rep <- function(seed) {
  set.seed(seed)
  sim <- sim_etas(true_params, T_window, win, m0 = m0, beta_gr = beta_gr)
  sim_df <- data.frame(x = sim$x, y = sim$y, t = sim$t,
                       mag = sim$mag, W = sim$W)
  n_pts <- nrow(sim_df)

  if (n_pts < 10) {
    cat(sprintf("  Rep %d: only %d events, skipping\n", seed, n_pts))
    return(list(seed = seed, n_pts = n_pts, par = rep(NA, 8), convergence = -1))
  }

  fit <- tryCatch(
    fit_etas(
      params_init = true_params,
      realiz = sim_df,
      windowT = T_window,
      windowS = win,
      m0 = m0,
      maxit = MAXIT,
      fixed_params = fixed_structural
    ),
    error = function(e) {
      cat(sprintf("  Rep %d: error - %s\n", seed, e$message))
      NULL
    }
  )

  if (!is.null(fit) && all(is.finite(fit$par))) {
    cat(sprintf("  Rep %d: n=%d, mu=%.3f A=%.3f alpha_m=%.3f (conv=%d)\n",
                seed, n_pts, fit$par["mu"], fit$par["A"], fit$par["alpha_m"],
                fit$convergence))
    return(list(seed = seed, n_pts = n_pts, par = fit$par,
                convergence = fit$convergence, value = fit$value))
  } else {
    return(list(seed = seed, n_pts = n_pts, par = rep(NA, 8), convergence = -1))
  }
}

# ---- Run replicates ----
seeds <- seq_len(N_REPS) + 2024L
t_start <- proc.time()[3]

if (N_CORES > 1) {
  results <- mclapply(seeds, run_one_rep, mc.cores = N_CORES)
} else {
  results <- lapply(seeds, run_one_rep)
}

elapsed <- proc.time()[3] - t_start
cat(sprintf("\nCompleted %d reps in %.1f seconds\n", N_REPS, elapsed))

# ---- Summarize ----
valid <- Filter(function(r) all(is.finite(r$par)), results)
cat(sprintf("Valid fits: %d / %d\n\n", length(valid), N_REPS))

if (length(valid) >= 2) {
  res_mat <- do.call(rbind, lapply(valid, function(r) r$par))
  colnames(res_mat) <- names(true_params)

  cat("--- Results (Mean +/- SD) ---\n")
  means <- colMeans(res_mat)
  sds <- apply(res_mat, 2, sd)
  true_vec <- unlist(true_params)

  for (nm in names(means)) {
    cat(sprintf("  %8s: True=%.3f  Est=%.3f (+/- %.3f)  Bias=%.3f\n",
                nm, true_vec[nm], means[nm], sds[nm], means[nm] - true_vec[nm]))
  }

  free_names <- setdiff(names(true_params), names(fixed_structural))
  free_means <- means[free_names]
  free_true  <- true_vec[free_names]
  all_close <- all(abs(free_means - free_true) / pmax(abs(free_true), 0.01) < 0.5)
  cat(sprintf("\nVERDICT: %s\n", if (all_close) "SUCCESS - free parameters recovered" else "WARNING - some parameters may be biased"))
}

# ---- Save ----
output <- list(
  results = results,
  true_params = true_params,
  fixed_structural = fixed_structural,
  m0 = m0, beta_gr = beta_gr,
  config = list(N_REPS = N_REPS, MAXIT = MAXIT, TEST_MODE = TEST_MODE,
                win = win, T_window = T_window),
  elapsed_s = elapsed
)
saveRDS(output, OUT_FILE)
cat(sprintf("Saved to: %s\n", OUT_FILE))
