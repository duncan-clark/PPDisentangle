
# Consistency Study for Hawkes Fit
library(PPDisentangle)
library(spatstat)
library(Rcpp)
library(data.table)
library(parallel)

# 1. Configuration
set.seed(42)
true_params <- list(mu = 0.5, alpha = 0.1, beta = 0.2, K = 0.4)
win <- owin(c(0, 10), c(0, 10))
T_window <- c(0, 200) 
N_reps <- 2

cat("Starting Consistency Study with", N_reps, "replications...\n")
cat(sprintf("True Parameters: mu=%.2f, alpha=%.2f, beta=%.2f, K=%.2f\n\n", 
            true_params$mu, true_params$alpha, true_params$beta, true_params$K))

results <- list()

for (i in 1:N_reps) {
  cat("Repetition", i, "...")
  
  # Simulate
  sim_data <- sim_hawkes_fast(params = true_params, windowT = T_window, windowS = win)
  sim_df <- as.data.frame(sim_data)
  
  # Fit (using true param as start to test local convergence)
  start_params <- unlist(true_params)
  
  fit <- tryCatch({
    fit_hawkes(
      params_init = start_params,
      realiz = sim_df,
      windowT = T_window,
      windowS = win,
      use_fast = TRUE,
      method = "Nelder-Mead",
      maxit = 1000
    )
  }, error = function(e) NULL)
  
  if (!is.null(fit)) {
    results[[i]] <- fit$par
    cat(" Done. Estimated K:", round(fit$par[4], 3), " mu:", round(fit$par[1], 3), " N_points:", nrow(sim_df), "\n")
  } else {
    cat(" Failed.\n")
  }
}

# 2. Summary
res_mat <- do.call(rbind, results)
colnames(res_mat) <- c("mu", "alpha", "beta", "K")

cat("\n--- Consistency Study Results (Mean +/- SD) ---\n")
means <- colMeans(res_mat)
sds <- apply(res_mat, 2, sd)

for (p in names(means)) {
  cat(sprintf("%s: True=%.2f, Est=%.3f (+/- %.3f)\n", 
              p, true_params[[p]], means[p], sds[p]))
}

# 3. Final Check
if (abs(means["K"] - true_params$K) < 0.15) {
  cat("\nVERDICT: SUCCESS - K is reasonably recovered.\n")
} else {
  cat("\nVERDICT: FAILURE - K recovery is biased.\n")
}
