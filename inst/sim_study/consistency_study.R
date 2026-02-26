
# Consistency Study for Hawkes Fit
library(PPDisentangle)
library(spatstat)

set.seed(42)

# Window large enough relative to spatial kernel SD (~0.7) to avoid edge effects.
true_params <- list(mu = 5, alpha = 1.0, beta = 1.0, K = 0.4)
win <- owin(c(0, 50), c(0, 50))
T_window <- c(0, 20)
N_reps <- 20

cat("Starting Consistency Study with", N_reps, "replications...\n")
cat(sprintf("True Parameters: mu=%.2f, alpha=%.2f, beta=%.2f, K=%.2f\n\n",
            true_params$mu, true_params$alpha, true_params$beta, true_params$K))

results <- list()

for (i in 1:N_reps) {
  cat("Repetition", i, "...")

  sim_data <- sim_hawkes_fast(params = true_params, windowT = T_window, windowS = win)
  sim_df <- as.data.frame(sim_data)

  start_params <- unlist(true_params)

  fit <- tryCatch({
    fit_hawkes(
      params_init = start_params,
      realiz = sim_df,
      windowT = T_window,
      windowS = win,
      use_fast = TRUE,
      method = "Nelder-Mead",
      maxit = 5000
    )
  }, error = function(e) NULL)

  if (!is.null(fit)) {
    results[[i]] <- fit$par
    cat(" Done. Estimated K:", round(fit$par[4], 3), " mu:", round(fit$par[1], 3), " N_points:", nrow(sim_df), "\n")
  } else {
    cat(" Failed.\n")
  }
}

# Summary
res_mat <- do.call(rbind, results)
colnames(res_mat) <- c("mu", "alpha", "beta", "K")

cat("\n--- Consistency Study Results (Mean +/- SD) ---\n")
means <- colMeans(res_mat)
sds <- apply(res_mat, 2, sd)

for (p in names(means)) {
  cat(sprintf("%s: True=%.2f, Est=%.3f (+/- %.3f)\n",
              p, true_params[[p]], means[p], sds[p]))
}

if (abs(means["K"] - true_params$K) < 0.15) {
  cat("\nVERDICT: SUCCESS - K is reasonably recovered.\n")
} else {
  cat("\nVERDICT: FAILURE - K recovery is biased.\n")
}
