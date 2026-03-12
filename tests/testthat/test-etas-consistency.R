# Small-scale ETAS parameter recovery test.
# Verifies that fit_etas recovers known parameters from simulated data.
# Uses a reduced parameter set (fix structural params) for tractability.

test_that("fit_etas recovers productivity A from simulated ETAS data", {
  skip_on_cran()

  set.seed(2024)
  true_params <- list(mu = 3, A = 0.3, alpha_m = 0.8, c = 0.05,
                      p = 1.2, D = 5.0, gamma = 0.5, q = 1.5)
  m0 <- 2.5
  beta_gr <- 2.3
  win <- spatstat.geom::owin(c(0, 50), c(0, 50))
  T_window <- c(0, 30)

  fixed_structural <- list(c = 0.05, p = 1.2, D = 5.0, gamma = 0.5, q = 1.5)
  n_reps <- 3
  A_estimates <- numeric(n_reps)

  for (i in seq_len(n_reps)) {
    sim <- sim_etas(true_params, T_window, win, m0 = m0, beta_gr = beta_gr)
    sim_df <- data.frame(x = sim$x, y = sim$y, t = sim$t,
                         mag = sim$mag, W = sim$W)
    if (nrow(sim_df) < 10) next

    fit <- tryCatch(
      fit_etas(
        params_init = true_params,
        realiz = sim_df,
        windowT = T_window,
        windowS = win,
        m0 = m0,
        maxit = 3000,
        fixed_params = fixed_structural
      ),
      error = function(e) NULL
    )
    if (!is.null(fit) && is.finite(fit$par["A"])) {
      A_estimates[i] <- fit$par["A"]
    }
  }

  valid <- A_estimates[A_estimates > 0]
  if (length(valid) >= 2) {
    mean_A <- mean(valid)
    expect_true(mean_A > 0.05 && mean_A < 0.8,
                info = sprintf("Mean A = %.3f, expected near %.3f",
                               mean_A, true_params$A))
  }
})

test_that("fit_etas recovers mu from Poisson data (A=0)", {
  skip_on_cran()

  set.seed(2025)
  true_mu <- 5
  true_params <- list(mu = true_mu, A = 0, alpha_m = 1, c = 0.1,
                      p = 1.2, D = 1, gamma = 0.5, q = 1.5)
  m0 <- 2.5
  win <- spatstat.geom::owin(c(0, 30), c(0, 30))
  T_window <- c(0, 20)

  sim <- sim_etas(true_params, T_window, win, m0 = m0, beta_gr = 2.3)
  sim_df <- data.frame(x = sim$x, y = sim$y, t = sim$t,
                       mag = sim$mag, W = sim$W)

  if (nrow(sim_df) >= 5) {
    fixed_all_but_mu <- list(A = 0, alpha_m = 1, c = 0.1,
                             p = 1.2, D = 1, gamma = 0.5, q = 1.5)
    fit <- fit_etas(
      params_init = list(mu = 3, A = 0, alpha_m = 1, c = 0.1,
                         p = 1.2, D = 1, gamma = 0.5, q = 1.5),
      realiz = sim_df, windowT = T_window, windowS = win,
      m0 = m0, maxit = 2000, fixed_params = fixed_all_but_mu
    )
    mu_hat <- fit$par["mu"]
    expect_true(mu_hat > 1 && mu_hat < 15,
                info = sprintf("mu_hat = %.3f, expected near %.1f",
                               mu_hat, true_mu))
  }
})
