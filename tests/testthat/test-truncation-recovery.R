test_that("parameter recovery is equivalent with high truncation vs no truncation", {
  # Set up parameters - use higher K and beta for better identifiability
  true_params <- list(mu = 1.0, alpha = 0.5, beta = 2.0, K = 0.6)
  windowT <- c(0, 50)
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  
  # 1. Simulate data without truncation
  set.seed(123)
  sim_data <- sim_hawkes_fast(true_params, windowT, windowS)
  realiz <- data.frame(x = sim_data$x, y = sim_data$y, t = sim_data$t, W = 1.0)
  
  # 1b. Direct likelihood equivalence test
  ll_no_trunc   <- loglik_hawk_fast(true_params, realiz, windowT, windowS, t_trunc = NULL)
  ll_high_trunc <- loglik_hawk_fast(true_params, realiz, windowT, windowS, t_trunc = 1e9)
  expect_equal(ll_no_trunc, ll_high_trunc, tolerance = 1e-8)

  # 2. Fit without truncation (t_trunc = NULL)
  set.seed(456)
  fit_no_trunc <- fit_hawkes(
    params_init = true_params,
    realiz = realiz,
    windowT = windowT,
    windowS = windowS,
    background_rate_var = "W",
    use_fast = TRUE,
    maxit = 500
  )
  
  # 3. Fit with very high truncation (t_trunc = 1e9)
  set.seed(456)
  fit_high_trunc <- fit_hawkes(
    params_init = true_params,
    realiz = realiz,
    windowT = windowT,
    windowS = windowS,
    background_rate_var = "W",
    use_fast = TRUE,
    t_trunc = 1e9,
    maxit = 500
  )
  
  # Compare parameter estimates
  # They should be very close to each other
  expect_equal(fit_no_trunc$par, fit_high_trunc$par, tolerance = 1e-5)
  
  # And they should be reasonably close to true parameters
  expect_equal(unname(fit_no_trunc$par["K"]), true_params$K, tolerance = 0.3)
})
