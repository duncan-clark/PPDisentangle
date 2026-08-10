# Background exposure and precomp consistency for the bivariate ETAS loglik.
#
# The treated background is switched on only from
# treated_background_zero_before onward, so the mu_1 compensator must be
# charged over that span only (not the full window), and precomputed
# background weights must give the same likelihood as the direct path.

.biv_test_params <- function(mu_1 = 1) {
  c(
    mu_0 = 2, mu_1 = mu_1,
    A_00 = 0.1, alpha_m_00 = 0.5,
    A_11 = 0.1, alpha_m_11 = 0.5,
    A_01 = 0.02, alpha_m_01 = 0.5,
    A_10 = 0.02, alpha_m_10 = 0.5,
    c = 0.1, p = 2.1, D = 1, gamma = 0.5, q = 1.6
  )
}

test_that("mu_1 compensator is charged only over the treated background exposure", {
  skip_if_not_installed("spatstat.geom")
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  ctrl <- spatstat.geom::owin(c(0, 5), c(0, 10))
  treat <- spatstat.geom::owin(c(5, 10), c(0, 10))
  # All events control-labelled so the event term never involves mu_1:
  # the loglik difference in mu_1 isolates the compensator charge.
  dat <- data.frame(
    x = c(1, 2, 3, 4, 2, 3),
    y = c(2, 4, 6, 8, 5, 7),
    t = c(-40, -20, -5, 10, 40, 90),
    mag = c(2.6, 2.8, 3.0, 2.7, 2.9, 3.1),
    inferred_process = "control"
  )
  ll_args <- list(
    realiz = dat, windowT = c(-50, 100), windowS = win, m0 = 2.5,
    control_state_space = ctrl, treated_state_space = treat,
    beta_gr = 2.3
  )

  # Treated background activates at t = 0: exposure is 100 days, not 150.
  ll_lo <- do.call(loglik_etas_bivariate, c(
    list(params = .biv_test_params(mu_1 = 1)), ll_args,
    list(treated_background_zero_before = 0)
  ))
  ll_hi <- do.call(loglik_etas_bivariate, c(
    list(params = .biv_test_params(mu_1 = 3)), ll_args,
    list(treated_background_zero_before = 0)
  ))
  expect_equal(ll_lo - ll_hi, 2 * 100, tolerance = 1e-8)

  # Without a cutoff the exposure is the full 150-day window.
  ll_lo_full <- do.call(loglik_etas_bivariate, c(
    list(params = .biv_test_params(mu_1 = 1)), ll_args
  ))
  ll_hi_full <- do.call(loglik_etas_bivariate, c(
    list(params = .biv_test_params(mu_1 = 3)), ll_args
  ))
  expect_equal(ll_lo_full - ll_hi_full, 2 * 150, tolerance = 1e-8)

  # Cutoff before the window start is a no-op; after the end zeroes exposure.
  ll_lo_early <- do.call(loglik_etas_bivariate, c(
    list(params = .biv_test_params(mu_1 = 1)), ll_args,
    list(treated_background_zero_before = -60)
  ))
  expect_equal(ll_lo_early, ll_lo_full, tolerance = 1e-10)
  ll_lo_late <- do.call(loglik_etas_bivariate, c(
    list(params = .biv_test_params(mu_1 = 1)), ll_args,
    list(treated_background_zero_before = 200)
  ))
  ll_hi_late <- do.call(loglik_etas_bivariate, c(
    list(params = .biv_test_params(mu_1 = 3)), ll_args,
    list(treated_background_zero_before = 200)
  ))
  expect_equal(ll_lo_late - ll_hi_late, 0, tolerance = 1e-8)
})

test_that("precomp path matches the direct path with covariate and cutoff", {
  skip_if_not_installed("spatstat.geom")
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  ctrl <- spatstat.geom::owin(c(0, 5), c(0, 10))
  treat <- spatstat.geom::owin(c(5, 10), c(0, 10))
  # Includes a treated-labelled pre-cutoff event so the event-side time mask
  # is exercised, plus a KDE-style covariate column W.
  dat <- data.frame(
    x = c(1, 6, 2, 7, 3, 8),
    y = c(2, 3, 5, 6, 8, 9),
    t = c(-30, -10, 5, 20, 60, 80),
    mag = c(2.6, 3.2, 2.8, 3.0, 2.7, 2.9),
    inferred_process = c("control", "treated", "control", "treated",
                         "control", "treated"),
    W = c(0.5, 1.5, 0.8, 1.2, 1.0, 0.9)
  )
  params <- .biv_test_params(mu_1 = 1.5)
  windowT <- c(-50, 100)
  cutoff <- 0

  ll_direct <- loglik_etas_bivariate(
    params = params, realiz = dat, windowT = windowT, windowS = win,
    m0 = 2.5, control_state_space = ctrl, treated_state_space = treat,
    background_rate_var = "W", treated_background_zero_before = cutoff,
    beta_gr = 2.3
  )

  # Precomputed weights built the same way fit_etas_bivariate builds them:
  # region masks x covariate, then the activation-time mask.
  W_0 <- as.numeric(!spatstat.geom::inside.owin(dat$x, dat$y, treat)) * dat$W
  W_1 <- as.numeric(!spatstat.geom::inside.owin(dat$x, dat$y, ctrl)) * dat$W
  W_1[dat$t < cutoff] <- 0
  precomp <- list(
    W_0 = W_0, W_1 = W_1,
    areaS_0 = spatstat.geom::area(ctrl),
    areaS_1 = spatstat.geom::area(treat)
  )
  ll_precomp <- loglik_etas_bivariate(
    params = params, realiz = dat, windowT = windowT, windowS = win,
    m0 = 2.5, control_state_space = ctrl, treated_state_space = treat,
    treated_background_zero_before = cutoff,
    beta_gr = 2.3, precomp = precomp
  )
  expect_equal(ll_precomp, ll_direct, tolerance = 1e-10)

  # fit_etas_bivariate's internal precomp must agree with the direct path
  # too: with maxit = 1 the reported value corresponds to (projected)
  # parameters evaluated through its precomp; just check it is finite and
  # that a one-step fit runs with the mask + covariate enabled.
  fit <- fit_etas_bivariate(
    params, dat, windowT, win, m0 = 2.5,
    control_state_space = ctrl, treated_state_space = treat,
    background_rate_var = "W", treated_background_zero_before = cutoff,
    beta_gr = 2.3, max_branching_radius = 0.98, maxit = 1
  )
  expect_true(is.finite(fit$value))
  expect_gt(fit$value, -1e14)
})

test_that("batched loglik applies the same mu_1 exposure as the single path", {
  skip_if_not_installed("spatstat.geom")
  set.seed(20260810)
  n <- 30
  t <- sort(runif(n, -20, 20))
  x <- runif(n, 0, 10)
  y <- runif(n, 0, 10)
  mag <- 2.5 + rexp(n, 1.5)
  windowT <- c(-20, 20)
  windowS <- spatstat.geom::owin(c(0, 10), c(0, 10))
  ctrl <- spatstat.geom::owin(c(0, 5), c(0, 10))
  treat <- spatstat.geom::owin(c(5, 10), c(0, 10))
  cutoff <- 0
  W0 <- as.numeric(!spatstat.geom::inside.owin(x, y, treat))
  W1 <- as.numeric(!spatstat.geom::inside.owin(x, y, ctrl))
  W1[t < cutoff] <- 0
  area0 <- spatstat.geom::area(ctrl)
  area1 <- spatstat.geom::area(treat)
  params <- .biv_test_params(mu_1 = 0.6)

  K <- 4L
  pid_mat <- matrix(rbinom(n * K, 1, 0.5), n, K)

  batch <- loglik_etas_bivariate_batch(
    params = params, t = t, x = x, y = y, mag = mag,
    process_ids = pid_mat,
    W0s = matrix(W0, n, K), W1s = matrix(W1, n, K),
    areaS_0 = area0, areaS_1 = area1, windowT = windowT,
    treated_background_zero_before = cutoff,
    t_already_shifted = FALSE, m0 = 2.5, beta_gr = 2.3,
    max_branching_radius = 0.98, t_trunc = 5
  )
  single <- vapply(seq_len(K), function(k) {
    realiz <- data.frame(
      x = x, y = y, t = t, mag = mag,
      inferred_process = ifelse(pid_mat[, k] == 1L, "treated", "control")
    )
    loglik_etas_bivariate(
      params = params, realiz = realiz, windowT = windowT, windowS = windowS,
      control_state_space = ctrl, treated_state_space = treat,
      treated_background_zero_before = cutoff,
      precomp = list(W_0 = W0, W_1 = W1, areaS_0 = area0, areaS_1 = area1),
      m0 = 2.5, beta_gr = 2.3, max_branching_radius = 0.98, t_trunc = 5
    )
  }, 0)
  expect_equal(batch, single, tolerance = 1e-9)

  # The batch wrapper requires a window anchor for the cutoff.
  expect_error(
    loglik_etas_bivariate_batch(
      params = params, t = t - windowT[1], x = x, y = y, mag = mag,
      process_ids = pid_mat,
      W0s = matrix(W0, n, K), W1s = matrix(W1, n, K),
      areaS_0 = area0, areaS_1 = area1, t_max = diff(windowT),
      treated_background_zero_before = cutoff,
      t_already_shifted = TRUE, m0 = 2.5, beta_gr = 2.3
    ),
    "requires windowT"
  )
})
