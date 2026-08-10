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

# ---------------------------------------------------------------------------
# Control background everywhere before the policy cutoff (spec: pre-treatment
# everywhere is control, so control events can be background everywhere).
# ---------------------------------------------------------------------------

# Pure-background manual reference (all A = 0): event terms are just the
# background intensities (floored at 1e-15 like the C++ kernel) and the
# compensator is mu_0*expo_0 + mu_1*expo_1.
.biv_manual_bg_ll <- function(params, dat, windowT, ctrl, treat,
                              cutoff0 = NULL, mass = NULL, cutoff1 = NULL) {
  a0 <- spatstat.geom::area(ctrl)
  a1 <- spatstat.geom::area(treat)
  in_treat <- spatstat.geom::inside.owin(dat$x, dat$y, treat)
  W0 <- rep(1, nrow(dat))
  W1 <- rep(1, nrow(dat))
  if (is.null(cutoff0)) W0[in_treat] <- 0 else W0[in_treat & dat$t >= cutoff0] <- 0
  W1[!in_treat] <- 0
  if (!is.null(cutoff1)) W1[dat$t < cutoff1] <- 0
  lam <- ifelse(dat$inferred_process == "treated",
                params[["mu_1"]] * W1 / a1,
                params[["mu_0"]] * W0 / a0)
  lam <- pmax(lam, 1e-15)
  expo1 <- if (is.null(cutoff1)) {
    diff(windowT)
  } else {
    max(0, windowT[2] - max(windowT[1], cutoff1))
  }
  M <- if (is.null(mass)) (a0 + a1) / a0 else mass
  expo0 <- if (is.null(cutoff0)) {
    diff(windowT)
  } else {
    ct <- min(max(cutoff0, windowT[1]), windowT[2])
    M * (ct - windowT[1]) + (windowT[2] - ct)
  }
  sum(log(lam)) - (params[["mu_0"]] * expo0 + params[["mu_1"]] * expo1)
}

.biv_test_params_bg_only <- function(mu_0 = 2, mu_1 = 1) {
  p <- .biv_test_params(mu_1 = mu_1)
  p[["mu_0"]] <- mu_0
  p[c("A_00", "A_11", "A_01", "A_10")] <- 0
  p
}

test_that("pre-cutoff control events in the treated region get background", {
  skip_if_not_installed("spatstat.geom")
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  ctrl <- spatstat.geom::owin(c(0, 5), c(0, 10))
  treat <- spatstat.geom::owin(c(5, 10), c(0, 10))
  # Two pre-cutoff control events INSIDE the treated region, plus regular
  # control-region and post-cutoff treated events.
  dat <- data.frame(
    x = c(7, 8, 2, 3, 7, 2),
    y = c(2, 6, 4, 8, 5, 6),
    t = c(-40, -10, -20, 30, 50, 80),
    mag = c(2.7, 2.9, 2.6, 3.0, 2.8, 3.1),
    inferred_process = c("control", "control", "control", "control",
                         "treated", "control")
  )
  windowT <- c(-50, 100)
  params <- .biv_test_params_bg_only(mu_0 = 2, mu_1 = 1)
  ll_args <- list(
    params = params, realiz = dat, windowT = windowT, windowS = win,
    m0 = 2.5, control_state_space = ctrl, treated_state_space = treat,
    beta_gr = 2.3, treated_background_zero_before = 0
  )

  # Old spec (no control cutoff): those events hit the 1e-15 intensity floor.
  ll_old <- do.call(loglik_etas_bivariate, ll_args)
  ll_old_ref <- .biv_manual_bg_ll(params, dat, windowT, ctrl, treat,
                                  cutoff0 = NULL, cutoff1 = 0)
  expect_equal(ll_old, ll_old_ref, tolerance = 1e-10)

  # New spec: control background is live everywhere before t = 0.
  ll_new <- do.call(loglik_etas_bivariate, c(
    ll_args, list(control_background_everywhere_before = 0)
  ))
  ll_ref <- .biv_manual_bg_ll(params, dat, windowT, ctrl, treat,
                              cutoff0 = 0, cutoff1 = 0)
  expect_equal(ll_new, ll_ref, tolerance = 1e-10)

  # With flat weights and the cutoff on, moving a pre-cutoff control event
  # across the region boundary must not change the likelihood at all.
  dat_moved <- dat
  dat_moved$x[1] <- 2
  ll_moved <- do.call(loglik_etas_bivariate, c(
    list(params = params,
         realiz = dat_moved, windowT = windowT, windowS = win, m0 = 2.5,
         control_state_space = ctrl, treated_state_space = treat,
         beta_gr = 2.3, treated_background_zero_before = 0,
         control_background_everywhere_before = 0)
  ))
  expect_equal(ll_moved, ll_new, tolerance = 1e-12)

  # Explicit (KDE-style) mass ratio changes only the mu_0 compensator charge.
  ll_mass <- do.call(loglik_etas_bivariate, c(
    ll_args, list(control_background_everywhere_before = 0,
                  control_background_pre_mass_ratio = 3.5)
  ))
  ll_mass_ref <- .biv_manual_bg_ll(params, dat, windowT, ctrl, treat,
                                   cutoff0 = 0, mass = 3.5, cutoff1 = 0)
  expect_equal(ll_mass, ll_mass_ref, tolerance = 1e-10)

  # A post-cutoff control event in the treated region must still hit the
  # intensity floor (control background stays off there after the cutoff).
  dat_bad <- dat
  dat_bad$t[1] <- 10
  ll_bad <- do.call(loglik_etas_bivariate, c(
    list(params = params,
         realiz = dat_bad, windowT = windowT, windowS = win, m0 = 2.5,
         control_state_space = ctrl, treated_state_space = treat,
         beta_gr = 2.3, treated_background_zero_before = 0,
         control_background_everywhere_before = 0)
  ))
  ll_bad_ref <- .biv_manual_bg_ll(params, dat_bad, windowT, ctrl, treat,
                                  cutoff0 = 0, cutoff1 = 0)
  expect_equal(ll_bad, ll_bad_ref, tolerance = 1e-10)
  expect_lt(ll_bad, ll_new - 20)
})

test_that("precomp/fit path matches the direct path with the control cutoff", {
  skip_if_not_installed("spatstat.geom")
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  ctrl <- spatstat.geom::owin(c(0, 5), c(0, 10))
  treat <- spatstat.geom::owin(c(5, 10), c(0, 10))
  # Pre-cutoff control event in the treated region carries its KDE weight.
  dat <- data.frame(
    x = c(7, 1, 2, 7, 3, 8),
    y = c(2, 3, 5, 6, 8, 9),
    t = c(-30, -10, 5, 20, 60, 80),
    mag = c(2.6, 3.2, 2.8, 3.0, 2.7, 2.9),
    inferred_process = c("control", "control", "control", "treated",
                         "control", "treated"),
    W = c(0.7, 0.5, 0.8, 1.2, 1.0, 0.9)
  )
  params <- .biv_test_params(mu_1 = 1.5)
  windowT <- c(-50, 100)
  cutoff <- 0
  mass <- 2.8

  ll_direct <- loglik_etas_bivariate(
    params = params, realiz = dat, windowT = windowT, windowS = win,
    m0 = 2.5, control_state_space = ctrl, treated_state_space = treat,
    background_rate_var = "W",
    treated_background_zero_before = cutoff,
    control_background_everywhere_before = cutoff,
    control_background_pre_mass_ratio = mass,
    beta_gr = 2.3
  )
  expect_gt(ll_direct, -1e14)

  # Precomputed weights with the time-aware control mask baked in.
  in_treat <- spatstat.geom::inside.owin(dat$x, dat$y, treat)
  W_0 <- dat$W
  W_0[in_treat & dat$t >= cutoff] <- 0
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
    control_background_everywhere_before = cutoff,
    control_background_pre_mass_ratio = mass,
    beta_gr = 2.3, precomp = precomp
  )
  expect_equal(ll_precomp, ll_direct, tolerance = 1e-10)

  # fit_etas_bivariate bakes the same mask into its internal precomp.
  fit <- fit_etas_bivariate(
    params, dat, windowT, win, m0 = 2.5,
    control_state_space = ctrl, treated_state_space = treat,
    background_rate_var = "W",
    treated_background_zero_before = cutoff,
    control_background_everywhere_before = cutoff,
    control_background_pre_mass_ratio = mass,
    beta_gr = 2.3, max_branching_radius = 0.98, maxit = 1
  )
  expect_true(is.finite(fit$value))
  expect_gt(fit$value, -1e14)
})

test_that("batched loglik applies the same mu_0 exposure as the single path", {
  skip_if_not_installed("spatstat.geom")
  set.seed(20260811)
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
  mass <- 2.2
  in_treat <- spatstat.geom::inside.owin(x, y, treat)
  W0 <- rep(1, n)
  W0[in_treat & t >= cutoff] <- 0
  W1 <- as.numeric(in_treat)
  W1[t < cutoff] <- 0
  area0 <- spatstat.geom::area(ctrl)
  area1 <- spatstat.geom::area(treat)
  params <- .biv_test_params(mu_1 = 0.6)

  K <- 4L
  # Pre-cutoff events must stay control-labelled (carry-over convention).
  pid_mat <- matrix(rbinom(n * K, 1, 0.5), n, K)
  pid_mat[t < cutoff, ] <- 0L

  batch <- loglik_etas_bivariate_batch(
    params = params, t = t, x = x, y = y, mag = mag,
    process_ids = pid_mat,
    W0s = matrix(W0, n, K), W1s = matrix(W1, n, K),
    areaS_0 = area0, areaS_1 = area1, windowT = windowT,
    treated_background_zero_before = cutoff,
    control_background_everywhere_before = cutoff,
    control_background_pre_mass_ratio = mass,
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
      control_background_everywhere_before = cutoff,
      control_background_pre_mass_ratio = mass,
      precomp = list(W_0 = W0, W_1 = W1, areaS_0 = area0, areaS_1 = area1),
      m0 = 2.5, beta_gr = 2.3, max_branching_radius = 0.98, t_trunc = 5
    )
  }, 0)
  expect_equal(batch, single, tolerance = 1e-9)

  expect_error(
    loglik_etas_bivariate_batch(
      params = params, t = t - windowT[1], x = x, y = y, mag = mag,
      process_ids = pid_mat,
      W0s = matrix(W0, n, K), W1s = matrix(W1, n, K),
      areaS_0 = area0, areaS_1 = area1, t_max = diff(windowT),
      control_background_everywhere_before = cutoff,
      t_already_shifted = TRUE, m0 = 2.5, beta_gr = 2.3
    ),
    "requires windowT"
  )
})

test_that("sim_etas_bivariate bg_ref_areas extends the background density", {
  skip_if_not_installed("spatstat.geom")
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  ctrl <- spatstat.geom::owin(c(0, 5), c(0, 10))
  params <- .biv_test_params_bg_only(mu_0 = 5, mu_1 = 0)
  windowT <- c(0, 20)

  # Same seed, default ref area: explicit own-area ref must reproduce the
  # legacy stream exactly.
  set.seed(42)
  sim_default <- sim_etas_bivariate(
    params = params, windowT = windowT, windowS = win,
    state_spaces = list(control = win, treated = NULL),
    m0 = 2.5, beta_gr = 2.3
  )
  set.seed(42)
  sim_own_ref <- sim_etas_bivariate(
    params = params, windowT = windowT, windowS = win,
    state_spaces = list(control = win, treated = NULL),
    m0 = 2.5, beta_gr = 2.3,
    bg_ref_areas = list(control = spatstat.geom::area(win))
  )
  expect_identical(sim_default, sim_own_ref)

  # Referencing the density to the control region (half of the domain)
  # scales the expected count by 2: mu_0*T = 100 -> 200.
  set.seed(43)
  sim_ext <- sim_etas_bivariate(
    params = params, windowT = windowT, windowS = win,
    state_spaces = list(control = win, treated = NULL),
    m0 = 2.5, beta_gr = 2.3,
    bg_ref_areas = list(control = spatstat.geom::area(ctrl))
  )
  expect_gt(nrow(sim_ext), 150)
  expect_lt(nrow(sim_ext), 260)
  set.seed(43)
  sim_base <- sim_etas_bivariate(
    params = params, windowT = windowT, windowS = win,
    state_spaces = list(control = win, treated = NULL),
    m0 = 2.5, beta_gr = 2.3
  )
  expect_gt(nrow(sim_base), 60)
  expect_lt(nrow(sim_base), 140)
})
