test_that("univariate GR branching constraint uses supplied beta_gr", {
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  dat <- data.frame(
    x = c(2, 4, 6, 8),
    y = c(2, 4, 6, 8),
    t = c(1, 2, 3, 4),
    mag = c(2.6, 2.8, 3.0, 3.2),
    W = 1
  )
  stable <- c(
    mu = 2, A = 0.2, alpha_m = 0.5,
    c = 0.1, p = 2.1, D = 1, gamma = 0.5, q = 1.6
  )
  explosive <- stable
  explosive[["A"]] <- 0.9

  expect_true(is.finite(loglik_etas(
    stable, dat, c(0, 5), win, m0 = 2.5,
    beta_gr = 2.3, max_branching_ratio = 0.98
  )))
  expect_lte(loglik_etas(
    explosive, dat, c(0, 5), win, m0 = 2.5,
    beta_gr = 2.3, max_branching_ratio = 0.98
  ), -1e10)
})

test_that("fit_etas projects its result inside the configured GR margin", {
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  dat <- data.frame(
    x = c(2, 4, 6, 8),
    y = c(2, 4, 6, 8),
    t = c(1, 2, 3, 4),
    mag = c(2.6, 2.8, 3.0, 3.2),
    W = 1
  )
  init <- c(
    mu = 2, A = 1.2, alpha_m = 2.2,
    c = 0.1, p = 2.1, D = 1, gamma = 0.5, q = 1.6
  )
  fit <- fit_etas(
    init, dat, c(0, 5), win, m0 = 2.5,
    beta_gr = 2.3, max_branching_ratio = 0.98, maxit = 1
  )

  eta <- PPDisentangle:::.etas_univ_branching_ratio(fit$par, 2.3)
  expect_lt(eta, 0.98)
  expect_lt(fit$par[["alpha_m"]], 2.3)
  expect_equal(fit$branching_ratio, eta)
})

test_that("fit_etas soft barrier keeps NM starts interior and can polish A", {
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  dat <- data.frame(
    x = c(2, 4, 6, 8),
    y = c(2, 4, 6, 8),
    t = c(1, 2, 3, 4),
    mag = c(2.6, 2.8, 3.0, 3.2),
    W = 1
  )
  init <- c(
    mu = 2, A = 0.9, alpha_m = 0.5,
    c = 0.1, p = 2.1, D = 1, gamma = 0.5, q = 1.6
  )
  eta0 <- PPDisentangle:::.etas_univ_branching_ratio(init, 2.3)
  expect_gt(eta0, 0.98)

  fit <- fit_etas(
    init, dat, c(0, 5), win, m0 = 2.5,
    beta_gr = 2.3, max_branching_ratio = 0.98,
    maxit = 20, log_transform = TRUE,
    soft_branching_barrier = TRUE,
    polish_productivity = TRUE,
    interior_restart = TRUE,
    init_branching_margin = 0.9
  )

  expect_true(isTRUE(fit$soft_branching_barrier))
  expect_true(isTRUE(fit$log_transform))
  expect_lt(fit$branching_ratio, 0.98)
  expect_true(is.finite(fit$value))
  expect_gt(fit$value, -1e10)
})

test_that("bivariate likelihood and fit enforce the supplied GR radius", {
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  ctrl <- spatstat.geom::owin(c(0, 5), c(0, 10))
  treat <- spatstat.geom::owin(c(5, 10), c(0, 10))
  dat <- data.frame(
    x = c(2, 3, 7, 8),
    y = c(2, 7, 3, 8),
    t = c(1, 2, 3, 4),
    mag = c(2.6, 2.8, 3.0, 3.2),
    inferred_process = c("control", "control", "treated", "treated")
  )
  stable <- c(
    mu_0 = 1, mu_1 = 1,
    A_00 = 0.1, alpha_m_00 = 0.5,
    A_11 = 0.1, alpha_m_11 = 0.5,
    A_01 = 0.01, alpha_m_01 = 0.5,
    A_10 = 0.01, alpha_m_10 = 0.5,
    c = 0.1, p = 2.1, D = 1, gamma = 0.5, q = 1.6
  )
  explosive <- stable
  explosive[["A_00"]] <- 1.2

  expect_true(is.finite(loglik_etas_bivariate(
    stable, dat, c(0, 5), win, m0 = 2.5,
    control_state_space = ctrl, treated_state_space = treat,
    beta_gr = 2.3, max_branching_radius = 0.98
  )))
  expect_lte(loglik_etas_bivariate(
    explosive, dat, c(0, 5), win, m0 = 2.5,
    control_state_space = ctrl, treated_state_space = treat,
    beta_gr = 2.3, max_branching_radius = 0.98
  ), -1e10)

  fit <- fit_etas_bivariate(
    explosive, dat, c(0, 5), win, m0 = 2.5,
    control_state_space = ctrl, treated_state_space = treat,
    beta_gr = 2.3, max_branching_radius = 0.98, maxit = 1
  )
  rho <- PPDisentangle:::.etas_biv_spectral_radius(fit$par, 2.3)
  expect_lt(rho, 0.98)
  expect_equal(fit$branching_radius, rho)
})

test_that("bivariate simulator supports matched inhomogeneous backgrounds", {
  set.seed(409)
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  params <- c(
    mu_0 = 100, mu_1 = 100,
    A_00 = 0, alpha_m_00 = 0.5,
    A_11 = 0, alpha_m_11 = 0.5,
    A_01 = 0, alpha_m_01 = 0.5,
    A_10 = 0, alpha_m_10 = 0.5,
    c = 0.1, p = 2.1, D = 1, gamma = 0.5, q = 1.6
  )
  sim <- sim_etas_bivariate(
    params, c(0, 3), win,
    state_spaces = list(control = win, treated = win),
    m0 = 2.5, beta_gr = 2.3,
    covariate_lookup = list(
      control = function(x, y) x / 5,
      treated = function(x, y) (10 - x) / 5
    )
  )

  expect_gt(sum(sim$process_id == 0L), 100)
  expect_gt(sum(sim$process_id == 1L), 100)
  expect_gt(mean(sim$x[sim$process_id == 0L]), 6)
  expect_lt(mean(sim$x[sim$process_id == 1L]), 4)
})
