test_that("bivariate ETAS kills offspring outside irregular support", {
  skip_if_not_installed("spatstat.geom")

  horizontal <- spatstat.geom::owin(c(0, 2), c(0, 1))
  vertical <- spatstat.geom::owin(c(0, 1), c(0, 2))
  support <- spatstat.geom::union.owin(horizontal, vertical)
  params <- c(
    mu_0 = 0, mu_1 = 0,
    A_00 = 0.5, alpha_m_00 = 0.2,
    A_11 = 0, alpha_m_11 = 0.2,
    A_01 = 0, alpha_m_01 = 0.2,
    A_10 = 0, alpha_m_10 = 0.2,
    c = 0.01, p = 2.2, D = 0.8, gamma = 0, q = 2.2
  )
  history <- data.frame(
    x = rep(0.5, 100), y = rep(0.5, 100),
    t = rep(-0.001, 100), mag = rep(2.5, 100),
    process_id = rep(0L, 100)
  )

  set.seed(20260812)
  sim <- sim_etas_bivariate(
    params = params,
    windowT = c(0, 1),
    windowS = support,
    state_spaces = list(control = support, treated = NULL),
    m0 = 2.5,
    beta_gr = 2.3,
    filtration = history,
    t_trunc = 1
  )

  expect_gt(nrow(sim), 0L)
  expect_true(all(spatstat.geom::inside.owin(sim$x, sim$y, support)))
})

test_that("bivariate post-treatment simulation conditions only on supported history", {
  skip_if_not_installed("spatstat.geom")

  horizontal <- spatstat.geom::owin(c(0, 2), c(0, 1))
  vertical <- spatstat.geom::owin(c(0, 1), c(0, 2))
  support <- spatstat.geom::union.owin(horizontal, vertical)
  params <- c(
    mu_0 = 0, mu_1 = 0,
    A_00 = 0.5, alpha_m_00 = 0.2,
    A_11 = 0, alpha_m_11 = 0.2,
    A_01 = 0, alpha_m_01 = 0.2,
    A_10 = 0, alpha_m_10 = 0.2,
    c = 0.01, p = 2.2, D = 0.2, gamma = 0, q = 2.2
  )
  inside_history <- data.frame(
    x = rep(0.5, 100), y = rep(0.5, 100),
    t = rep(-0.001, 100), mag = rep(2.5, 100),
    process_id = rep(0L, 100)
  )
  outside_history <- transform(inside_history, x = 1.5, y = 1.5)
  common <- list(
    params = params,
    windowT = c(0, 1),
    windowS = support,
    state_spaces = list(control = support, treated = NULL),
    m0 = 2.5,
    beta_gr = 2.3,
    t_trunc = 1
  )

  set.seed(20260812)
  conditioned <- do.call(
    sim_etas_bivariate, c(common, list(filtration = inside_history))
  )
  set.seed(20260812)
  unconditioned <- do.call(
    sim_etas_bivariate, c(common, list(filtration = NULL))
  )
  set.seed(20260812)
  unsupported <- do.call(
    sim_etas_bivariate, c(common, list(filtration = outside_history))
  )

  expect_gt(nrow(conditioned), 0L)
  expect_equal(nrow(unconditioned), 0L)
  expect_equal(nrow(unsupported), 0L)
  expect_true(all(conditioned$t >= 0))
})

test_that("univariate ETAS kills unsupported history before branching", {
  horizontal <- spatstat.geom::owin(c(0, 2), c(0, 1))
  vertical <- spatstat.geom::owin(c(0, 1), c(1, 2))
  support <- spatstat.geom::union.owin(horizontal, vertical)
  params <- list(
    mu = 0, A = 0.5, alpha_m = 0.2,
    c = 0.01, p = 2.2, D = 0.2, gamma = 0, q = 2.2
  )
  inside_history <- data.frame(
    x = rep(0.5, 100), y = rep(0.5, 100),
    t = rep(-0.001, 100), mag = rep(2.5, 100)
  )
  outside_history <- transform(inside_history, x = 1.5, y = 1.5)
  common <- list(
    params = params, windowT = c(0, 1), windowS = support,
    m0 = 2.5, beta_gr = 2.3, t_trunc = 1
  )

  set.seed(20260812)
  conditioned <- do.call(sim_etas, c(common, list(filtration = inside_history)))
  set.seed(20260812)
  unsupported <- do.call(sim_etas, c(common, list(filtration = outside_history)))

  expect_gt(length(conditioned$t), 0L)
  expect_equal(length(unsupported$t), 0L)
  expect_true(all(spatstat.geom::inside.owin(
    conditioned$x, conditioned$y, support
  )))
})

test_that("partitioned ETAS generators use tile union during branching", {
  horizontal <- spatstat.geom::owin(c(0, 2), c(0, 1))
  vertical <- spatstat.geom::owin(c(0, 1), c(1, 2))
  support <- spatstat.geom::union.owin(horizontal, vertical)
  omega <- spatstat.geom::owin(c(0, 2), c(0, 2))
  partition <- spatstat.geom::tess(
    tiles = list(lower = horizontal, upper_left = vertical),
    window = omega
  )
  partition_processes <- c(lower = "control", upper_left = "control")
  outside_history <- data.frame(
    x = rep(1.5, 100), y = rep(1.5, 100),
    t = rep(-0.001, 100), mag = rep(2.5, 100),
    location_process = "control", process_id = 0L
  )
  univ_params <- list(
    control = list(
      mu = 0, A = 0.5, alpha_m = 0.2,
      c = 0.01, p = 2.2, D = 0.2, gamma = 0, q = 2.2
    )
  )
  biv_params <- c(
    mu_0 = 0, mu_1 = 0,
    A_00 = 0.5, alpha_m_00 = 0.2,
    A_11 = 0, alpha_m_11 = 0.2,
    A_01 = 0, alpha_m_01 = 0.2,
    A_10 = 0, alpha_m_10 = 0.2,
    c = 0.01, p = 2.2, D = 0.2, gamma = 0, q = 2.2
  )

  set.seed(20260812)
  univ <- generate_inhomogeneous_etas(
    Omega = omega, partition = partition, time_window = c(0, 1),
    partition_processes = partition_processes,
    etas_params = univ_params, m0 = 2.5, beta_gr = 2.3,
    state_spaces = list(control = support),
    filtration = outside_history, t_trunc = 1
  )
  set.seed(20260812)
  biv <- generate_inhomogeneous_etas_bivariate(
    Omega = omega, partition = partition, time_window = c(0, 1),
    partition_processes = partition_processes,
    etas_bivariate_params = biv_params, m0 = 2.5, beta_gr = 2.3,
    state_spaces = list(control = support, treated = NULL),
    filtration = outside_history, t_trunc = 1
  )

  expect_equal(nrow(univ), 0L)
  expect_equal(nrow(biv), 0L)
})
