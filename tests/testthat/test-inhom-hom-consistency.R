cons_windowT <- c(0, 10)
cons_windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))

cons_hawkes <- data.frame(
  x = c(1.0, 2.5, 4.2, 7.8, 8.5),
  y = c(1.5, 3.0, 6.5, 7.2, 2.4),
  t = c(0.4, 1.8, 3.6, 6.9, 9.3),
  W = rep(1.0, 5)
)

cons_etas <- data.frame(
  x = c(1.0, 2.5, 4.2, 7.8, 8.5),
  y = c(1.5, 3.0, 6.5, 7.2, 2.4),
  t = c(0.4, 1.8, 3.6, 6.9, 9.3),
  mag = c(3.0, 2.8, 3.2, 3.4, 2.9),
  W = rep(1.0, 5)
)

test_that("hawkes inhomogeneous path matches homogeneous when W is constant 1", {
  params <- c(mu = 12, alpha = 0.25, beta = 3.5, K = 0.35)
  ll_with_W <- loglik_hawk_fast(params, cons_hawkes, cons_windowT, cons_windowS)
  ll_without_W <- loglik_hawk_fast(params, cons_hawkes[, c("x", "y", "t")], cons_windowT, cons_windowS)
  expect_equal(ll_with_W, ll_without_W, tolerance = 1e-10)
})

test_that("etas inhomogeneous path matches homogeneous when W is constant 1", {
  params <- c(mu = 4, A = 0.3, alpha_m = 1.1, c = 0.08, p = 1.25, D = 1.0, gamma = 0.5, q = 1.6)
  ll_with_W <- loglik_etas(params, cons_etas, cons_windowT, cons_windowS, m0 = 2.5)
  ll_without_W <- loglik_etas(params, cons_etas[, c("x", "y", "t", "mag")], cons_windowT, cons_windowS, m0 = 2.5)
  expect_equal(ll_with_W, ll_without_W, tolerance = 1e-10)
})

test_that("hawkes filtration cpp equals inhom cpp with no pre-history", {
  params <- c(mu = 12, alpha = 0.25, beta = 3.5, K = 0.35)
  shifted_t <- cons_hawkes$t - cons_windowT[1]
  t_max <- cons_windowT[2] - cons_windowT[1]
  areaS <- spatstat.geom::area(cons_windowS)

  ll_inhom <- hawkes_loglik_inhom_cpp(
    t = shifted_t,
    x = cons_hawkes$x,
    y = cons_hawkes$y,
    W_val = cons_hawkes$W,
    mu = params["mu"],
    alpha = params["alpha"],
    beta = params["beta"],
    K = params["K"],
    areaS = areaS,
    t_max = t_max,
    t_trunc = -1
  )

  ll_filtration <- hawkes_loglik_inhom_filtration_cpp(
    post_t = cons_hawkes$t,
    post_x = cons_hawkes$x,
    post_y = cons_hawkes$y,
    W_val = cons_hawkes$W,
    parent_t = cons_hawkes$t,
    parent_x = cons_hawkes$x,
    parent_y = cons_hawkes$y,
    mu = params["mu"],
    alpha = params["alpha"],
    beta = params["beta"],
    K = params["K"],
    areaS = areaS,
    t_start = cons_windowT[1],
    t_end = cons_windowT[2],
    adjust_factor = 1,
    t_trunc = -1
  )

  expect_equal(ll_filtration, ll_inhom, tolerance = 1e-10)
})

test_that("hawkes filtration cpp matches active-area formulation under zero-background mask", {
  params <- c(mu = 9, alpha = 0.18, beta = 2.4, K = 0.42)
  zbr <- spatstat.geom::owin(xrange = c(0, 5), yrange = c(0, 10))
  total_area <- spatstat.geom::area(cons_windowS)
  active_area <- total_area - spatstat.geom::area(zbr)
  in_zero <- spatstat.geom::inside.owin(cons_hawkes$x, cons_hawkes$y, zbr)
  W_masked <- cons_hawkes$W
  W_masked[in_zero] <- 0
  shifted_t <- cons_hawkes$t - cons_windowT[1]

  ll_active <- hawkes_loglik_inhom_cpp(
    t = shifted_t,
    x = cons_hawkes$x,
    y = cons_hawkes$y,
    W_val = W_masked,
    mu = params["mu"],
    alpha = params["alpha"],
    beta = params["beta"],
    K = params["K"],
    areaS = active_area,
    t_max = cons_windowT[2] - cons_windowT[1],
    t_trunc = -1
  )

  ll_filtration <- hawkes_loglik_inhom_filtration_cpp(
    post_t = cons_hawkes$t,
    post_x = cons_hawkes$x,
    post_y = cons_hawkes$y,
    W_val = W_masked,
    parent_t = cons_hawkes$t,
    parent_x = cons_hawkes$x,
    parent_y = cons_hawkes$y,
    mu = params["mu"],
    alpha = params["alpha"],
    beta = params["beta"],
    K = params["K"],
    areaS = active_area,
    t_start = cons_windowT[1],
    t_end = cons_windowT[2],
    adjust_factor = 1,
    t_trunc = -1
  )

  expect_equal(ll_filtration, ll_active, tolerance = 1e-10)
})
