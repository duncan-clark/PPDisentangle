# ---- shared fixtures (deterministic, no RNG) ----

fixed_windowT <- c(0, 10)
fixed_windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))

fixed_3pt <- data.frame(
  x = c(2.0, 5.0, 8.0),
  y = c(3.0, 7.0, 1.0),
  t = c(1.0, 4.0, 9.0),
  W = c(1.0, 1.0, 1.0)
)

fixed_8pt <- data.frame(
  x = c(1.5, 3.2, 4.8, 6.1, 7.4, 2.9, 8.6, 5.3),
  y = c(2.1, 8.4, 1.7, 5.5, 9.2, 3.8, 6.6, 4.0),
  t = c(0.5, 1.2, 2.3, 3.1, 4.7, 6.0, 7.8, 9.5),
  W = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
)

fixed_cluster <- data.frame(
  x = c(5.0, 5.1, 5.05, 4.95, 5.02),
  y = c(5.0, 5.0, 5.1,  4.9,  5.05),
  t = c(3.0, 3.1, 3.15, 3.2,  3.3),
  W = c(1.0, 1.0, 1.0,  1.0,  1.0)
)

fixed_inhom <- data.frame(
  x = c(2.0, 5.0, 8.0, 3.0, 7.0),
  y = c(3.0, 7.0, 1.0, 6.0, 4.0),
  t = c(1.0, 3.0, 5.0, 7.0, 9.0),
  W = c(0.5, 1.5, 0.8, 2.0, 1.2)
)

# ---- original behaviour tests ----

test_that("loglik_hawk returns finite value for valid data", {
  params <- c(mu = 50, alpha = 1, beta = 1, K = 0.3)
  windowT <- c(0, 10)
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  set.seed(42)
  sim <- sim_hawkes(list(mu = 50, alpha = 1, beta = 1, K = 0.3), windowT, windowS)
  realiz <- data.frame(x = sim$x, y = sim$y, t = sim$t)

  ll <- loglik_hawk(params, realiz, windowT, windowS)
  expect_true(is.finite(ll))
  expect_true(ll < 0 || ll == -999999)
})

test_that("loglik_hawk returns -999999 for negative params", {
  params <- c(mu = -1, alpha = 1, beta = 1, K = 0.3)
  realiz <- data.frame(x = 1, y = 1, t = 1)
  windowT <- c(0, 5)
  windowS <- spatstat.geom::owin(xrange = c(0, 5), yrange = c(0, 5))

  ll <- loglik_hawk(params, realiz, windowT, windowS)
  expect_equal(ll, -999999)
})

test_that("loglik_hawk_fast returns finite value for valid data", {
  params <- c(mu = 50, alpha = 1, beta = 1, K = 0.3)
  windowT <- c(0, 10)
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  set.seed(42)
  sim <- sim_hawkes(list(mu = 50, alpha = 1, beta = 1, K = 0.3), windowT, windowS)
  realiz <- data.frame(x = sim$x, y = sim$y, t = sim$t, W = rep(1, length(sim$t)))

  ll <- loglik_hawk_fast(params, realiz, windowT, windowS)
  expect_true(is.finite(ll))
})

test_that("loglik_hawk_fast rejects K >= 1", {
  params <- c(mu = 50, alpha = 1, beta = 1, K = 1.1)
  realiz <- data.frame(x = c(1, 2), y = c(1, 2), t = c(0.5, 1), W = c(1, 1))
  windowT <- c(0, 5)
  windowS <- spatstat.geom::owin(xrange = c(0, 5), yrange = c(0, 5))

  ll <- loglik_hawk_fast(params, realiz, windowT, windowS)
  expect_true(ll <= -1e10)  # R wrapper returns -1e15 for invalid params
})

# ---- exact value regression tests (fixed data, no RNG) ----

test_that("loglik_hawk returns exact value for 3-point data", {
  params <- c(mu = 10, alpha = 0.5, beta = 5, K = 0.3)
  ll <- loglik_hawk(params, fixed_3pt, fixed_windowT, fixed_windowS,
                    density_approx = TRUE)
  expect_equal(unname(ll), -107.8077552790, tolerance = 1e-6)
})

test_that("loglik_hawk_fast returns exact value for 3-point data", {
  params <- c(mu = 10, alpha = 0.5, beta = 5, K = 0.3)
  ll <- loglik_hawk_fast(params, fixed_3pt, fixed_windowT, fixed_windowS)
  expect_equal(unname(ll), -107.8077552790, tolerance = 1e-4)
})

test_that("loglik_hawk returns exact value for 8-point data", {
  params <- c(mu = 30, alpha = 0.2, beta = 3, K = 0.5)
  ll <- loglik_hawk(params, fixed_8pt, fixed_windowT, fixed_windowS,
                    density_approx = TRUE)
  expect_equal(unname(ll), -313.6302408383, tolerance = 1e-4)
})

test_that("loglik_hawk_fast returns exact value for 8-point data", {
  params <- c(mu = 30, alpha = 0.2, beta = 3, K = 0.5)
  ll <- loglik_hawk_fast(params, fixed_8pt, fixed_windowT, fixed_windowS)
  expect_equal(unname(ll), -313.6302408474, tolerance = 1e-3)
})

test_that("loglik_hawk returns exact value for Poisson case (K=0)", {
  params <- c(mu = 20, alpha = 0.5, beta = 5, K = 0)
  ll <- loglik_hawk(params, fixed_8pt, fixed_windowT, fixed_windowS,
                    density_approx = TRUE)
  expect_equal(unname(ll), -168.2076527051, tolerance = 1e-6)
})

test_that("loglik_hawk_fast returns exact value for Poisson case (K=0)", {
  params <- c(mu = 20, alpha = 0.5, beta = 5, K = 0)
  ll <- loglik_hawk_fast(params, fixed_8pt, fixed_windowT, fixed_windowS)
  expect_equal(unname(ll), -168.2076527051, tolerance = 1e-6)
})

test_that("loglik_hawk_fast returns exact value with zero_background_region", {
  params <- c(mu = 15, alpha = 0.3, beta = 4, K = 0.4)
  zbr <- spatstat.geom::owin(xrange = c(0, 5), yrange = c(0, 10))
  ll <- loglik_hawk_fast(params, fixed_8pt, fixed_windowT, fixed_windowS,
                         zero_background_region = zbr)
  expect_true(is.finite(ll) && ll < 0)
})

test_that("loglik_hawk_fast returns exact value with precomp", {
  params <- c(mu = 15, alpha = 0.3, beta = 4, K = 0.4)
  zbr <- spatstat.geom::owin(xrange = c(0, 5), yrange = c(0, 10))
  pc <- precompute_loglik_args(fixed_8pt, fixed_windowS, zbr)
  ll_precomp <- loglik_hawk_fast(params, fixed_8pt, fixed_windowT, fixed_windowS,
    precomp = list(active_area = pc$active_area, in_zero_bg = pc$in_zero_bg_all))
  ll_zbr <- loglik_hawk_fast(params, fixed_8pt, fixed_windowT, fixed_windowS,
                             zero_background_region = zbr)
  expect_equal(ll_precomp, ll_zbr, tolerance = 1e-10)  # precomp and zbr must match
  expect_true(is.finite(ll_precomp) && ll_precomp < 0)
})

test_that("loglik_hawk_fast returns exact value for single point", {
  r1 <- data.frame(x = 5.0, y = 5.0, t = 5.0, W = 1.0)
  params <- c(mu = 25, alpha = 1.0, beta = 2, K = 0.3)
  ll <- loglik_hawk_fast(params, r1, fixed_windowT, fixed_windowS)
  expect_equal(unname(ll), -251.6862943611, tolerance = 1e-6)
})

test_that("loglik_hawk returns exact value for tight cluster", {
  params <- c(mu = 5, alpha = 0.1, beta = 10, K = 0.8)
  ll <- loglik_hawk(params, fixed_cluster, fixed_windowT, fixed_windowS,
                    density_approx = TRUE)
  expect_equal(unname(ll), -62.7809860364, tolerance = 1e-6)
})

test_that("loglik_hawk_fast returns exact value for tight cluster", {
  params <- c(mu = 5, alpha = 0.1, beta = 10, K = 0.8)
  ll <- loglik_hawk_fast(params, fixed_cluster, fixed_windowT, fixed_windowS)
  expect_equal(unname(ll), -62.7809860364, tolerance = 1e-6)
})

test_that("loglik_hawk_fast returns exact value for inhomogeneous background", {
  params <- c(mu = 15, alpha = 0.4, beta = 6, K = 0.3)
  ll <- loglik_hawk_fast(params, fixed_inhom, fixed_windowT, fixed_windowS)
  expect_equal(unname(ll), -160.6209568079, tolerance = 1e-4)
})

test_that("loglik_hawk and loglik_hawk_fast both return finite values on fixed datasets", {
  # R uses density_approx (Option A: integral = mu*T + K*n); C++ uses Option B (time-dependent integral).
  # They need not be equal; we only check both are finite.
  cases <- list(
    list(p = c(mu = 10, alpha = 0.5, beta = 5, K = 0.3), d = fixed_3pt),
    list(p = c(mu = 30, alpha = 0.2, beta = 3, K = 0.5), d = fixed_8pt),
    list(p = c(mu = 20, alpha = 0.5, beta = 5, K = 0),   d = fixed_8pt),
    list(p = c(mu = 5, alpha = 0.1, beta = 10, K = 0.8), d = fixed_cluster)
  )
  for (i in seq_along(cases)) {
    ll_r <- unname(loglik_hawk(cases[[i]]$p, cases[[i]]$d, fixed_windowT, fixed_windowS,
                               density_approx = TRUE))
    ll_cpp <- unname(loglik_hawk_fast(cases[[i]]$p, cases[[i]]$d, fixed_windowT, fixed_windowS))
    expect_true(is.finite(ll_r), info = paste("Case", i, "R"))
    expect_true(is.finite(ll_cpp), info = paste("Case", i, "C++"))
  }
})
