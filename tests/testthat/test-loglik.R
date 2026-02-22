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
  expect_equal(ll, -Inf)
})
