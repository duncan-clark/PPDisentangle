test_that("sim_hawkes produces events inside the window", {
  params <- list(mu = 50, alpha = 1, beta = 1, K = 0.3)
  windowT <- c(0, 10)
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  set.seed(42)
  result <- sim_hawkes(params, windowT, windowS, optimized = TRUE)

  expect_true(length(result$t) > 0)
  expect_true(all(result$t >= windowT[1] & result$t <= windowT[2]))
  expect_true(all(result$x >= 0 & result$x <= 10))
  expect_true(all(result$y >= 0 & result$y <= 10))
  expect_true(any(result$background))
})

test_that("sim_hawkes with K=0 is a Poisson process (no children)", {
  params <- list(mu = 100, alpha = 1, beta = 1, K = 0)
  windowT <- c(0, 5)
  windowS <- spatstat.geom::owin(xrange = c(0, 5), yrange = c(0, 5))
  set.seed(123)
  result <- sim_hawkes(params, windowT, windowS)

  expect_true(all(result$background))
})

test_that("sim_hawkes_fast produces valid output", {
  params <- list(mu = 50, alpha = 1, beta = 1, K = 0.3)
  windowT <- c(0, 10)
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  set.seed(42)
  result <- sim_hawkes_fast(params, windowT, windowS)

  expect_true(length(result$t) > 0)
  expect_true(all(result$t >= windowT[1] & result$t <= windowT[2]))
  expect_equal(length(result$x), length(result$t))
  expect_equal(length(result$y), length(result$t))
  expect_equal(length(result$background), length(result$t))
})
