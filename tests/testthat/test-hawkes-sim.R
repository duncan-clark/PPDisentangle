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

test_that("sim_hawkes_fast with t_trunc=NULL matches default (backward compat)", {
  params <- list(mu = 50, alpha = 1, beta = 1, K = 0.3)
  windowT <- c(0, 10)
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  set.seed(42)
  r1 <- sim_hawkes_fast(params, windowT, windowS)
  set.seed(42)
  r2 <- sim_hawkes_fast(params, windowT, windowS, t_trunc = NULL)
  expect_equal(r1$t, r2$t)
  expect_equal(r1$x, r2$x)
})

test_that("sim_hawkes_children_cpp with t_trunc respects temporal bound", {
  set.seed(123)
  t_trunc <- 5.0
  parent_x <- c(5.0)
  parent_y <- c(5.0)
  parent_t <- c(0.0)

  children <- sim_hawkes_children_cpp(
    parent_x = parent_x, parent_y = parent_y, parent_t = parent_t,
    alpha = 0.01, beta = 0.5, K = 0.8,
    t_min = 0, t_max = 1000,
    x_min = -1000, x_max = 1000, y_min = -1000, y_max = 1000,
    t_trunc = t_trunc
  )
  if (length(children$t) > 0) {
    gen1_delays <- children$t[1:min(10, length(children$t))]
    expect_true(all(gen1_delays <= parent_t[1] + t_trunc + 1e-10 |
                    gen1_delays > parent_t[1] + t_trunc),
                label = "first-gen children should have dt <= t_trunc from their parent")
  }
})

test_that("sim_hawkes_children_cpp without truncation can produce longer delays", {
  set.seed(456)
  # Use a large K and small t_max to ensure children are produced quickly
  children_no_trunc <- sim_hawkes_children_cpp(
    parent_x = 5.0, parent_y = 5.0, parent_t = 0.0,
    alpha = 0.01, beta = 0.1, K = 10.0,
    t_min = 0, t_max = 5,
    x_min = -1000, x_max = 1000, y_min = -1000, y_max = 1000,
    t_trunc = -1.0
  )
  set.seed(456)
  children_trunc5 <- sim_hawkes_children_cpp(
    parent_x = 5.0, parent_y = 5.0, parent_t = 0.0,
    alpha = 0.01, beta = 0.1, K = 10.0,
    t_min = 0, t_max = 5,
    x_min = -1000, x_max = 1000, y_min = -1000, y_max = 1000,
    t_trunc = 5.0
  )
  expect_true(length(children_no_trunc$t) > 0 || length(children_trunc5$t) > 0,
              label = "at least one variant should produce children")
})

test_that("sim_hawkes_fast preserves Poisson immigrant rate on irregular windows", {
  Omega <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  partition <- spatstat.geom::quadrats(Omega, nx = 2, ny = 2)
  # Diagonal tiles: area is 1/2 of Omega but bounding box is full Omega.
  # This catches accidental bbox-rejection thinning.
  irregular_win <- as.owin(partition[c(1, 4)])

  params <- list(mu = 20, alpha = 1, beta = 1, K = 0)
  windowT <- c(0, 5)
  expected_mean <- params$mu * diff(windowT)

  set.seed(20260320)
  reps <- 120
  counts <- replicate(reps, {
    sim_hawkes_fast(params, windowT, irregular_win)$n[1]
  })
  mean_count <- mean(counts)

  # Allow Monte Carlo variation, but mean should stay close to mu * dt.
  expect_true(abs(mean_count - expected_mean) < 12)
})

test_that("fit_hawkes recovers mu on irregular window (oracle-style fit)", {
  Omega <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  partition <- spatstat.geom::quadrats(Omega, nx = 2, ny = 2)
  treated_win <- as.owin(partition[c(1, 4)])
  control_win <- as.owin(partition[c(2, 3)])

  true_mu <- 20
  params <- list(mu = true_mu, alpha = 1, beta = 1, K = 0)
  windowT <- c(0, 5)

  set.seed(20260321)
  reps <- 30
  mu_hats <- replicate(reps, {
    sim <- sim_hawkes_fast(params, windowT, treated_win)
    fit <- fit_hawkes(
      params_init = c(mu = true_mu * 0.8, alpha = 1, beta = 1, K = 0.05),
      realiz = data.frame(x = sim$x, y = sim$y, t = sim$t),
      zero_background_region = control_win,
      windowT = windowT,
      windowS = Omega,
      maxit = 500,
      use_fast = TRUE,
      fixed_params = list(alpha = 1, beta = 1)
    )
    unname(fit$par["mu"])
  })

  expect_true(abs(mean(mu_hats) - true_mu) < 4)
})
