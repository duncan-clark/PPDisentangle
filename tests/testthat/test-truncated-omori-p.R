etas_omori_pts <- data.frame(
  x   = c(2.0, 5.0, 8.0, 3.0, 7.0),
  y   = c(3.0, 7.0, 1.0, 6.0, 4.0),
  t   = c(1.0, 3.0, 5.0, 7.0, 9.0),
  mag = c(3.0, 3.5, 2.8, 3.2, 3.0)
)
etas_omori_windowT <- c(0, 10)
etas_omori_windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))

etas_omori_params <- function(p) {
  c(mu = 5, A = 0.3, alpha_m = 0.4, c = 0.1, p = p,
    D = 1, gamma = 0, q = 1.6)
}

test_that("truncated Omori loglik is finite for p < 1 and p = 1", {
  skip_if_not_installed("spatstat.geom")
  for (p in c(0.5, 0.8, 1.0, 1.2)) {
    ll <- loglik_etas(
      etas_omori_params(p), etas_omori_pts, etas_omori_windowT, etas_omori_windowS,
      m0 = 2.5, t_trunc = 5, enforce_finite_trigger_moments = FALSE,
      enforce_alpha_subcritical = FALSE, beta_gr = 1.5
    )
    expect_true(is.finite(ll), info = paste("p =", p))
    expect_gt(ll, -1e10)
  }
})

test_that("truncated Omori loglik is continuous across p = 1", {
  skip_if_not_installed("spatstat.geom")
  ll_at <- function(p) {
    loglik_etas(
      etas_omori_params(p), etas_omori_pts, etas_omori_windowT, etas_omori_windowS,
      m0 = 2.5, t_trunc = 8, enforce_finite_trigger_moments = FALSE,
      enforce_alpha_subcritical = FALSE, beta_gr = 1.5
    )
  }
  ll0 <- ll_at(1)
  expect_true(is.finite(ll0))
  expect_lt(abs(ll_at(1 - 1e-6) - ll0), 1e-3)
  expect_lt(abs(ll_at(1 + 1e-6) - ll0), 1e-3)
})

test_that("truncated bivariate Omori accepts p < 1", {
  skip_if_not_installed("spatstat.geom")
  realiz <- data.frame(
    x = etas_omori_pts$x, y = etas_omori_pts$y,
    t = etas_omori_pts$t, mag = etas_omori_pts$mag,
    location_process = c("control", "treated", "control", "treated", "control")
  )
  params <- c(
    mu_0 = 0.8, mu_1 = 0.6,
    A_00 = 0.15, alpha_m_00 = 0.4,
    A_11 = 0.12, alpha_m_11 = 0.4,
    A_01 = 0.02, alpha_m_01 = 0.3,
    A_10 = 0.02, alpha_m_10 = 0.3,
    c = 0.2, p = 0.7, D = 1.0, gamma = 0, q = 1.6
  )
  ll <- loglik_etas_bivariate(
    params, realiz, windowT = etas_omori_windowT, windowS = etas_omori_windowS,
    m0 = 2.5, beta_gr = 1.5, max_branching_radius = 0.98, t_trunc = 6,
    enforce_finite_trigger_moments = FALSE
  )
  expect_true(is.finite(ll))
  expect_gt(ll, -1e10)
  ll_untrunc <- loglik_etas_bivariate(
    params, realiz, windowT = etas_omori_windowT, windowS = etas_omori_windowS,
    m0 = 2.5, beta_gr = 1.5, max_branching_radius = 0.98, t_trunc = NULL,
    enforce_finite_trigger_moments = FALSE
  )
  expect_true(ll_untrunc <= -1e10)
})

test_that("truncated Omori simulation with p < 1 stays inside t_trunc", {
  skip_if_not_installed("spatstat.geom")
  t_trunc <- 4
  params <- c(mu = 0, A = 0.8, alpha_m = 0.2, c = 0.05, p = 0.6,
              D = 0.5, gamma = 0, q = 2.0)
  history <- data.frame(x = 5, y = 5, t = -0.01, mag = 3.0)
  set.seed(20260812)
  sim <- sim_etas(
    params = params,
    windowT = c(0, 8),
    windowS = spatstat.geom::owin(c(0, 10), c(0, 10)),
    m0 = 2.5, beta_gr = 2.3,
    filtration = history,
    t_trunc = t_trunc
  )
  expect_true(is.list(sim))
  expect_true(is.numeric(sim$t))
  if (length(sim$t) > 0L) {
    expect_true(all(is.finite(sim$t)))
    expect_true(all(sim$t >= 0))
    expect_true(all(sim$t <= 8))
  }
})
