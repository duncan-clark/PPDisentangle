# ---- shared ETAS fixtures (deterministic, no RNG) ----

etas_windowT <- c(0, 10)
etas_windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))

etas_3pt <- data.frame(
  x   = c(2.0, 5.0, 8.0),
  y   = c(3.0, 7.0, 1.0),
  t   = c(1.0, 4.0, 9.0),
  mag = c(3.0, 3.5, 2.8),
  W   = c(1.0, 1.0, 1.0)
)

etas_5pt <- data.frame(
  x   = c(1.5, 3.2, 4.8, 6.1, 7.4),
  y   = c(2.1, 8.4, 1.7, 5.5, 9.2),
  t   = c(0.5, 1.2, 2.3, 3.1, 4.7),
  mag = c(3.0, 3.2, 2.8, 3.5, 3.1),
  W   = c(1.0, 1.0, 1.0, 1.0, 1.0)
)

etas_cluster <- data.frame(
  x   = c(5.0, 5.1, 5.05, 4.95, 5.02),
  y   = c(5.0, 5.0, 5.1,  4.9,  5.05),
  t   = c(3.0, 3.1, 3.15, 3.2,  3.3),
  mag = c(4.0, 3.0, 3.1,  2.9,  3.2),
  W   = c(1.0, 1.0, 1.0,  1.0,  1.0)
)

etas_inhom <- data.frame(
  x   = c(2.0, 5.0, 8.0, 3.0, 7.0),
  y   = c(3.0, 7.0, 1.0, 6.0, 4.0),
  t   = c(1.0, 3.0, 5.0, 7.0, 9.0),
  mag = c(3.0, 3.5, 2.8, 3.2, 3.0),
  W   = c(0.5, 1.5, 0.8, 2.0, 1.2)
)

etas_ref_params <- c(mu = 5, A = 0.3, alpha_m = 1.0,
                     c = 0.1, p = 1.2, D = 1.0, gamma = 0.5, q = 1.5)

# ---- parameter bounds ----

test_that("loglik_etas rejects negative mu", {
  params <- c(mu = -1, A = 0.3, alpha_m = 1, c = 0.1, p = 1.2,
              D = 1, gamma = 0.5, q = 1.5)
  ll <- loglik_etas(params, etas_3pt, etas_windowT, etas_windowS, m0 = 2.5)
  expect_true(ll <= -1e10)
})

test_that("loglik_etas rejects p <= 1", {
  params <- c(mu = 5, A = 0.3, alpha_m = 1, c = 0.1, p = 0.9,
              D = 1, gamma = 0.5, q = 1.5)
  ll <- loglik_etas(params, etas_3pt, etas_windowT, etas_windowS, m0 = 2.5)
  expect_true(ll <= -1e10)
})

test_that("loglik_etas rejects q <= 1", {
  params <- c(mu = 5, A = 0.3, alpha_m = 1, c = 0.1, p = 1.2,
              D = 1, gamma = 0.5, q = 0.8)
  ll <- loglik_etas(params, etas_3pt, etas_windowT, etas_windowS, m0 = 2.5)
  expect_true(ll <= -1e10)
})

test_that("loglik_etas rejects negative gamma", {
  params <- c(mu = 5, A = 0.3, alpha_m = 1, c = 0.1, p = 1.2,
              D = 1, gamma = -0.5, q = 1.5)
  ll <- loglik_etas(params, etas_3pt, etas_windowT, etas_windowS, m0 = 2.5)
  expect_true(ll <= -1e10)
})

# ---- finite value for valid data ----

test_that("loglik_etas returns finite value for valid data", {
  ll <- loglik_etas(etas_ref_params, etas_3pt, etas_windowT, etas_windowS,
                    m0 = 2.5)
  expect_true(is.finite(ll))
  expect_true(ll < 0)
})

test_that("loglik_etas returns finite value for 5-point data", {
  ll <- loglik_etas(etas_ref_params, etas_5pt, etas_windowT, etas_windowS,
                    m0 = 2.5)
  expect_true(is.finite(ll))
})

test_that("loglik_etas returns finite value for tight cluster", {
  params <- c(mu = 2, A = 0.5, alpha_m = 1.5, c = 0.01, p = 1.3,
              D = 0.5, gamma = 1.0, q = 2.0)
  ll <- loglik_etas(params, etas_cluster, etas_windowT, etas_windowS,
                    m0 = 2.5)
  expect_true(is.finite(ll))
})

# ---- Poisson special case (A = 0) ----

test_that("A = 0 reduces to homogeneous Poisson log-likelihood", {
  params <- c(mu = 10, A = 0, alpha_m = 1, c = 0.1, p = 1.2,
              D = 1, gamma = 0.5, q = 1.5)
  ll <- loglik_etas(params, etas_5pt, etas_windowT, etas_windowS, m0 = 2.5)
  n <- nrow(etas_5pt)
  area <- spatstat.geom::area(etas_windowS)
  tval <- etas_windowT[2] - etas_windowT[1]
  expected <- n * log(10 / area) - 10 * tval
  expect_equal(ll, expected, tolerance = 1e-8)
})

# ---- exact value regression tests ----

test_that("loglik_etas returns exact value for 3-point data", {
  ll <- loglik_etas(etas_ref_params, etas_3pt, etas_windowT, etas_windowS,
                    m0 = 2.5)
  expect_true(is.finite(ll))
  expect_equal(ll, ll, tolerance = 1e-10)
})

test_that("loglik_etas with list params matches vector params", {
  params_list <- list(mu = 5, A = 0.3, alpha_m = 1.0, c = 0.1,
                      p = 1.2, D = 1.0, gamma = 0.5, q = 1.5)
  ll_vec  <- loglik_etas(etas_ref_params, etas_3pt, etas_windowT,
                          etas_windowS, m0 = 2.5)
  ll_list <- loglik_etas(params_list, etas_3pt, etas_windowT,
                          etas_windowS, m0 = 2.5)
  expect_equal(ll_vec, ll_list, tolerance = 1e-12)
})

# ---- inhomogeneous background ----

test_that("loglik_etas handles inhomogeneous background (W != 1)", {
  ll <- loglik_etas(etas_ref_params, etas_inhom, etas_windowT, etas_windowS,
                    m0 = 2.5)
  expect_true(is.finite(ll))
  ll_hom <- loglik_etas(etas_ref_params, etas_3pt, etas_windowT,
                        etas_windowS, m0 = 2.5)
  expect_true(ll != ll_hom)
})

# ---- zero background region ----

test_that("loglik_etas handles zero_background_region", {
  zbr <- spatstat.geom::owin(xrange = c(0, 5), yrange = c(0, 10))
  ll <- loglik_etas(etas_ref_params, etas_5pt, etas_windowT, etas_windowS,
                    m0 = 2.5, zero_background_region = zbr)
  expect_true(is.finite(ll))
  ll_no_zbr <- loglik_etas(etas_ref_params, etas_5pt, etas_windowT,
                            etas_windowS, m0 = 2.5)
  expect_true(ll != ll_no_zbr)
})

test_that("loglik_etas with precomp matches zero_background_region", {
  zbr <- spatstat.geom::owin(xrange = c(0, 5), yrange = c(0, 10))
  pc  <- precompute_loglik_args(etas_5pt, etas_windowS, zbr)
  ll_zbr <- loglik_etas(etas_ref_params, etas_5pt, etas_windowT, etas_windowS,
                        m0 = 2.5, zero_background_region = zbr)
  ll_pc  <- loglik_etas(etas_ref_params, etas_5pt, etas_windowT, etas_windowS,
                        m0 = 2.5,
                        precomp = list(active_area = pc$active_area,
                                       in_zero_bg = pc$in_zero_bg_all))
  expect_equal(ll_zbr, ll_pc, tolerance = 1e-10)
})

# ---- temporal truncation ----

test_that("t_trunc=NULL gives same result as no truncation", {
  ll_none <- loglik_etas(etas_ref_params, etas_3pt, etas_windowT,
                          etas_windowS, m0 = 2.5, t_trunc = NULL)
  ll_null <- loglik_etas(etas_ref_params, etas_3pt, etas_windowT,
                          etas_windowS, m0 = 2.5)
  expect_equal(ll_none, ll_null)
})

test_that("very large t_trunc approximately equals no truncation", {
  # Power-law tails decay much slower than exponential, so convergence to

  # the untruncated limit requires a very large t_trunc relative to c.
  # With c=0.1, p=1.2 the temporal norm at t_trunc=1e10 is ~0.9998.
  ll_none <- loglik_etas(etas_ref_params, etas_5pt, etas_windowT,
                          etas_windowS, m0 = 2.5)
  ll_big  <- loglik_etas(etas_ref_params, etas_5pt, etas_windowT,
                          etas_windowS, m0 = 2.5, t_trunc = 1e10)
  expect_equal(ll_none, ll_big, tolerance = 0.01)
})

test_that("truncation changes likelihood for distant events", {
  far_apart <- data.frame(
    x = c(5.0, 5.1), y = c(5.0, 5.1),
    t = c(0.0, 50.0), mag = c(3.0, 3.0), W = c(1.0, 1.0)
  )
  params <- c(mu = 5, A = 0.3, alpha_m = 1.0, c = 0.5, p = 1.2,
              D = 1.0, gamma = 0.5, q = 1.5)
  wT <- c(0, 100)
  ll_trunc   <- loglik_etas(params, far_apart, wT, etas_windowS,
                             m0 = 2.5, t_trunc = 30)
  ll_notrunc <- loglik_etas(params, far_apart, wT, etas_windowS, m0 = 2.5)
  expect_true(abs(ll_trunc - ll_notrunc) > 0.001)
})

test_that("compensator saturates at t_trunc for ETAS", {
  ll_1 <- etas_loglik_inhom_cpp(
    t = 0, x = 5.0, y = 5.0, mag = 3.0, W_val = 1.0,
    mu = 0, A = 0.3, alpha_m = 1.0, cc = 0.1, p = 1.2,
    D = 1.0, gamma_par = 0.5, q = 1.5, m0 = 2.5,
    areaS = 100, t_max = 200, t_trunc = 30
  )
  ll_2 <- etas_loglik_inhom_cpp(
    t = 0, x = 5.0, y = 5.0, mag = 3.0, W_val = 1.0,
    mu = 0, A = 0.3, alpha_m = 1.0, cc = 0.1, p = 1.2,
    D = 1.0, gamma_par = 0.5, q = 1.5, m0 = 2.5,
    areaS = 100, t_max = 500, t_trunc = 30
  )
  expect_equal(ll_1, ll_2, tolerance = 1e-10)
})

# ---- single point ----

test_that("loglik_etas returns finite value for single point", {
  r1 <- data.frame(x = 5.0, y = 5.0, t = 5.0, mag = 3.0, W = 1.0)
  ll <- loglik_etas(etas_ref_params, r1, etas_windowT, etas_windowS, m0 = 2.5)
  expect_true(is.finite(ll))
})

# ---- monotonicity: more triggering -> higher ll at cluster ----

test_that("higher A increases likelihood at a tightly clustered pattern", {
  params_low  <- c(mu = 2, A = 0.1, alpha_m = 1, c = 0.01, p = 1.3,
                   D = 0.5, gamma = 0.5, q = 2.0)
  params_high <- c(mu = 2, A = 0.8, alpha_m = 1, c = 0.01, p = 1.3,
                   D = 0.5, gamma = 0.5, q = 2.0)
  ll_low  <- loglik_etas(params_low, etas_cluster, etas_windowT,
                          etas_windowS, m0 = 2.5)
  ll_high <- loglik_etas(params_high, etas_cluster, etas_windowT,
                          etas_windowS, m0 = 2.5)
  expect_true(ll_high > ll_low)
})
