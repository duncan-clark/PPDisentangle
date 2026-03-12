# ---- ETAS simulation tests ----

etas_sim_windowT <- c(0, 20)
etas_sim_windowS <- spatstat.geom::owin(xrange = c(0, 50), yrange = c(0, 50))
etas_sim_m0 <- 2.5

test_that("sim_etas returns all points inside the spatial/temporal window", {
  set.seed(101)
  params <- list(mu = 3, A = 0.2, alpha_m = 0.8, c = 0.1, p = 1.3,
                 D = 5.0, gamma = 0.5, q = 1.5)
  sim <- sim_etas(params, etas_sim_windowT, etas_sim_windowS,
                  m0 = etas_sim_m0, beta_gr = 2.3)
  n <- length(sim$t)
  if (n > 0) {
    expect_true(all(sim$t >= etas_sim_windowT[1]))
    expect_true(all(sim$t <= etas_sim_windowT[2]))
    expect_true(all(sim$x >= 0 & sim$x <= 50))
    expect_true(all(sim$y >= 0 & sim$y <= 50))
    expect_true(all(!is.na(sim$mag)))
    expect_true(all(sim$mag >= etas_sim_m0))
  }
})

test_that("sim_etas with A=0 gives only background events", {
  set.seed(102)
  params <- list(mu = 5, A = 0, alpha_m = 1, c = 0.1, p = 1.2,
                 D = 1, gamma = 0.5, q = 1.5)
  sim <- sim_etas(params, etas_sim_windowT, etas_sim_windowS,
                  m0 = etas_sim_m0, beta_gr = 2.3)
  if (length(sim$t) > 0) {
    expect_true(all(sim$background))
  }
})

test_that("sim_etas produces more events with larger A", {
  set.seed(103)
  params_low  <- list(mu = 3, A = 0.1, alpha_m = 0.5, c = 0.1, p = 1.3,
                      D = 5, gamma = 0.3, q = 1.5)
  params_high <- list(mu = 3, A = 0.6, alpha_m = 0.5, c = 0.1, p = 1.3,
                      D = 5, gamma = 0.3, q = 1.5)
  n_reps <- 10
  counts_low <- numeric(n_reps)
  counts_high <- numeric(n_reps)
  for (i in 1:n_reps) {
    counts_low[i]  <- length(sim_etas(params_low, etas_sim_windowT,
                                       etas_sim_windowS, m0 = etas_sim_m0,
                                       beta_gr = 2.3)$t)
    counts_high[i] <- length(sim_etas(params_high, etas_sim_windowT,
                                       etas_sim_windowS, m0 = etas_sim_m0,
                                       beta_gr = 2.3)$t)
  }
  expect_true(mean(counts_high) > mean(counts_low))
})

test_that("sim_etas offspring magnitudes follow GR distribution", {
  set.seed(104)
  beta_gr <- 2.3
  params <- list(mu = 5, A = 0.4, alpha_m = 0.8, c = 0.05, p = 1.3,
                 D = 10.0, gamma = 0.5, q = 1.5)
  sim <- sim_etas(params, etas_sim_windowT, etas_sim_windowS,
                  m0 = etas_sim_m0, beta_gr = beta_gr)
  offspring_mags <- sim$mag[!sim$background]
  if (length(offspring_mags) > 30) {
    ks <- ks.test(offspring_mags - etas_sim_m0, "pexp", rate = beta_gr)
    expect_true(ks$p.value > 0.01,
                label = "offspring magnitudes should follow GR distribution")
  }
})

test_that("sim_etas with mag_pool resamples observed magnitudes", {
  set.seed(105)
  pool <- c(3.0, 3.5, 4.0, 2.8, 3.2)
  params <- list(mu = 3, A = 0.3, alpha_m = 0.5, c = 0.1, p = 1.3,
                 D = 5, gamma = 0.3, q = 1.5)
  sim <- sim_etas(params, etas_sim_windowT, etas_sim_windowS,
                  m0 = etas_sim_m0, mag_pool = pool)
  offspring_mags <- sim$mag[!sim$background]
  if (length(offspring_mags) > 0) {
    expect_true(all(offspring_mags %in% pool))
  }
})

test_that("sim_etas with truncation produces no child beyond t_trunc of parent", {
  set.seed(106)
  t_trunc <- 5.0
  params <- list(mu = 2, A = 0.3, alpha_m = 0.5, c = 0.05, p = 1.2,
                 D = 5, gamma = 0.3, q = 1.5)
  for (rep in 1:5) {
    sim <- sim_etas(params, c(0, 100), etas_sim_windowS,
                    m0 = etas_sim_m0, beta_gr = 2.3, t_trunc = t_trunc)
    if (length(sim$t) > 1) {
      bg_times <- sim$t[sim$background]
      ch_times <- sim$t[!sim$background]
      if (length(ch_times) > 0 && length(bg_times) > 0) {
        for (ct in ch_times) {
          min_delay <- min(ct - bg_times[bg_times < ct])
          if (is.finite(min_delay)) {
            expect_true(min_delay <= t_trunc + 1e-10)
          }
        }
      }
    }
  }
})

test_that("sim_etas with background_realization works", {
  set.seed(107)
  bg <- list(x = c(10, 20, 30), y = c(10, 20, 30), t = c(1, 5, 10),
             mag = c(3.0, 3.5, 3.2))
  params <- list(mu = 0, A = 0.3, alpha_m = 0.5, c = 0.1, p = 1.3,
                 D = 5, gamma = 0.3, q = 1.5)
  sim <- sim_etas(params, etas_sim_windowT, etas_sim_windowS,
                  m0 = etas_sim_m0, beta_gr = 2.3,
                  background_realization = bg)
  expect_true(length(sim$t) >= 3)
  expect_true(all(bg$t %in% sim$t))
})

test_that("sim_etas_fast is an alias for sim_etas", {
  expect_identical(sim_etas_fast, sim_etas)
})

test_that("sim_etas_children_cpp returns correct types", {
  ch <- sim_etas_children_cpp(
    parent_x = 25, parent_y = 25, parent_t = 5, parent_mag = 3.5,
    A = 0.5, alpha_m = 1.0, cc = 0.1, p = 1.3,
    D = 5.0, gamma_par = 0.5, q = 1.5, m0 = 2.5, beta_gr = 2.3,
    t_min = 0, t_max = 100, x_min = 0, x_max = 50, y_min = 0, y_max = 50
  )
  expect_true(is.numeric(ch$x))
  expect_true(is.numeric(ch$y))
  expect_true(is.numeric(ch$t))
  expect_true(is.numeric(ch$mag))
  expect_equal(length(ch$x), length(ch$t))
  expect_equal(length(ch$x), length(ch$mag))
})
