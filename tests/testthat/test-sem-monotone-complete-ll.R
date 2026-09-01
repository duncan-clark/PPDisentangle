# SEM with require_monotone_complete_ll must not walk downhill.

pkg_root <- NULL
for (cand in c(".", "..", "../..", file.path("..", "..", ".."))) {
  if (file.exists(file.path(cand, "DESCRIPTION"))) {
    pkg_root <- normalizePath(cand, winslash = "/", mustWork = FALSE)
    break
  }
}
if (!is.null(pkg_root) && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_root, quiet = TRUE, export_all = TRUE, helpers = FALSE,
                    attach_testthat = FALSE)
} else if (!exists("em_style_labelling", mode = "function")) {
  library(PPDisentangle)
}

test_that("monotone SEM does not decrease complete-data likelihood", {
  skip_if_not_installed("spatstat.geom")
  set.seed(7)
  n_pre <- 30L
  n_post <- 80L
  x <- c(runif(n_pre, 0, 5), runif(n_post, 0, 10))
  y <- runif(n_pre + n_post, 0, 10)
  t <- c(sort(runif(n_pre, -8, -0.01)), sort(runif(n_post, 0, 25)))
  mag <- 2.5 + rexp(n_pre + n_post, 1.5)
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  partition <- spatstat.geom::tess(xgrid = c(0, 5, 10), ygrid = c(0, 10))
  partition_processes <- c("control", "treated")
  state_spaces <- list(
    control = spatstat.geom::as.owin(partition[c(TRUE, FALSE)]),
    treated = spatstat.geom::as.owin(partition[c(FALSE, TRUE)])
  )
  pp <- data.frame(
    x = x, y = y, t = t, mag = mag,
    location_process = ifelse(x < 5 | t < 0, "control", "treated"),
    stringsAsFactors = FALSE
  )
  pp$inferred_process <- pp$location_process
  init <- c(
    mu_0 = 0.9, mu_1 = 0.7,
    A_00 = 0.15, alpha_m_00 = 0.4,
    A_11 = 0.18, alpha_m_11 = 0.45,
    A_01 = 0.01, alpha_m_01 = 0.3,
    A_10 = 0.01, alpha_m_10 = 0.3,
    c = 0.4, p = 2.05, D = 1.2, gamma = 0, q = 1.55
  )
  res <- em_style_labelling(
    pp_data = pp,
    partition = partition,
    partition_processes = partition_processes,
    statespace = win,
    state_spaces = state_spaces,
    time_window = c(0, 25),
    treatment_time = 0,
    hawkes_params_control = list(mu = 0.9, A = 0.15, alpha_m = 0.4, c = 0.4,
                                 p = 2.05, D = 1.2, gamma = 0, q = 1.55),
    hawkes_params_treated = list(mu = 0.7, A = 0.18, alpha_m = 0.45, c = 0.4,
                                 p = 2.05, D = 1.2, gamma = 0, q = 1.55),
    update_control_params = TRUE,
    param_update_cadence = 1L,
    proposal_update_cadence = 1,
    update_starting_data = TRUE,
    include_starting_data = TRUE,
    include_starting_first_n = 6L,
    optim_method = "sample_weighted",
    selection_temperature = 5,
    metric_name = "post_likelihood",
    require_monotone_complete_ll = TRUE,
    iter = 6L, n_props = 8L, change_factor = 0.35,
    stagnation_trigger_every = 50,
    proposal_method = "simulation",
    fixed_params = list(gamma = 0),
    verbose = FALSE,
    model_type = "etas_bivariate",
    outer_maxit_biv = 60L,
    m0 = 2.5, beta_gr = 1.5,
    etas_bivariate_params = init,
    hard_subcritical = TRUE,
    max_branching_radius = 0.98,
    alpha_beta_gap_min = 1e-4,
    t_trunc = 5
  )
  expect_true(!is.null(res$complete_ll_trace))
  expect_equal(res$max_metric_flips[[1L]], 0)
  ll <- as.numeric(res$complete_ll_trace)
  finite_ll <- ll[is.finite(ll)]
  expect_true(length(finite_ll) >= 1L)
  expect_true(all(diff(finite_ll) >= -1e-6))
  expect_true(is.finite(res$n_downhill_rejected))
  last_fit <- res$fits[[length(res$fits)]]
  if (!is.null(last_fit$complete_ll) && is.finite(last_fit$complete_ll)) {
    expect_equal(
      as.numeric(last_fit$complete_ll),
      max(finite_ll),
      tolerance = 1e-6
    )
  }
})

test_that("single_flip_from_iter accepts at most one flip after the switch", {
  skip_if_not_installed("spatstat.geom")
  set.seed(11)
  n_pre <- 20L
  n_post <- 24L
  x <- c(runif(n_pre, 0, 5), runif(n_post, 0, 10))
  y <- runif(n_pre + n_post, 0, 10)
  t <- c(sort(runif(n_pre, -6, -0.01)), sort(runif(n_post, 0, 18)))
  mag <- 2.5 + rexp(n_pre + n_post, 1.5)
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  partition <- spatstat.geom::tess(xgrid = c(0, 5, 10), ygrid = c(0, 10))
  partition_processes <- c("control", "treated")
  state_spaces <- list(
    control = spatstat.geom::as.owin(partition[c(TRUE, FALSE)]),
    treated = spatstat.geom::as.owin(partition[c(FALSE, TRUE)])
  )
  pp <- data.frame(
    x = x, y = y, t = t, mag = mag,
    location_process = ifelse(x < 5 | t < 0, "control", "treated"),
    stringsAsFactors = FALSE
  )
  pp$inferred_process <- pp$location_process
  init <- c(
    mu_0 = 0.9, mu_1 = 0.7,
    A_00 = 0.15, alpha_m_00 = 0.4,
    A_11 = 0.18, alpha_m_11 = 0.45,
    A_01 = 0.01, alpha_m_01 = 0.3,
    A_10 = 0.01, alpha_m_10 = 0.3,
    c = 0.4, p = 2.05, D = 1.2, gamma = 0, q = 1.55
  )
  res <- em_style_labelling(
    pp_data = pp,
    partition = partition,
    partition_processes = partition_processes,
    statespace = win,
    state_spaces = state_spaces,
    time_window = c(0, 18),
    treatment_time = 0,
    hawkes_params_control = list(mu = 0.9, A = 0.15, alpha_m = 0.4, c = 0.4,
                                 p = 2.05, D = 1.2, gamma = 0, q = 1.55),
    hawkes_params_treated = list(mu = 0.7, A = 0.18, alpha_m = 0.45, c = 0.4,
                                 p = 2.05, D = 1.2, gamma = 0, q = 1.55),
    update_control_params = TRUE,
    param_update_cadence = 2L,
    update_starting_data = TRUE,
    include_starting_data = TRUE,
    include_starting_first_n = 5L,
    optim_method = "max",
    metric_name = "post_likelihood",
    require_monotone_complete_ll = TRUE,
    iter = 5L, n_props = 4L, change_factor = 0.2,
    proposal_method = "simulation",
    single_flip_from_iter = 3L,
    fixed_params = list(gamma = 0),
    verbose = FALSE,
    model_type = "etas_bivariate",
    outer_maxit_biv = 40L,
    m0 = 2.5, beta_gr = 1.5,
    etas_bivariate_params = init,
    hard_subcritical = TRUE,
    max_branching_radius = 0.98,
    t_trunc = 5
  )
  flips <- as.numeric(res$max_metric_flips)
  expect_equal(length(flips), 5L)
  expect_true(all(flips[3:5] <= 1))
})
