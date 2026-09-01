test_that("SEM last-labelling vs AOI-boundary loglik uses a common AOI geometry", {
  skip_if_not_installed("spatstat.geom")
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(quiet = TRUE, export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)
  } else if (!exists("loglik_etas_bivariate", mode = "function")) {
    skip_if_not_installed("PPDisentangle")
    library(PPDisentangle)
  }
  geom_candidates <- c(
    file.path("inst", "oklahoma", "oklahoma_geometry.R"),
    file.path("..", "..", "inst", "oklahoma", "oklahoma_geometry.R")
  )
  geom_file <- geom_candidates[file.exists(geom_candidates)][1]
  if (is.na(geom_file) || !nzchar(as.character(geom_file))) {
    skip("oklahoma_geometry.R is not available in this check")
  }
  sys.source(geom_file, envir = environment())
  support_candidates <- c(
    file.path("inst", "oklahoma", "oklahoma_report_support.R"),
    file.path("..", "..", "inst", "oklahoma", "oklahoma_report_support.R")
  )
  support_file <- support_candidates[file.exists(support_candidates)][1]
  if (is.na(support_file) || !nzchar(as.character(support_file))) {
    skip("oklahoma_report_support.R is not available in this check")
  }
  sys.source(support_file, envir = environment())

  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  aoi <- spatstat.geom::owin(c(0, 5), c(0, 10))
  aoi_ss <- oklahoma_report_aoi_state_spaces(win, aoi)
  expect_true(spatstat.geom::inside.owin(2, 5, aoi_ss$treated))
  expect_true(spatstat.geom::inside.owin(8, 5, aoi_ss$control))

  set.seed(20260827)
  holdout <- data.frame(
    x = runif(40, 0.2, 9.8),
    y = runif(40, 0.2, 9.8),
    t = seq(-20, -1, length.out = 40),
    mag = rep(3.0, 40)
  )
  kde_info <- list(bw_sigma = 2.0, bw_sigma_xy = 2.0)
  lambda_im <- oklahoma_report_rebuild_lambda_im(holdout, win, kde_info)
  expect_true(!is.null(lambda_im))

  pp_pre <- data.frame(
    x = c(2.0, 8.0),
    y = c(5.0, 5.0),
    t = c(-4.0, -2.0),
    mag = c(3.1, 3.0),
    location_process = c("treated", "control")
  )
  pp_post <- data.frame(
    x = c(2.5, 7.5, 3.0),
    y = c(4.0, 6.0, 7.0),
    t = c(1.0, 2.0, 3.0),
    mag = c(3.2, 2.9, 3.0),
    location_process = c("treated", "control", "treated")
  )
  aoi_labs <- oklahoma_report_aoi_boundary_labels(rbind(pp_pre, pp_post), aoi_ss$treated)
  expect_equal(
    aoi_labs$inferred_process,
    c("control", "control", "treated", "control", "treated")
  )

  last_lab <- pp_post
  last_lab$inferred_process <- c("treated", "treated", "control")
  catalog <- oklahoma_report_aoi_weighted_catalog(pp_pre, pp_post, aoi_ss, lambda_im)
  last_realiz <- oklahoma_report_apply_last_labels(catalog, last_lab)
  expect_equal(last_realiz$inferred_process[last_realiz$t < 0], c("control", "control"))
  expect_equal(
    last_realiz$inferred_process[last_realiz$t >= 0],
    c("treated", "treated", "control")
  )

  params <- c(
    mu_0 = 0.8, mu_1 = 0.6,
    A_00 = 0.12, alpha_m_00 = 0.35,
    A_11 = 0.10, alpha_m_11 = 0.35,
    A_01 = 0.02, alpha_m_01 = 0.25,
    A_10 = 0.02, alpha_m_10 = 0.25,
    c = 0.2, p = 1.2, D = 1.0, gamma = 0, q = 1.6
  )
  cfg <- list(
    ETAS_M0 = 2.5,
    BETA_GR = 1.5,
    ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS = FALSE,
    ETAS_P_LOWER_BOUND = 0,
    ETAS_Q_LOWER_BOUND = 1,
    ETAS_BRANCHING_MAX = 0.98,
    SEM_T_TRUNC_DAYS = 8,
    windowT_post = c(0, 10)
  )
  mass_ratio <- oklahoma_report_kde_pre_mass_ratio(lambda_im, aoi_ss$control)
  windowT <- c(-4, 10)
  aoi_realiz <- oklahoma_report_aoi_boundary_labels(catalog, aoi_ss$treated)
  ll_last <- oklahoma_report_eval_biv_kde_loglik(
    params, last_realiz, windowT, win, aoi_ss, mass_ratio, cfg
  )
  ll_aoi <- oklahoma_report_eval_biv_kde_loglik(
    params, aoi_realiz, windowT, win, aoi_ss, mass_ratio, cfg
  )
  expect_true(is.finite(ll_last))
  expect_true(is.finite(ll_aoi))
  expect_false(isTRUE(all.equal(ll_last, ll_aoi)))

  same_realiz <- oklahoma_report_apply_last_labels(
    catalog,
    data.frame(
      x = pp_post$x, y = pp_post$y, t = pp_post$t,
      inferred_process = c("treated", "control", "treated")
    )
  )
  ll_same <- oklahoma_report_eval_biv_kde_loglik(
    params, same_realiz, windowT, win, aoi_ss, mass_ratio, cfg
  )
  expect_equal(ll_same, ll_aoi, tolerance = 1e-8)

  tbl <- oklahoma_report_partition_sem_aoi_loglik(
    partition_results = list(
      grid_1.0R = list(
        label = "grid_1.0R",
        F_params = params,
        pp_post_sem = last_lab
      ),
      aoi_region = list(
        label = "aoi_region",
        F_params = params,
        pp_post_sem = data.frame(
          x = pp_post$x, y = pp_post$y, t = pp_post$t,
          inferred_process = c("treated", "control", "treated")
        )
      )
    ),
    primary_id = "grid_1.0R",
    primary_info = NULL,
    pp_pre = pp_pre,
    pp_post = pp_post,
    pp_pre_holdout = holdout,
    win_km = win,
    aoi_owin = aoi,
    kde_info = kde_info,
    cfg = cfg,
    format_label = identity
  )
  expect_equal(nrow(tbl), 2L)
  expect_equal(as.character(tbl$PartitionID[1]), "grid_1.0R")
  expect_true(is.finite(tbl$Last_labelling_LL[1]))
  expect_true(is.finite(tbl$AOI_boundary_LL[1]))
  expect_gt(abs(tbl$Delta_last_minus_AOI[1]), 1e-8)
  expect_equal(tbl$Post_last_ne_AOI[1], 2L)
  expect_equal(tbl$Post_last_ne_AOI[2], 0L)
  expect_equal(tbl$Delta_last_minus_AOI[2], 0, tolerance = 1e-8)
})
