.ate_kde_cache_params <- c(
  mu_0 = 2.0, mu_1 = 1.5,
  A_00 = 0.20, alpha_m_00 = 0.40,
  A_11 = 0.18, alpha_m_11 = 0.45,
  A_01 = 0.03, alpha_m_01 = 0.30,
  A_10 = 0.02, alpha_m_10 = 0.35,
  c = 0.2, p = 2.1, D = 1.0, gamma = 0, q = 1.6
)
.ate_kde_cache_full <- spatstat.geom::owin(
  xrange = c(0, 20), yrange = c(0, 10)
)
.ate_kde_cache_spaces <- list(
  control = spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10)),
  treated = spatstat.geom::owin(xrange = c(10, 20), yrange = c(0, 10))
)
.ate_kde_cache_lookup <- list(
  control = function(x, y) 1 + 0.20 * sin(x / 4) * cos(y / 3),
  treated = function(x, y) 1 + 0.15 * cos(x / 5) * sin(y / 2)
)
.ate_kde_cache_filtration <- data.frame(
  x = c(4, 16), y = c(5, 5), t = c(-2, -1),
  mag = c(2.8, 3.0), process_id = c(0L, 1L)
)

test_that("cached KDE bounds preserve the full bivariate simulation stream", {
  ref_areas <- lapply(.ate_kde_cache_spaces, spatstat.geom::area)
  bg_lmax <- list(
    control = PPDisentangle:::.etas_biv_background_lmax(
      .ate_kde_cache_params[["mu_0"]], c(0, 20),
      .ate_kde_cache_spaces$control, .ate_kde_cache_lookup$control,
      ref_areas$control
    ),
    treated = PPDisentangle:::.etas_biv_background_lmax(
      .ate_kde_cache_params[["mu_1"]], c(0, 20),
      .ate_kde_cache_spaces$treated, .ate_kde_cache_lookup$treated,
      ref_areas$treated
    )
  )
  expect_true(all(is.finite(unlist(bg_lmax))))

  common <- list(
    params = .ate_kde_cache_params,
    windowT = c(0, 20),
    windowS = .ate_kde_cache_full,
    state_spaces = .ate_kde_cache_spaces,
    m0 = 2.5,
    beta_gr = 2.3,
    filtration = .ate_kde_cache_filtration,
    covariate_lookup = .ate_kde_cache_lookup,
    bg_ref_areas = ref_areas,
    t_trunc = 5
  )
  set.seed(20260811)
  uncached <- do.call(sim_etas_bivariate, c(common, list(bg_lmax = NULL)))
  set.seed(20260811)
  cached <- do.call(sim_etas_bivariate, c(common, list(bg_lmax = bg_lmax)))

  expect_identical(cached, uncached)
})

test_that("cached KDE bounds preserve complete ATE payloads", {
  ate_env <- new.env(parent = globalenv())
  sys.source(
    testthat::test_path("..", "..", "inst", "oklahoma", "ate_bivariate.R"),
    envir = ate_env
  )
  common <- list(
    biv_params = .ate_kde_cache_params,
    windowT = c(0, 20),
    windowS = .ate_kde_cache_full,
    state_spaces_obs = .ate_kde_cache_spaces,
    label = "KDE cache check",
    n_sims = 12L,
    n_cores = 1L,
    m0 = 2.5,
    beta_gr = 2.3,
    filtration_history = .ate_kde_cache_filtration,
    t_trunc = 5,
    n_tiles = 4L,
    crn_base_seed = 20260811L,
    use_crn = TRUE,
    crn_pair = TRUE,
    quiet = TRUE,
    covariate_lookup = .ate_kde_cache_lookup
  )

  for (contrast in c("all_or_nothing", "observed")) {
    args <- c(common, list(contrast = contrast))
    uncached <- do.call(
      ate_env$ate_estim_bivariate,
      c(args, list(cache_bg_lmax = FALSE))
    )
    cached <- do.call(
      ate_env$ate_estim_bivariate,
      c(args, list(cache_bg_lmax = TRUE))
    )
    expect_identical(cached, uncached, info = contrast)
  }
})
