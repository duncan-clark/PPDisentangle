test_that("bivariate ATE uses a supplied reusable cluster", {
  ate_env <- new.env(parent = globalenv())
  sys.source(
    testthat::test_path("..", "..", "inst", "oklahoma", "ate_bivariate.R"),
    envir = ate_env
  )

  calls <- new.env(parent = emptyenv())
  calls$reused <- 0L
  calls$fresh <- 0L
  cluster_token <- list(name = "test-cluster")

  ate_env$run_parallel_on_cluster <- function(cl, X, FUN, label = "job") {
    expect_identical(cl, cluster_token)
    calls$reused <- calls$reused + 1L
    lapply(X, FUN)
  }
  ate_env$run_parallel <- function(X, FUN, cores, label = "job") {
    calls$fresh <- calls$fresh + 1L
    lapply(X, FUN)
  }
  ate_env$sim_etas_bivariate <- function(
      params, windowT, windowS, state_spaces, m0, beta_gr, filtration,
      covariate_lookup, bg_ref_areas, bg_lmax, t_trunc) {
    n <- sample.int(4L, 1L) - 1L
    data.frame(t = seq_len(n))
  }

  full <- spatstat.geom::owin(xrange = c(0, 2), yrange = c(0, 1))
  state_spaces <- list(
    control = spatstat.geom::owin(xrange = c(0, 1), yrange = c(0, 1)),
    treated = spatstat.geom::owin(xrange = c(1, 2), yrange = c(0, 1))
  )
  params <- c(
    mu_0 = 1, mu_1 = 1,
    A_00 = 0.15, alpha_m_00 = 0.3,
    A_11 = 0.15, alpha_m_11 = 0.3,
    A_01 = 0.01, alpha_m_01 = 0.3,
    A_10 = 0.01, alpha_m_10 = 0.3,
    c = 0.2, p = 2.1, D = 1, gamma = 0, q = 1.6
  )
  common <- list(
    biv_params = params,
    windowT = c(0, 2),
    windowS = full,
    state_spaces_obs = state_spaces,
    n_sims = 4L,
    n_cores = 2L,
    m0 = 2.5,
    beta_gr = 2.3,
    crn_base_seed = 100L,
    use_crn = TRUE,
    crn_pair = TRUE,
    quiet = TRUE,
    contrast = "all_or_nothing"
  )

  reused <- do.call(
    ate_env$ate_estim_bivariate,
    c(common, list(parallel_cluster = cluster_token))
  )
  expect_equal(calls$reused, 1L)
  expect_equal(calls$fresh, 0L)

  fresh <- do.call(
    ate_env$ate_estim_bivariate,
    c(common, list(parallel_cluster = NULL))
  )
  expect_equal(calls$reused, 1L)
  expect_equal(calls$fresh, 1L)
  expect_identical(reused, fresh)
})
