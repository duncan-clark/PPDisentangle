# Sequential MAP proposal + Hamming cap / one-shot teleport.

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

test_that("Hamming cap clips unless teleport is allowed", {
  cur <- data.frame(inferred_process = c("control", "treated", "treated", "treated"),
                    stringsAsFactors = FALSE)
  prop <- data.frame(inferred_process = c("control", "control", "control", "control"),
                     stringsAsFactors = FALSE)
  capped <- .sem_maybe_cap_labelling(
    cur, prop, max_flips_per_step = 1L,
    require_monotone_complete_ll = TRUE,
    allow_improving_teleport = FALSE,
    current_metric = -10, proposed_metric = -5
  )
  expect_false(capped$teleported)
  expect_equal(capped$n_flip, 1L)
  expect_equal(sum(capped$labelling$inferred_process != cur$inferred_process), 1L)

  teleported <- .sem_maybe_cap_labelling(
    cur, prop, max_flips_per_step = 1L,
    require_monotone_complete_ll = TRUE,
    allow_improving_teleport = TRUE,
    current_metric = -10, proposed_metric = -5
  )
  expect_true(teleported$teleported)
  expect_equal(teleported$n_flip, 3L)
  expect_equal(teleported$labelling$inferred_process, prop$inferred_process)

  downhill <- .sem_maybe_cap_labelling(
    cur, prop, max_flips_per_step = 1L,
    require_monotone_complete_ll = TRUE,
    allow_improving_teleport = TRUE,
    current_metric = -5, proposed_metric = -10
  )
  expect_false(downhill$teleported)
  expect_equal(downhill$n_flip, 0L)
})

test_that("one-shot MAP is offered only after the baseline lock", {
  expect_false(.sem_should_include_map_proposal(FALSE, TRUE, FALSE, FALSE))
  expect_true(.sem_should_include_map_proposal(TRUE, FALSE, TRUE, TRUE))
  expect_false(.sem_should_include_map_proposal(TRUE, TRUE, FALSE, TRUE))
  expect_true(.sem_should_include_map_proposal(TRUE, TRUE, FALSE, FALSE))
  expect_false(.sem_should_include_map_proposal(TRUE, TRUE, TRUE, FALSE))
})

test_that("sequential MAP keeps pre-treatment control and assigns by intensity", {
  skip_if_not_installed("spatstat.geom")
  t <- c(-2, -1, 1, 2, 3)
  x <- c(1, 2, 1.1, 8, 8.1)
  y <- c(1, 1, 1.1, 8, 8.1)
  mag <- rep(3, 5)
  W0 <- c(1, 1, 1, 0, 0)
  W1 <- c(0, 0, 0, 1, 1)
  par <- c(
    mu_0 = 0.8, mu_1 = 0.7,
    A_00 = 0.2, alpha_m_00 = 0.4,
    A_11 = 0.2, alpha_m_11 = 0.4,
    A_01 = 0.01, alpha_m_01 = 0.3,
    A_10 = 0.01, alpha_m_10 = 0.3,
    c = 0.3, p = 1.2, D = 1, gamma = 0, q = 1.6
  )
  pid0 <- c(0L, 0L, 1L, 1L, 1L)
  assignable <- c(0L, 0L, 1L, 1L, 1L)
  z <- sequential_map_etas_bivariate(
    par, t, x, y, mag, assignable, pid0, W0, W1, 10, 10, m0 = 2.5, t_trunc = 8
  )
  expect_equal(z[1:2], c(0L, 0L))
  expect_true(all(z %in% c(0L, 1L)))
  expect_equal(length(z), 5L)
})

test_that("clipped SEM stays inside the Hamming cap; one-shot MAP teleports at most once", {
  skip_if_not_installed("spatstat.geom")
  set.seed(21)
  n_pre <- 25L
  n_post <- 40L
  x <- c(runif(n_pre, 0, 5), runif(n_post, 0, 10))
  y <- runif(n_pre + n_post, 0, 10)
  t <- c(sort(runif(n_pre, -8, -0.01)), sort(runif(n_post, 0, 20)))
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
  step_cap <- max(1L, as.integer(ceiling(n_post / n_post)))
  sem_args <- list(
    pp_data = pp,
    partition = partition,
    partition_processes = partition_processes,
    statespace = win,
    state_spaces = state_spaces,
    time_window = c(0, 20),
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
    max_relabel_step_frac = 1 / n_post,
    max_map_relabel_step_frac = 1 / n_post,
    optim_method = "max",
    metric_name = "post_likelihood",
    require_monotone_complete_ll = TRUE,
    include_sequential_map_proposal = TRUE,
    iter = 4L, n_props = 2L, change_factor = 0.01,
    stagnation_trigger_every = 50,
    proposal_method = "simulation",
    fixed_params = list(gamma = 0),
    verbose = FALSE,
    model_type = "etas_bivariate",
    outer_maxit_biv = 40L,
    m0 = 2.5, beta_gr = 1.5,
    etas_bivariate_params = init,
    hard_subcritical = TRUE,
    max_branching_radius = 0.98,
    t_trunc = 8
  )
  blocked <- do.call(em_style_labelling, c(sem_args, list(
    allow_improving_teleport = FALSE,
    map_teleport_once = FALSE
  )))
  once <- do.call(em_style_labelling, c(sem_args, list(
    allow_improving_teleport = FALSE,
    map_teleport_once = TRUE
  )))
  expect_equal(blocked$n_teleports, 0L)
  expect_true(all(blocked$max_metric_flips <= step_cap))
  expect_true(once$n_teleports <= 1L)
  if (once$n_teleports == 0L) {
    expect_true(all(once$max_metric_flips <= step_cap))
  } else {
    expect_equal(sum(once$max_metric_flips > step_cap), 1L)
  }
})

test_that("label_init aliases normalize", {
  expect_equal(.sem_normalize_label_init("current"), "current")
  expect_equal(.sem_normalize_label_init("bernoulli"), "sequential_bernoulli")
  expect_equal(.sem_normalize_label_init("sequential_map"), "sequential_map")
  expect_error(.sem_normalize_label_init("icm"), "Unknown label_init")
})

test_that("sequential Bernoulli keeps pre-treatment control and is reproducible", {
  skip_if_not_installed("spatstat.geom")
  t <- c(-2, -1, 1, 2, 3)
  x <- c(1, 2, 1.1, 8, 8.1)
  y <- c(1, 1, 1.1, 8, 8.1)
  mag <- rep(3, 5)
  W0 <- c(1, 1, 1, 0, 0)
  W1 <- c(0, 0, 0, 1, 1)
  par <- c(
    mu_0 = 0.8, mu_1 = 0.7,
    A_00 = 0.2, alpha_m_00 = 0.4,
    A_11 = 0.2, alpha_m_11 = 0.4,
    A_01 = 0.01, alpha_m_01 = 0.3,
    A_10 = 0.01, alpha_m_10 = 0.3,
    c = 0.3, p = 1.2, D = 1, gamma = 0, q = 1.6
  )
  pid0 <- c(0L, 0L, 1L, 1L, 1L)
  assignable <- c(0L, 0L, 1L, 1L, 1L)
  set.seed(11)
  z1 <- sequential_bernoulli_etas_bivariate(
    par, t, x, y, mag, assignable, pid0, W0, W1, 10, 10, m0 = 2.5, t_trunc = 8
  )
  set.seed(11)
  z2 <- sequential_bernoulli_etas_bivariate(
    par, t, x, y, mag, assignable, pid0, W0, W1, 10, 10, m0 = 2.5, t_trunc = 8
  )
  expect_equal(z1[1:2], c(0L, 0L))
  expect_true(all(z1 %in% c(0L, 1L)))
  expect_identical(z1, z2)
})

test_that("SEM Bernoulli init locks on sampled labels, not location", {
  skip_if_not_installed("spatstat.geom")
  set.seed(21)
  n_pre <- 25L
  n_post <- 40L
  x <- c(runif(n_pre, 0, 5), runif(n_post, 0, 10))
  y <- runif(n_pre + n_post, 0, 10)
  t <- c(sort(runif(n_pre, -8, -0.01)), sort(runif(n_post, 0, 20)))
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
  set.seed(44)
  res <- em_style_labelling(
    pp_data = pp,
    partition = partition,
    partition_processes = partition_processes,
    statespace = win,
    state_spaces = state_spaces,
    time_window = c(0, 20),
    treatment_time = 0,
    hawkes_params_control = list(mu = 0.9, A = 0.15, alpha_m = 0.4, c = 0.4,
                                 p = 2.05, D = 1.2, gamma = 0, q = 1.55),
    hawkes_params_treated = list(mu = 0.7, A = 0.18, alpha_m = 0.45, c = 0.4,
                                 p = 2.05, D = 1.2, gamma = 0, q = 1.55),
    update_control_params = TRUE,
    param_update_cadence = 1L,
    include_starting_data = TRUE,
    include_starting_first_n = 4L,
    max_relabel_step_frac = 0.05,
    optim_method = "max",
    metric_name = "post_likelihood",
    require_monotone_complete_ll = TRUE,
    include_sequential_map_proposal = FALSE,
    map_teleport_once = FALSE,
    label_init = "sequential_bernoulli",
    iter = 3L, n_props = 2L, change_factor = 0.01,
    stagnation_trigger_every = 50,
    proposal_method = "simulation",
    fixed_params = list(gamma = 0),
    verbose = FALSE,
    model_type = "etas_bivariate",
    outer_maxit_biv = 40L,
    m0 = 2.5, beta_gr = 1.5,
    etas_bivariate_params = init,
    hard_subcritical = TRUE,
    max_branching_radius = 0.98,
    t_trunc = 8
  )
  expect_equal(res$label_init, "sequential_bernoulli")
  expect_true(is.finite(res$n_label_init_flips))
  expect_equal(res$max_metric_flips[[1]], 0)
  expect_equal(res$n_teleports, 0L)
})
