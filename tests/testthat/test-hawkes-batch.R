test_that("batched Hawkes filtration loglik matches per-labelling calls", {
  skip_if_not_installed("spatstat.geom")
  set.seed(20260810)
  windowT <- c(0, 40)
  windowS <- spatstat.geom::owin(c(0, 10), c(0, 10))
  treat_time <- 10
  params <- list(mu = 20, alpha = 0.4, beta = 3, K = 0.35)
  sim <- sim_hawkes_fast(params, windowT, windowS)
  dat <- data.frame(x = sim$x, y = sim$y, t = sim$t)
  dat <- dat[order(dat$t), , drop = FALSE]
  n <- nrow(dat)
  skip_if(n < 30)

  set.seed(11)
  base_lab <- ifelse(runif(n) < 0.55, "control", "treated")
  # Pre-treatment events are always control history in SEM.
  base_lab[dat$t < treat_time] <- "control"
  K <- 6L
  labs <- vector("list", K)
  for (k in seq_len(K)) {
    lab <- base_lab
    post <- which(dat$t >= treat_time)
    flip <- sample(post, size = max(1L, round(0.08 * length(post))))
    lab[flip] <- ifelse(lab[flip] == "control", "treated", "control")
    labs[[k]] <- data.frame(dat, inferred_process = lab, stringsAsFactors = FALSE)
  }

  mem <- PPDisentangle:::.hawkes_batch_membership(labs, "control", treat_time, TRUE)
  W <- rep(1.0, length(mem$t))
  areaS <- spatstat.geom::area(windowS)

  batch <- loglik_hawk_filtration_batch(
    params = params,
    t = mem$t, x = mem$x, y = mem$y,
    is_observed = mem$is_observed, member = mem$member,
    W_val = W, areaS = areaS,
    t_start = treat_time, t_end = windowT[2],
    t_trunc = 8
  )

  single <- vapply(seq_len(K), function(k) {
    y <- labs[[k]]
    post <- y[y$t >= treat_time & y$inferred_process == "control", , drop = FALSE]
    filt <- y[y$t < treat_time & y$inferred_process == "control", , drop = FALSE]
    if (nrow(post) < 1L) return(-1e15)
    post <- post[order(post$t), , drop = FALSE]
    filt <- filt[order(filt$t), , drop = FALSE]
    hawkes_loglik_inhom_filtration_cpp(
      post_t = post$t, post_x = post$x, post_y = post$y,
      W_val = rep(1, nrow(post)),
      parent_t = c(filt$t, post$t),
      parent_x = c(filt$x, post$x),
      parent_y = c(filt$y, post$y),
      mu = params$mu, alpha = params$alpha, beta = params$beta, K = params$K,
      areaS = areaS, t_start = treat_time, t_end = windowT[2],
      adjust_factor = 1.0, t_trunc = 8,
      kernel_type = 0L, cc = 1.0, p = 2.0,
      spatial_kernel_type = 0L, spatial_q = 2.0, spatial_d = -1.0
    )
  }, 0)

  expect_equal(batch, single, tolerance = 1e-9)

  # Power-law temporal + spatial branch
  params_pl <- list(mu = 15, alpha = 0.5, c = 0.4, p = 2.1, K = 0.3,
                    kernel = "power_law", spatial_kernel = "power_law",
                    spatial_q = 2.0, spatial_d = 1.2)
  batch_pl <- loglik_hawk_filtration_batch(
    params = params_pl,
    t = mem$t, x = mem$x, y = mem$y,
    is_observed = mem$is_observed, member = mem$member,
    W_val = W, areaS = areaS,
    t_start = treat_time, t_end = windowT[2],
    t_trunc = 8,
    kernel = "power_law", spatial_kernel = "power_law",
    spatial_q = 2.0, spatial_d = 1.2
  )
  single_pl <- vapply(seq_len(K), function(k) {
    y <- labs[[k]]
    post <- y[y$t >= treat_time & y$inferred_process == "control", , drop = FALSE]
    filt <- y[y$t < treat_time & y$inferred_process == "control", , drop = FALSE]
    if (nrow(post) < 1L) return(-1e15)
    post <- post[order(post$t), , drop = FALSE]
    filt <- filt[order(filt$t), , drop = FALSE]
    hawkes_loglik_inhom_filtration_cpp(
      post_t = post$t, post_x = post$x, post_y = post$y,
      W_val = rep(1, nrow(post)),
      parent_t = c(filt$t, post$t),
      parent_x = c(filt$x, post$x),
      parent_y = c(filt$y, post$y),
      mu = params_pl$mu, alpha = params_pl$alpha, beta = 1.0, K = params_pl$K,
      areaS = areaS, t_start = treat_time, t_end = windowT[2],
      adjust_factor = 1.0, t_trunc = 8,
      kernel_type = 1L, cc = params_pl$c, p = params_pl$p,
      spatial_kernel_type = 1L, spatial_q = 2.0, spatial_d = 1.2
    )
  }, 0)
  expect_equal(batch_pl, single_pl, tolerance = 1e-9)
})
