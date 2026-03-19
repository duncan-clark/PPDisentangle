#' Estimate the one-flip treatment effect for cell i
#'
#' Simulates under z_plus (cell i treated) and z_minus (cell i control),
#' returns the difference in expected counts in cell i.
#'
#' @param i Cell index
#' @param partition A spatstat tess object
#' @param treated_partitions Tilenames of treated partitions
#' @param statespace Full observation window
#' @param windowT Numeric vector c(start, end)
#' @param control_pp Control Hawkes parameters
#' @param treated_pp Treated Hawkes parameters
#' @param n_sim Number of simulations
#' @return Scalar estimate of tau_i
#' @export
tau_i <- function(i, partition, treated_partitions, statespace,
                  windowT, control_pp, treated_pp, n_sim) {
  treated_idx <- tilenames(partition) %in% treated_partitions

  z_plus <- treated_idx
  z_plus[i] <- TRUE
  z_minus <- treated_idx
  z_minus[i] <- FALSE

  partition_plus <- rep("control", partition$n)
  partition_plus[z_plus] <- "treated"
  partition_minus <- rep("control", partition$n)
  partition_minus[z_minus] <- "treated"

  plus_control_state_space <- as.owin(partition[!z_plus])
  plus_treated_state_space <- as.owin(partition[z_plus])
  minus_control_state_space <- as.owin(partition[!z_minus])
  minus_treated_state_space <- as.owin(partition[z_minus])

  if (partition_plus[1] == "treated") {
    state_spaces_plus <- list(plus_treated_state_space, plus_control_state_space)
  } else {
    state_spaces_plus <- list(plus_control_state_space, plus_treated_state_space)
  }

  if (partition_minus[1] == "treated") {
    state_spaces_minus <- list(minus_treated_state_space, minus_control_state_space)
  } else {
    state_spaces_minus <- list(minus_control_state_space, minus_treated_state_space)
  }

  plus_total <- 0
  minus_total <- 0
  for (j in seq_len(n_sim)) {
    pp <- generate_inhomogeneous_hawkes(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_plus,
      hawkes_params = list(control = control_pp, treated = treated_pp),
      state_spaces = state_spaces_plus, space_triggering = FALSE
    )
    plus_total <- plus_total + sum(as.numeric(tileindex(pp$x, pp$y, partition)) == i)

    pp <- generate_inhomogeneous_hawkes(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_minus,
      hawkes_params = list(control = control_pp, treated = treated_pp),
      state_spaces = state_spaces_minus, space_triggering = FALSE
    )
    minus_total <- minus_total + sum(as.numeric(tileindex(pp$x, pp$y, partition)) == i)
  }
  return(plus_total / n_sim - minus_total / n_sim)
}

#' Non-parametric ATE estimation from labeled data
#'
#' @param statespace Full observation window
#' @param partition A spatstat tess object
#' @param observed_data Data frame with inferred_process column
#' @param treated_partitions Tilenames of treated partitions
#' @param pp_rhs_form Formula for intensity (unused, kept for API compat)
#' @param windowT Numeric vector c(start, end)
#' @return List with ATE_non_param
#' @export
ATE_estim <- function(statespace, partition, observed_data, treated_partitions,
                      pp_rhs_form = "~1", windowT = c(0, 1)) {
  treated_idx <- tilenames(partition) %in% treated_partitions
  control_points_per_tile <- sum(observed_data$inferred_process == "control") / sum(!treated_idx)
  treated_points_per_tile <- sum(observed_data$inferred_process == "treated") / sum(treated_idx)
  ATE_naive <- treated_points_per_tile - control_points_per_tile
  return(list(ATE_non_param = ATE_naive, ATE_non_param_z = ATE_naive))
}

#' Parametric ATE estimation using fitted Hawkes processes
#'
#' Fits Hawkes models to control/treated subsets, then estimates ATE via
#' simulation (all-or-nothing and one-flip).
#'
#' @param statespace Full observation window
#' @param partition A spatstat tess object
#' @param observed_data Data frame with labeled point process data
#' @param treated_partitions Tilenames of treated partitions
#' @param hawkes_params Optional pre-fitted params (list with control, treated)
#' @param n_tau_sims Number of sims per tau_i estimate
#' @param n_tau_i Number of tau_i estimates
#' @param n_sims Number of all-or-nothing simulations
#' @param windowT Numeric vector c(start, end)
#' @param windowS Observation window for fitting
#' @param maxit Max optim iterations
#' @param poisson_flags List with control and treated Poisson flags
#' @param filtration_data Optional pre-treatment history rows used as Hawkes
#'   filtration for control-parameter fitting from labeled post data.
#' @return List with all-or-nothing sims, tau estimates, and decomposed ATEs
#' @export
fit_hawkes_with_filtration <- function(params_init,
                                       realiz,
                                       filtration = NULL,
                                       windowT = c(0, 1),
                                       windowS = c(0, 1, 0, 1),
                                       maxit = 1000,
                                       poisson_flag = FALSE,
                                       zero_background_region = NULL) {
  if (!inherits(windowS, "owin")) windowS <- as.owin(windowS)
  if (inherits(params_init, "list")) params_init <- unlist(params_init)
  params_init <- as.numeric(params_init[c("mu", "alpha", "beta", "K")])
  names(params_init) <- c("mu", "alpha", "beta", "K")
  params_init[!is.finite(params_init)] <- c(1, 0.01, 1, 0.01)[!is.finite(params_init)]

  realiz <- as.data.frame(realiz)
  if (nrow(realiz) > 0) {
    realiz <- realiz[realiz$t >= windowT[1] & realiz$t <= windowT[2], , drop = FALSE]
    realiz <- realiz[order(realiz$t), , drop = FALSE]
  }
  if (nrow(realiz) < 1L) {
    return(list(par = list(mu = 0, alpha = 0.01, beta = 1, K = 0), converged = TRUE))
  }

  if (is.null(filtration)) {
    filtration <- realiz[0, c("x", "y", "t"), drop = FALSE]
  } else {
    filtration <- as.data.frame(filtration)
    filtration <- filtration[filtration$t < windowT[1], , drop = FALSE]
    filtration <- filtration[order(filtration$t), , drop = FALSE]
  }

  total_area <- spatstat.geom::area(windowS)
  active_area <- total_area
  in_zero <- rep(FALSE, nrow(realiz))
  if (!is.null(zero_background_region)) {
    if (!inherits(zero_background_region, "owin")) {
      zero_background_region <- as.owin(zero_background_region)
    }
    zero_area <- spatstat.geom::area(zero_background_region)
    active_area <- max(1e-12, total_area - zero_area)
    in_zero <- inside.owin(realiz$x, realiz$y, w = zero_background_region)
  }

  if (isTRUE(poisson_flag)) {
    dt <- windowT[2] - windowT[1]
    mu_hat <- if (dt > 0) nrow(realiz) / dt else 0
    return(list(par = list(mu = mu_hat, alpha = 1, beta = 1, K = 0), converged = TRUE))
  }

  W_post <- if ("W" %in% names(realiz)) as.numeric(realiz$W) else rep(1, nrow(realiz))
  W_post[!is.finite(W_post)] <- 0
  W_post[in_zero] <- 0

  # Use pre-treatment filtration as parent-history only. We evaluate a
  # conditional post-treatment likelihood to avoid fitting pre+post regimes
  # with a single stationary parameter vector.
  post_t_fit <- as.numeric(realiz$t)
  post_x_fit <- as.numeric(realiz$x)
  post_y_fit <- as.numeric(realiz$y)
  W_fit <- as.numeric(W_post)
  parent_x <- c(filtration$x, realiz$x)
  parent_y <- c(filtration$y, realiz$y)
  parent_t <- c(filtration$t, realiz$t)
  t_start_fit <- windowT[1]
  t_end_fit <- windowT[2]
  adjust_factor_fit <- 1

  dt <- windowT[2] - windowT[1]
  min_trigger_sd <- 0.001 * sqrt(total_area)
  alpha_max <- 1 / min_trigger_sd^2
  beta_min <- if (dt > 0) 0.1 / dt else 0

  obj_fn <- function(par4) {
    mu <- par4[1]; alpha <- par4[2]; beta <- par4[3]; K <- par4[4]
    if (!is.finite(mu) || !is.finite(alpha) || !is.finite(beta) || !is.finite(K)) return(-1e15)
    if (mu < 0 || alpha < 0 || beta <= 0 || K < 0 || K >= 1) return(-1e15)
    if (alpha > alpha_max || beta < beta_min) return(-1e15)
    loglik <- hawkes_loglik_inhom_filtration_cpp(
      post_t = post_t_fit,
      post_x = post_x_fit,
      post_y = post_y_fit,
      W_val = W_fit,
      parent_t = parent_t,
      parent_x = parent_x,
      parent_y = parent_y,
      mu = mu,
      alpha = alpha,
      beta = beta,
      K = K,
      areaS = active_area,
      t_start = t_start_fit,
      t_end = t_end_fit,
      adjust_factor = adjust_factor_fit,
      t_trunc = -1
    )
    if (!is.finite(loglik)) return(-1e15)
    loglik
  }

  fit <- tryCatch(
    optim(
      par = params_init,
      fn = obj_fn,
      method = "Nelder-Mead",
      control = list(fnscale = -1, trace = 0, maxit = maxit)
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || is.null(fit$par) || any(!is.finite(fit$par))) {
    return(list(par = as.list(params_init), converged = FALSE))
  }
  pv <- as.numeric(fit$par)
  names(pv) <- c("mu", "alpha", "beta", "K")
  list(par = as.list(pv), converged = is.null(fit$convergence) || fit$convergence == 0)
}

ATE_estim_hawkes <- function(statespace, partition, observed_data, treated_partitions,
                             hawkes_params = NULL, n_tau_sims = 10, n_tau_i = 10,
                             n_sims = 100, windowT = c(0, 1), windowS = c(0, 1, 0, 1),
                             maxit = 1000, poisson_flags = list(control = FALSE, treated = FALSE),
                             filtration_data = NULL, explosive_K_threshold = 0.98,
                             control_filtration_aware = TRUE) {
  treated_idx <- tilenames(partition) %in% treated_partitions
  control_state_space <- as.owin(partition[!treated_idx])
  treated_state_space <- as.owin(partition[treated_idx])

  control_mean <- sum(observed_data$location_process == "control") / sum(!treated_idx)
  treated_mean <- sum(observed_data$location_process == "treated") / sum(treated_idx)
  control_points_in_control_cells <- sum(observed_data$inferred_process == "control" & observed_data$location_process == "control") / sum(!treated_idx)
  treated_points_in_treated_cells <- sum(observed_data$inferred_process == "treated" & observed_data$location_process == "treated") / sum(treated_idx)

  ATE_total <- treated_mean - control_points_in_control_cells
  ATE_treatment <- treated_mean - control_mean
  ATE_spillover <- control_mean - control_points_in_control_cells
  ATE_naive <- ATE_treatment

  if (is.null(hawkes_params)) {
    filt <- if (is.null(filtration_data)) {
      observed_data[0, , drop = FALSE]
    } else {
      as.data.frame(filtration_data)
    }
    if (!("inferred_process" %in% names(filt)) && ("location_process" %in% names(filt))) {
      filt$inferred_process <- filt$location_process
    }
    if (!("inferred_process" %in% names(filt)) && nrow(filt) > 0) {
      filt$inferred_process <- "control"
    }

    ctrl_realiz <- observed_data[observed_data$inferred_process == "control", , drop = FALSE]
    treat_realiz <- observed_data[observed_data$inferred_process == "treated", , drop = FALSE]
    dt_fit <- windowT[2] - windowT[1]
    if (!is.finite(dt_fit) || dt_fit <= 0) dt_fit <- 1

    ctrl_init <- list(
      mu = max(1e-8, nrow(ctrl_realiz) / dt_fit),
      alpha = 0,
      beta = dt_fit / 100,
      K = 0.01
    )
    treat_init <- list(
      mu = max(1e-8, nrow(treat_realiz) / dt_fit),
      alpha = 0,
      beta = dt_fit / 100,
      K = 0.01
    )

    if (isTRUE(control_filtration_aware)) {
      # Control fit can optionally use full pre-treatment filtration.
      control_pp <- fit_hawkes_with_filtration(
        params_init = ctrl_init,
        realiz = ctrl_realiz,
        filtration = filt,
        windowT = windowT,
        windowS = windowS,
        maxit = maxit,
        poisson_flag = poisson_flags$control,
        zero_background_region = treated_state_space
      )$par

      # Guard against near-critical/near-zero-baseline filtration fits:
      # if the implied stationary rate is far below observed post-treatment
      # control counts, fall back to the legacy post-only control fit.
      empirical_rate <- if (dt_fit > 0) nrow(ctrl_realiz) / dt_fit else Inf
      fitted_rate <- as.numeric(control_pp$mu) / max(1e-6, 1 - as.numeric(control_pp$K))
      degenerate_fit <- (!is.finite(fitted_rate) || !is.finite(empirical_rate) ||
                         (as.numeric(control_pp$K) >= 0.98) ||
                         (fitted_rate < 0.2 * empirical_rate))
      if (isTRUE(degenerate_fit)) {
        legacy <- fit_hawkes(
          unlist(ctrl_init),
          realiz = ctrl_realiz,
          zero_background_region = treated_state_space,
          windowT = windowT,
          windowS = windowS,
          trace = 0,
          maxit = maxit,
          density_approx = FALSE,
          numeric_integral = FALSE,
          poisson_flag = poisson_flags$control,
          t_trunc = -1
        )$par
        control_pp <- as.list(legacy)
        names(control_pp) <- c("mu", "alpha", "beta", "K")
      }
    } else {
      # Legacy behavior: control fit uses only post-treatment observed control points.
      control_pp <- fit_hawkes(
        unlist(ctrl_init),
        realiz = ctrl_realiz,
        zero_background_region = treated_state_space,
        windowT = windowT,
        windowS = windowS,
        trace = 0,
        maxit = maxit,
        density_approx = FALSE,
        numeric_integral = FALSE,
        poisson_flag = poisson_flags$control,
        t_trunc = -1
      )$par
      control_pp <- as.list(control_pp)
      names(control_pp) <- c("mu", "alpha", "beta", "K")
    }

    # Keep treated fitting behavior aligned with the historical (pre-filtration) path.
    treated_pp <- fit_hawkes(
      unlist(treat_init),
      realiz = treat_realiz,
      zero_background_region = control_state_space,
      windowT = windowT,
      windowS = windowS,
      trace = 0,
      maxit = maxit,
      density_approx = FALSE,
      numeric_integral = FALSE,
      poisson_flag = poisson_flags$treated,
      t_trunc = -1
    )$par
    treated_pp <- as.list(treated_pp)
    names(treated_pp) <- c("mu", "alpha", "beta", "K")
  } else {
    control_pp <- hawkes_params$control
    treated_pp <- hawkes_params$treated
  }

  is_explosive <- function(pp, k_thr) {
    if (is.null(pp) || is.null(pp$K) || !is.finite(pp$K)) return(TRUE)
    as.numeric(pp$K) > as.numeric(k_thr)
  }
  if (is_explosive(control_pp, explosive_K_threshold) ||
      is_explosive(treated_pp, explosive_K_threshold)) {
    all_nothing_theory <- data.frame(c_mean = NA_real_, t_mean = NA_real_, ATE = NA_real_)
    return(list(
      all_nothing_sim = NULL,
      all_nothing_theory = all_nothing_theory,
      tau_1_estim = NA_real_,
      ATE_total = ATE_total, ATE_treatment = ATE_treatment,
      ATE_spillover = ATE_spillover, ATE_naive = ATE_naive,
      treated_pp = treated_pp, control_pp = control_pp,
      skipped_explosive = TRUE,
      explosive_K_threshold = explosive_K_threshold
    ))
  }

  partition_process <- rep("control", partition$n)
  partition_process[treated_idx] <- "treated"

  tau_i_estim <- vapply(seq_len(n_tau_i), function(j) {
    tau_i(
      sample(length(partition_process), 1),
      partition = partition, treated_partitions = treated_partitions,
      statespace = statespace, windowT = windowT,
      control_pp = control_pp, treated_pp = treated_pp, n_sim = n_tau_sims
    )
  }, numeric(1))
  tau_1_estim <- mean(tau_i_estim)

  c_counts <- numeric(n_sims)
  t_counts <- numeric(n_sims)
  for (s in seq_len(n_sims)) {
    c_counts[s] <- sim_hawkes(control_pp, windowT, windowS = statespace)$n[1]
    t_counts[s] <- sim_hawkes(treated_pp, windowT, windowS = statespace)$n[1]
  }
  all_nothing_sim <- data.frame(
    c_mean = c_counts / partition$n,
    t_mean = t_counts / partition$n,
    ATE = (t_counts - c_counts) / partition$n
  )

  all_nothing_theory <- data.frame(
    c_mean = control_pp$mu * (windowT[2] - windowT[1]) * (1 / (1 - control_pp$K)) / partition$n,
    t_mean = treated_pp$mu * (windowT[2] - windowT[1]) * (1 / (1 - treated_pp$K)) / partition$n,
    ATE = (treated_pp$mu * (windowT[2] - windowT[1]) * (1 / (1 - treated_pp$K)) -
             control_pp$mu * (windowT[2] - windowT[1]) * (1 / (1 - control_pp$K))) / partition$n
  )

  return(list(
    all_nothing_sim = all_nothing_sim,
    all_nothing_theory = all_nothing_theory,
    tau_1_estim = tau_1_estim,
    ATE_total = ATE_total, ATE_treatment = ATE_treatment,
    ATE_spillover = ATE_spillover, ATE_naive = ATE_naive,
    treated_pp = treated_pp, control_pp = control_pp
  ))
}


#' Estimate one-flip treatment effect for cell i (ETAS model)
#'
#' Analogue of \code{tau_i} for the ETAS model.
#'
#' @inheritParams tau_i
#' @param control_pp  Control ETAS parameter list.
#' @param treated_pp  Treated ETAS parameter list.
#' @param m0  Reference magnitude.
#' @param beta_gr  Gutenberg-Richter beta (or NULL for resampling).
#' @param mag_pool  Magnitude pool for resampling.
#' @return Scalar estimate of tau_i.
#' @export
tau_i_etas <- function(i, partition, treated_partitions, statespace,
                       windowT, control_pp, treated_pp, n_sim,
                       m0, beta_gr = NULL, mag_pool = NULL) {
  treated_idx <- tilenames(partition) %in% treated_partitions

  z_plus <- treated_idx;  z_plus[i]  <- TRUE
  z_minus <- treated_idx; z_minus[i] <- FALSE

  partition_plus  <- ifelse(z_plus, "treated", "control")
  partition_minus <- ifelse(z_minus, "treated", "control")

  make_ss <- function(z) {
    pr <- unique(ifelse(z, "treated", "control"))
    ss <- list()
    for (p in pr) ss[[p]] <- as.owin(partition[ifelse(z, "treated", "control") == p])
    ss
  }

  plus_total <- 0; minus_total <- 0
  for (j in seq_len(n_sim)) {
    pp <- generate_inhomogeneous_etas(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_plus,
      etas_params = list(control = control_pp, treated = treated_pp),
      m0 = m0, beta_gr = beta_gr, mag_pool = mag_pool
    )
    plus_total <- plus_total + sum(as.numeric(tileindex(pp$x, pp$y, partition)) == i)

    pp <- generate_inhomogeneous_etas(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_minus,
      etas_params = list(control = control_pp, treated = treated_pp),
      m0 = m0, beta_gr = beta_gr, mag_pool = mag_pool
    )
    minus_total <- minus_total + sum(as.numeric(tileindex(pp$x, pp$y, partition)) == i)
  }
  return(plus_total / n_sim - minus_total / n_sim)
}


#' Parametric ATE estimation using fitted ETAS processes
#'
#' ETAS analogue of \code{ATE_estim_hawkes}.  Fits ETAS models to
#' control/treated subsets, then estimates ATE via simulation.
#'
#' @inheritParams ATE_estim_hawkes
#' @param etas_params  Optional pre-fitted ETAS params
#'   (list with control, treated).
#' @param m0  Reference magnitude.
#' @param beta_gr  Gutenberg-Richter beta (or NULL for resampling).
#' @param mag_pool  Numeric vector of magnitudes for resampling.
#' @param fixed_params  Named list of ETAS params to hold fixed.
#' @return List with all-or-nothing sims, tau estimates, decomposed ATEs,
#'   and fitted parameter lists.
#' @export
ATE_estim_etas <- function(statespace, partition, observed_data,
                           treated_partitions,
                           etas_params = NULL,
                           n_tau_sims = 10, n_tau_i = 10,
                           n_sims = 100,
                           windowT = c(0, 1),
                           windowS = NULL,
                           maxit = 3000,
                           m0 = NULL,
                           beta_gr = NULL,
                           mag_pool = NULL,
                           fixed_params = NULL) {
  if (is.null(windowS)) windowS <- statespace
  treated_idx <- tilenames(partition) %in% treated_partitions
  control_state_space <- as.owin(partition[!treated_idx])
  treated_state_space <- as.owin(partition[treated_idx])

  if (is.null(m0)) m0 <- min(observed_data$mag, na.rm = TRUE)

  control_mean <- sum(observed_data$location_process == "control") / sum(!treated_idx)
  treated_mean <- sum(observed_data$location_process == "treated") / sum(treated_idx)
  control_in_ctrl <- sum(observed_data$inferred_process == "control" &
                         observed_data$location_process == "control") / sum(!treated_idx)
  treated_in_treat <- sum(observed_data$inferred_process == "treated" &
                          observed_data$location_process == "treated") / sum(treated_idx)

  ATE_total     <- treated_mean - control_in_ctrl
  ATE_treatment <- treated_mean - control_mean
  ATE_spillover <- control_mean - control_in_ctrl
  ATE_naive     <- ATE_treatment

  if (is.null(etas_params)) {
    dt <- windowT[2] - windowT[1]
    init_ctrl <- list(mu = sum(observed_data$inferred_process == "control") / dt,
                      A = 0.1, alpha_m = 0.5, c = 0.05, p = 1.2,
                      D = 1.0, gamma = 0.3, q = 1.5)
    init_treat <- list(mu = sum(observed_data$inferred_process == "treated") / dt,
                       A = 0.1, alpha_m = 0.5, c = 0.05, p = 1.2,
                       D = 1.0, gamma = 0.3, q = 1.5)

    ctrl_fit <- fit_etas(
      params_init = init_ctrl,
      realiz = observed_data[observed_data$inferred_process == "control", ],
      windowT = windowT, windowS = windowS, m0 = m0,
      maxit = maxit, fixed_params = fixed_params,
      zero_background_region = treated_state_space
    )
    control_pp <- as.list(ctrl_fit$par)
    names(control_pp) <- .etas_par_names

    treat_fit <- fit_etas(
      params_init = init_treat,
      realiz = observed_data[observed_data$inferred_process == "treated", ],
      windowT = windowT, windowS = windowS, m0 = m0,
      maxit = maxit, fixed_params = fixed_params,
      zero_background_region = control_state_space
    )
    treated_pp <- as.list(treat_fit$par)
    names(treated_pp) <- .etas_par_names
  } else {
    control_pp <- etas_params$control
    treated_pp <- etas_params$treated
  }

  partition_process <- rep("control", partition$n)
  partition_process[treated_idx] <- "treated"

  tau_i_estim <- vapply(seq_len(n_tau_i), function(j) {
    tau_i_etas(
      sample(length(partition_process), 1),
      partition = partition, treated_partitions = treated_partitions,
      statespace = statespace, windowT = windowT,
      control_pp = control_pp, treated_pp = treated_pp, n_sim = n_tau_sims,
      m0 = m0, beta_gr = beta_gr, mag_pool = mag_pool
    )
  }, numeric(1))
  tau_1_estim <- mean(tau_i_estim)

  c_counts <- numeric(n_sims)
  t_counts <- numeric(n_sims)
  for (s in seq_len(n_sims)) {
    c_sim <- sim_etas(control_pp, windowT, windowS = statespace,
                      m0 = m0, beta_gr = beta_gr, mag_pool = mag_pool)
    t_sim <- sim_etas(treated_pp, windowT, windowS = statespace,
                      m0 = m0, beta_gr = beta_gr, mag_pool = mag_pool)
    c_counts[s] <- length(c_sim$t)
    t_counts[s] <- length(t_sim$t)
  }
  all_nothing_sim <- data.frame(
    c_mean = c_counts / partition$n,
    t_mean = t_counts / partition$n,
    ATE    = (t_counts - c_counts) / partition$n
  )

  analytic <- ATE_analytic_etas(control_pp, treated_pp,
                                windowT = windowT, n_tiles = partition$n,
                                beta_gr = beta_gr, m0 = m0)

  return(list(
    all_nothing_sim = all_nothing_sim,
    tau_1_estim = tau_1_estim,
    ATE_total = ATE_total, ATE_treatment = ATE_treatment,
    ATE_spillover = ATE_spillover, ATE_naive = ATE_naive,
    treated_pp = treated_pp, control_pp = control_pp,
    analytic = analytic
  ))
}


#' Analytic expected count for an ETAS process
#'
#' For a univariate ETAS with Gutenberg-Richter magnitudes (exponential
#' rate beta_gr), the expected total count on the full spatial domain
#' over a time window of length T is  mu * T / (1 - eta),
#' where eta = A * beta_gr / (beta_gr - alpha_m) is the branching ratio.
#'
#' When eta >= 1 (supercritical) the infinite-time expectation diverges.
#' In this case we report both the background-only count (mu * T) and
#' mark the total as NA (must use simulation).
#'
#' @param params  Named list with at least mu, A, alpha_m.
#' @param windowT Numeric c(start, end).
#' @param beta_gr Gutenberg-Richter beta.
#' @param m0      Reference magnitude (unused in the formula but kept for API).
#' @return List with E_N (expected count, NA if supercritical), E_N_bg
#'   (background-only count), and eta (branching ratio).
#' @export
etas_expected_count <- function(params, windowT, beta_gr, m0 = NULL) {
  if (!is.list(params)) params <- as.list(params)
  mu      <- params$mu
  A       <- params$A
  alpha_m <- params$alpha_m
  dt <- windowT[2] - windowT[1]
  E_N_bg <- mu * dt

  if (beta_gr <= alpha_m) {
    eta <- NA_real_
    E_N <- NA_real_
  } else {
    eta <- A * beta_gr / (beta_gr - alpha_m)
    if (eta < 1) {
      E_N <- mu * dt / (1 - eta)
    } else {
      E_N <- NA_real_
    }
  }
  list(E_N = E_N, E_N_bg = E_N_bg, eta = eta, mu = mu, A = A, alpha_m = alpha_m)
}


#' Analytic ATE for independent ETAS models (all-or-nothing)
#'
#' Computes the per-tile difference in expected counts between the
#' treated and control ETAS processes analytically. When the branching
#' ratio eta < 1 (subcritical), uses the full formula mu*T/(1-eta).
#' Always reports the background-rate ATE (difference in mu*T) which
#' is valid regardless of criticality.
#'
#' @param control_pp  Control ETAS parameters (list).
#' @param treated_pp  Treated ETAS parameters (list).
#' @param windowT     Numeric c(start, end).
#' @param n_tiles     Number of spatial tiles.
#' @param beta_gr     Gutenberg-Richter beta.
#' @param m0          Reference magnitude.
#' @return List with ATE (NA if supercritical), ATE_bg (background-rate
#'   ATE, always valid), E_N_ctrl, E_N_treat, eta_ctrl, eta_treat.
#' @export
ATE_analytic_etas <- function(control_pp, treated_pp, windowT,
                              n_tiles, beta_gr, m0 = NULL) {
  ctrl <- etas_expected_count(control_pp, windowT, beta_gr, m0)
  treat <- etas_expected_count(treated_pp, windowT, beta_gr, m0)
  ATE <- (treat$E_N - ctrl$E_N) / n_tiles
  ATE_bg <- (treat$E_N_bg - ctrl$E_N_bg) / n_tiles
  list(ATE = ATE, ATE_bg = ATE_bg,
       E_N_ctrl = ctrl$E_N, E_N_treat = treat$E_N,
       E_N_bg_ctrl = ctrl$E_N_bg, E_N_bg_treat = treat$E_N_bg,
       eta_ctrl = ctrl$eta, eta_treat = treat$eta)
}


#' Louis-method standard errors for SEM ATE estimates
#'
#' Implements the Louis (1982) observed information identity for the SEM.
#' Given K importance-weighted labellings with log-likelihoods and
#' parameter estimates, computes the variance of the ATE estimator
#' accounting for both sampling variability and label uncertainty.
#'
#' The complete-data score variance (missing information) is estimated
#' from the weighted variance of the per-labelling ATE estimates.
#' The observed information combines the curvature at the MLE with
#' this missing information correction.
#'
#' @param labellings   List of labelling data frames from the SEM.
#' @param weights      Numeric vector of normalized importance weights.
#' @param ate_fn       Function(labelling) -> scalar ATE for one labelling.
#' @param loglik_fn    Optional function(labelling) -> scalar loglik.
#' @param n_boot       Number of weighted bootstrap resamples for SE.
#' @return List with ate_mean, ate_se, ate_values, ci_lower, ci_upper.
#' @export
louis_se_ate <- function(labellings, weights, ate_fn, loglik_fn = NULL,
                         n_boot = 1000) {
  K <- length(labellings)
  if (length(weights) != K) stop("weights must have same length as labellings")

  ate_vals <- vapply(labellings, ate_fn, numeric(1))
  ate_mean <- sum(weights * ate_vals)

  # Missing information: weighted variance of ATE across labellings
  ate_var_missing <- sum(weights * (ate_vals - ate_mean)^2)

  # Weighted bootstrap for total SE (accounts for effective sample size)
  boot_means <- replicate(n_boot, {
    idx <- sample(K, K, replace = TRUE, prob = weights)
    mean(ate_vals[idx])
  })
  ate_se <- sd(boot_means)

  ess <- 1 / sum(weights^2)
  ate_se <- max(ate_se, sqrt(ate_var_missing))

  list(ate_mean = ate_mean, ate_se = ate_se,
       ate_values = ate_vals, weights = weights,
       ci_lower = ate_mean - 1.96 * ate_se,
       ci_upper = ate_mean + 1.96 * ate_se,
       ate_var_missing = ate_var_missing, ess = ess)
}
