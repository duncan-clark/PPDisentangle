#' Adaptive Stochastic EM algorithm for point process disentangling
#'
#' Outer loop of the SEM: runs adaptive labeling proposals, computes
#' importance weights, and iteratively re-estimates Hawkes parameters.
#'
#' @param pp_data Data frame with columns x, y, t, process, location_process
#' @param partition A spatstat tess object
#' @param partition_processes Character vector of process names per tile
#' @param statespace Full observation window (owin)
#' @param time_window Numeric vector c(start, end) (used for internal naming)
#' @param treatment_time Scalar treatment time
#' @param hawkes_params_control Initial control Hawkes parameters
#' @param hawkes_params_treated Initial treated Hawkes parameters
#' @param N_labellings Number of labeling proposals per outer iteration
#' @param N_iter Number of outer EM iterations
#' @param reset Logical; if TRUE, the outer loop counter resets to 0 whenever the adaptive step is triggered by weight concentration.
#' @param verbose Print progress
#' @param adaptive_control List of adaptive control parameters
#' @param ... Additional arguments passed through (e.g. background_rate_var, covariate_lookup)
#' @return List with hawkes_params_control, hawkes_params_treated, t_params,
#'   labellings, adaptive results, and timing
#' @export
adaptive_SEM <- function(pp_data,
                         partition,
                         partition_processes,
                         statespace,
                         time_window,
                         treatment_time,
                         hawkes_params_control,
                         hawkes_params_treated,
                         N_labellings,
                         N_iter,
                         reset = FALSE,
                         verbose = TRUE,
                         adaptive_control = list(
                           param_update_cadence = 20,
                           proposal_update_cadence = 1,
                           state_spaces = NULL,
                           iter = 100,
                           n_props = 10,
                           change_factor = 0.1,
                           verbose = FALSE,
                           update_starting_data = TRUE,
                           include_starting_data = FALSE
                         ),
                         model_type = "hawkes",
                         ...) {
  t_global <- proc.time()[3]
  N_labellings <- suppressWarnings(as.integer(N_labellings))
  if (!is.finite(N_labellings) || is.na(N_labellings) || N_labellings < 0L) {
    stop("N_labellings must be a non-negative integer.")
  }
  dots <- list(...)
  background_rate_var <- if ("background_rate_var" %in% names(dots)) dots$background_rate_var else NULL
  t_trunc <- if ("t_trunc" %in% names(dots)) dots$t_trunc else NULL
  use_pre_history_for_biv <- isTRUE(dots$use_pre_history_for_biv)
  hawkes_use_filtration_history <- if ("hawkes_use_filtration_history" %in% names(dots)) {
    isTRUE(dots$hawkes_use_filtration_history)
  } else {
    TRUE
  }
  etas_use_filtration_history <- if ("etas_use_filtration_history" %in% names(dots)) {
    isTRUE(dots$etas_use_filtration_history)
  } else {
    TRUE
  }
  treated_background_zero_before <- if ("treated_background_zero_before" %in% names(dots)) {
    as.numeric(dots$treated_background_zero_before)
  } else {
    NULL
  }
  control_background_everywhere_before <- if ("control_background_everywhere_before" %in% names(dots)) {
    .etas_bg_cutoff(dots$control_background_everywhere_before)
  } else {
    NULL
  }
  control_background_pre_mass_ratio <- if ("control_background_pre_mass_ratio" %in% names(dots)) {
    suppressWarnings(as.numeric(dots$control_background_pre_mass_ratio))
  } else {
    NULL
  }
  dots_no_trunc <- dots
  dots_no_trunc$t_trunc <- NULL

  is_etas <- identical(model_type, "etas")
  is_biv_etas <- identical(model_type, "etas_bivariate")
  etas_beta_eff <- .etas_resolve_beta_gr(
    dots$beta_gr,
    realiz = pp_data,
    m0 = dots$m0
  )
  etas_gap_min <- suppressWarnings(as.numeric(dots$alpha_beta_gap_min))
  if (length(etas_gap_min) != 1L || !is.finite(etas_gap_min) ||
      is.na(etas_gap_min) || etas_gap_min < 0) {
    etas_gap_min <- 1e-4
  }
  etas_eta_max <- suppressWarnings(as.numeric(dots$max_branching_ratio))
  if (length(etas_eta_max) != 1L || !is.finite(etas_eta_max) ||
      is.na(etas_eta_max) || etas_eta_max <= 0) {
    etas_eta_max <- 0.98
  }
  etas_rho_max <- suppressWarnings(as.numeric(dots$max_branching_radius))
  if (length(etas_rho_max) != 1L || !is.finite(etas_rho_max) ||
      is.na(etas_rho_max) || etas_rho_max <= 0) {
    etas_rho_max <- 0.98
  }
  etas_alpha_lo <- suppressWarnings(as.numeric(dots$alpha_m_lower_bound))
  if (length(etas_alpha_lo) != 1L || is.na(etas_alpha_lo)) {
    etas_alpha_lo <- 0
  }
  etas_init_margin <- suppressWarnings(as.numeric(dots$init_branching_margin))
  if (length(etas_init_margin) != 1L || !is.finite(etas_init_margin) ||
      is.na(etas_init_margin) || etas_init_margin <= 0 || etas_init_margin > 1) {
    etas_init_margin <- 0.9
  }
  etas_soft_barrier <- !isFALSE(dots$soft_branching_barrier)
  etas_log_transform <- !isFALSE(dots$log_transform)
  etas_near_cap_frac <- suppressWarnings(as.numeric(dots$near_cap_frac))
  if (length(etas_near_cap_frac) != 1L || !is.finite(etas_near_cap_frac) ||
      etas_near_cap_frac <= 0 || etas_near_cap_frac > 1) {
    etas_near_cap_frac <- 0.99
  }
  univ_barrier_ctrl <- .etas_soft_barrier_controls(
    etas_eta_max,
    start = dots$stability_barrier_start,
    weight = if (!is.null(dots$stability_barrier_weight) &&
                 is.finite(suppressWarnings(as.numeric(dots$stability_barrier_weight))) &&
                 as.numeric(dots$stability_barrier_weight) > 0) {
      dots$stability_barrier_weight
    } else {
      5000
    },
    power = dots$stability_barrier_power
  )
  biv_barrier_ctrl <- .etas_soft_barrier_controls(
    etas_rho_max,
    start = dots$stability_barrier_start,
    weight = if (!is.null(dots$stability_barrier_weight) &&
                 is.finite(suppressWarnings(as.numeric(dots$stability_barrier_weight))) &&
                 as.numeric(dots$stability_barrier_weight) > 0) {
      dots$stability_barrier_weight
    } else {
      5000
    },
    power = dots$stability_barrier_power
  )
  hard_subcritical <- !isFALSE(dots$hard_subcritical)
  if ((is_etas || is_biv_etas) && !is.finite(etas_beta_eff)) {
    stop("adaptive_SEM requires a finite positive beta_gr for stable ETAS fitting.")
  }
  etas_init_eta_cap <- etas_eta_max * etas_init_margin
  etas_init_rho_cap <- etas_rho_max * etas_init_margin
  if (is_biv_etas && !is.null(dots$etas_bivariate_params)) {
    dots$etas_bivariate_params <- .etas_project_subcritical_biv(
      dots$etas_bivariate_params,
      beta_gr = etas_beta_eff,
      gap_min = etas_gap_min,
      rho_max = etas_init_rho_cap,
      alpha_m_lower_bound = etas_alpha_lo
    )
  }
  if (is_etas) {
    hawkes_params_control <- as.list(.etas_project_subcritical(
      hawkes_params_control,
      beta_gr = etas_beta_eff,
      gap_min = etas_gap_min,
      eta_max = etas_init_eta_cap,
      alpha_m_lower_bound = etas_alpha_lo
    ))
    hawkes_params_treated <- as.list(.etas_project_subcritical(
      hawkes_params_treated,
      beta_gr = etas_beta_eff,
      gap_min = etas_gap_min,
      eta_max = etas_init_eta_cap,
      alpha_m_lower_bound = etas_alpha_lo
    ))
  }
  dots$beta_gr <- etas_beta_eff
  dots$alpha_beta_gap_min <- etas_gap_min
  dots$max_branching_ratio <- etas_eta_max
  dots$max_branching_radius <- etas_rho_max
  dots$alpha_m_lower_bound <- etas_alpha_lo
  dots$hard_subcritical <- hard_subcritical
  dots$init_branching_margin <- etas_init_margin
  dots$soft_branching_barrier <- etas_soft_barrier
  dots$log_transform <- etas_log_transform
  dots_no_trunc <- dots
  dots_no_trunc$t_trunc <- NULL
  hawkes_kernel <- normalize_hawkes_kernel(dots$kernel, hawkes_params_control)
  hawkes_spatial_kernel <- normalize_hawkes_spatial_kernel(dots$spatial_kernel, hawkes_params_control)
  hawkes_spatial_q <- if (!is.null(dots$spatial_q)) dots$spatial_q else hawkes_params_control$spatial_q
  hawkes_spatial_d <- dots$spatial_d
  loglik_fn <- if (is_biv_etas) loglik_etas_bivariate
               else if (is_etas) loglik_etas
               else loglik_hawk_fast

  if (partition$type != "mask") {
    if (verbose) message("Converting partition to raster mask for speed...")
    statespace <- spatstat.geom::as.mask(as.owin(statespace), dimyx = c(512, 512))
    tiles_mask <- lapply(tiles(partition), function(w) spatstat.geom::as.mask(w, xy = statespace))
    partition <- tess(tiles = tiles_mask)
    if (!is.null(adaptive_control$state_spaces)) {
      adaptive_control$state_spaces <- lapply(adaptive_control$state_spaces,
        function(w) spatstat.geom::as.mask(w, xy = statespace))
      if (is.null(names(adaptive_control$state_spaces))) {
        procs <- unique(partition_processes)
        if (length(procs) == length(adaptive_control$state_spaces)) {
          names(adaptive_control$state_spaces) <- procs
        }
      }
    }
  }

  background_rate_var <- if ("background_rate_var" %in% names(dots)) dots$background_rate_var else NULL

  treated_idx <- (partition_processes == "treated")
  treated_state_space <- as.owin(partition[treated_idx])
  control_state_space <- as.owin(partition[!treated_idx])

  if (is.null(pp_data$inferred_process)) {
    pp_data$inferred_process <- pp_data$location_process
  }
  starting_data <- pp_data

  if (is.null(hawkes_params_treated)) {
    hawkes_params_treated <- hawkes_params_control
  }

  ac_defaults <- list(
    update_control_params = FALSE, param_update_cadence = 20,
    proposal_update_cadence = 1, update_starting_data = TRUE,
    include_starting_data = FALSE, include_starting_first_n = 50, iter = 100, n_props = 10,
    max_relabel_step_frac = 1.0, force_param_update_flip_frac = 1.0,
    optim_method = "max", selection_temperature = 0.15,
    change_factor_min_mult = 0.2, change_factor_max_mult = 2.0,
    change_factor = 0.1, stagnation_trigger_every = 10, proposal_method = "simulation",
    temporal_weight = 0, temporal_scale_days = NULL,
    fixed_params = NULL, verbose = FALSE, state_spaces = NULL,
    outer_maxit = 500, outer_maxit_biv = NULL, param_refit_cadence = 1,
    biv_n_threads = 1L
  )
  for (nm in names(ac_defaults)) {
    if (is.null(adaptive_control[[nm]])) adaptive_control[[nm]] <- ac_defaults[[nm]]
  }

  adaptive_step <- function(starting_data) {
    t_adapt_start <- proc.time()[3]
    result <- do.call(
      em_style_labelling,
      c(
        list(
          pp_data = starting_data,
          partition = partition,
          partition_processes = partition_processes,
          statespace = statespace,
          state_spaces = adaptive_control$state_spaces,
          time_window = c(treatment_time, max(starting_data$t)),
          treatment_time = treatment_time,
          hawkes_params_control = hawkes_params_control,
          hawkes_params_treated = hawkes_params_treated,
          update_control_params = adaptive_control$update_control_params,
          param_update_cadence = adaptive_control$param_update_cadence,
          proposal_update_cadence = adaptive_control$proposal_update_cadence,
          update_starting_data = adaptive_control$update_starting_data,
          include_starting_data = adaptive_control$include_starting_data,
          include_starting_first_n = adaptive_control$include_starting_first_n,
          max_relabel_step_frac = adaptive_control$max_relabel_step_frac,
          force_param_update_flip_frac = adaptive_control$force_param_update_flip_frac,
          metric_name = "post_likelihood",
          optim_method = adaptive_control$optim_method,
          selection_temperature = adaptive_control$selection_temperature,
          change_factor_min_mult = adaptive_control$change_factor_min_mult,
          change_factor_max_mult = adaptive_control$change_factor_max_mult,
          iter = adaptive_control$iter,
          n_props = adaptive_control$n_props,
          change_factor = adaptive_control$change_factor,
          stagnation_trigger_every = adaptive_control$stagnation_trigger_every,
          MCMC_style = FALSE,
          proposal_method = adaptive_control$proposal_method,
          temporal_weight = adaptive_control$temporal_weight,
          temporal_scale_days = adaptive_control$temporal_scale_days,
          biv_n_threads = adaptive_control$biv_n_threads,
          fixed_params = adaptive_control$fixed_params,
          # Respect caller-configured SEM verbosity instead of forcing full trace.
          verbose = isTRUE(adaptive_control$verbose),
          model_type = model_type
        ),
        dots
      )
    )
    t_adapt_end <- proc.time()[3]
    return(list(
      adaptive_labelling = result$labelling,
      treated_par = result$treated_par,
      control_par = result$control_par,
      accuracies = result$accuracies,
      average_flips = result$average_flips,
      max_metric_flips = result$max_metric_flips,
      metrics = result$metrics,
      all_accuracies = result$all_accuracies,
      all_metrics = result$all_metrics,
      class_results = result$class_results,
      fits = result$fits,
      etas_bivariate_params = result$etas_bivariate_params,
      etas_bivariate_convergence = result$etas_bivariate_convergence,
      etas_bivariate_value = result$etas_bivariate_value,
      time_taken = t_adapt_end - t_adapt_start
    ))
  }

  weights <- rep(1, N_labellings)
  em_iter <- 0
  counter <- 0
  adaptive_counter <- 0
  adaptive_history <- list()
  baseline_adaptive_labelling <- NULL
  last_biv_optim <- NULL

  pre <- as.data.frame(starting_data) %>% dplyr::filter(.data$t < treatment_time)
  post <- as.data.frame(starting_data) %>% dplyr::filter(.data$t >= treatment_time)
  post <- post[order(post$t), ]
  max_data_t <- max(starting_data$t)
  sem_windowT <- c(treatment_time, max_data_t)
  hawkes_bg_var <- if (!is.null(background_rate_var)) background_rate_var else "W"

  hawkes_loglik_with_filtration <- function(params, post_realiz, filt_realiz, zero_bg_region) {
    par_obj <- as_hawkes_params(params, hawkes_kernel, hawkes_spatial_kernel, hawkes_spatial_q, hawkes_spatial_d)
    mu <- par_obj$mu
    alpha <- par_obj$alpha
    beta <- par_obj$beta
    K <- par_obj$K
    cc <- par_obj[["c"]]
    p <- par_obj$p
    q_spatial <- if (is.null(par_obj$spatial_q)) 2.0 else as.numeric(par_obj$spatial_q)
    d_spatial <- if (is.null(par_obj$spatial_d)) NA_real_ else as.numeric(par_obj$spatial_d)
    if (!is.finite(mu) || !is.finite(alpha) || !is.finite(K)) return(-Inf)
    if (mu < 0 || alpha < 0 || K < 0 || K >= 1) return(-Inf)
    if (identical(hawkes_kernel, "power_law")) {
      if (!is.finite(cc) || !is.finite(p) || cc <= 0 || p <= 1) return(-Inf)
    } else {
      if (!is.finite(beta) || beta <= 0) return(-Inf)
    }
    if (!inherits(statespace, "owin")) statespace <- as.owin(statespace)
    total_area <- spatstat.geom::area(statespace)
    dt <- max_data_t - treatment_time
    if (!is.finite(dt) || dt <= 0) return(-Inf)

    active_area <- total_area
    in_zero <- rep(FALSE, nrow(post_realiz))
    if (!is.null(zero_bg_region) && nrow(post_realiz) > 0) {
      if (!inherits(zero_bg_region, "owin")) zero_bg_region <- as.owin(zero_bg_region)
      zero_area <- spatstat.geom::area(zero_bg_region)
      active_area <- max(1e-12, total_area - zero_area)
      in_zero <- inside.owin(post_realiz$x, post_realiz$y, w = zero_bg_region)
    }

    W_post <- if (!is.null(hawkes_bg_var) && hawkes_bg_var %in% names(post_realiz)) {
      as.numeric(post_realiz[[hawkes_bg_var]])
    } else {
      rep(1, nrow(post_realiz))
    }
    W_post[!is.finite(W_post)] <- 0
    W_post[in_zero] <- 0

    if (!is.null(filt_realiz) && nrow(filt_realiz) > 0) {
      filt_realiz <- as.data.frame(filt_realiz)
      filt_realiz <- filt_realiz[filt_realiz$t < treatment_time, , drop = FALSE]
      filt_realiz <- filt_realiz[order(filt_realiz$t), , drop = FALSE]
    } else {
      filt_realiz <- post_realiz[0, c("x", "y", "t"), drop = FALSE]
    }
    post_realiz <- post_realiz[order(post_realiz$t), , drop = FALSE]

    # Keep SEM filtration consistent with ATE fitting:
    # evaluate conditional post-treatment likelihood, using pre-treatment
    # filtration only as parent history (not as observed events).
    post_t_fit <- as.numeric(post_realiz$t)
    post_x_fit <- as.numeric(post_realiz$x)
    post_y_fit <- as.numeric(post_realiz$y)
    W_fit <- as.numeric(W_post)
    parent_x <- c(filt_realiz$x, post_realiz$x)
    parent_y <- c(filt_realiz$y, post_realiz$y)
    parent_t <- c(filt_realiz$t, post_realiz$t)
    t_start_fit <- treatment_time
    t_end_fit <- max_data_t
    adjust_factor_fit <- 1
    if (length(parent_t) < 1L || length(post_t_fit) < 1L) return(-Inf)

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
      t_trunc = if (!is.null(t_trunc)) t_trunc else -1,
      kernel_type = hawkes_kernel_type(hawkes_kernel),
      cc = if (is.null(cc)) 1.0 else as.numeric(cc),
      p = if (is.null(p)) 2.0 else as.numeric(p),
      spatial_kernel_type = hawkes_spatial_kernel_type(hawkes_spatial_kernel),
      spatial_q = q_spatial,
      spatial_d = d_spatial
    )
    if (!is.finite(loglik)) return(-Inf)
    loglik
  }

  calculate_weights <- function(labellings, treat_par, control_par, ...) {
    if (is_biv_etas) {
      biv_par <- if (is.null(dots$etas_bivariate_params)) {
        init_bivariate_from_independent(control_par, treat_par)
      } else {
        dots$etas_bivariate_params
      }
      biv_wT <- if (use_pre_history_for_biv) {
        c(min(starting_data$t), max_data_t)
      } else {
        sem_windowT
      }
      K_w <- length(labellings)
      if (K_w == 0L) return(numeric(0))
      geom0 <- if (use_pre_history_for_biv) {
        labellings[[1]]
      } else {
        labellings[[1]][labellings[[1]]$t >= treatment_time, , drop = FALSE]
      }
      geom0 <- geom0[order(geom0$t), , drop = FALSE]
      nn <- nrow(geom0)
      aS0 <- spatstat.geom::area(control_state_space)
      aS1 <- spatstat.geom::area(treated_state_space)
      if (aS0 <= 0) aS0 <- 1; if (aS1 <= 0) aS1 <- 1
      W0 <- rep(1.0, nn); W1 <- rep(1.0, nn)
      in_treated_w <- inside.owin(geom0$x, geom0$y, treated_state_space)
      if (is.null(control_background_everywhere_before)) {
        W0[in_treated_w] <- 0
      } else {
        # Control background covers the whole domain before the policy cutoff.
        W0[in_treated_w & geom0$t >= control_background_everywhere_before] <- 0
      }
      W1[inside.owin(geom0$x, geom0$y, control_state_space)] <- 0
      if (!is.null(treated_background_zero_before)) {
        W1[geom0$t < treated_background_zero_before] <- 0
      }
      if (!is.null(background_rate_var) && background_rate_var %in% names(geom0)) {
        W_cov <- as.numeric(geom0[[background_rate_var]])
        W_cov[!is.finite(W_cov)] <- 0
        min_pos <- suppressWarnings(min(W_cov[W_cov > 0], na.rm = TRUE))
        if (!is.finite(min_pos)) min_pos <- 1e-12
        W_cov[W_cov <= 0] <- min_pos
        W0 <- W0 * W_cov; W1 <- W1 * W_cov
      }
      pid_mat <- matrix(0L, nn, K_w)
      for (kk in seq_len(K_w)) {
        r <- if (use_pre_history_for_biv) {
          labellings[[kk]]
        } else {
          labellings[[kk]][labellings[[kk]]$t >= treatment_time, , drop = FALSE]
        }
        r <- r[order(r$t), , drop = FALSE]
        if (nrow(r) != nn) {
          stop("SEM weight labellings must share event geometry for batch likelihood.")
        }
        proc_col <- if ("inferred_process" %in% names(r)) r$inferred_process else r$location_process
        pid_mat[, kk] <- if (is.character(proc_col)) {
          as.integer(proc_col == "treated")
        } else {
          as.integer(proc_col)
        }
      }
      # Keep weight scoring under the same kernel constraints as the M-step /
      # proposal LL (including Oklahoma's relaxed p>1, q>1 regime).
      biv_weight_extra_names <- intersect(
        names(dots),
        c(
          "m0", "beta_gr", "enforce_finite_trigger_moments",
          "p_lower_bound", "q_lower_bound", "finite_moment_soft_width",
          "finite_moment_soft_weight", "finite_moment_soft_power",
          "enforce_alpha_subcritical", "alpha_beta_gap_min",
          "alpha_m_lower_bound",
          "alpha_beta_soft_gap", "alpha_beta_soft_weight",
          "alpha_beta_soft_power", "max_branching_radius",
          "stability_barrier_start", "stability_barrier_weight",
          "stability_barrier_power"
        )
      )
      biv_weight_extra <- dots[biv_weight_extra_names]
      if (is.null(biv_weight_extra$max_branching_radius)) {
        biv_weight_extra$max_branching_radius <- 0.98
      }
      if (is.null(biv_weight_extra$alpha_beta_gap_min)) {
        biv_weight_extra$alpha_beta_gap_min <- 1e-4
      }
      if (is.null(biv_weight_extra$enforce_alpha_subcritical)) {
        biv_weight_extra$enforce_alpha_subcritical <- TRUE
      }
      return(do.call(
        loglik_etas_bivariate_batch,
        c(
          list(
            params = biv_par,
            t = geom0$t - biv_wT[1], x = geom0$x, y = geom0$y, mag = geom0$mag,
            process_ids = pid_mat,
            W0s = matrix(W0, nn, K_w), W1s = matrix(W1, nn, K_w),
            areaS_0 = aS0, areaS_1 = aS1,
            t_max = biv_wT[2] - biv_wT[1],
            windowT = biv_wT,
            treated_background_zero_before = treated_background_zero_before,
            control_background_everywhere_before = control_background_everywhere_before,
            control_background_pre_mass_ratio = control_background_pre_mass_ratio,
            t_trunc = t_trunc, t_already_shifted = TRUE
          ),
          biv_weight_extra
        )
      ))
    }
    if (!is_etas) {
      # Hawkes: batch control and treated filtration likelihoods over shared geometry.
      K_w <- length(labellings)
      if (K_w == 0L) return(numeric(0))
      ctrl_mem <- .hawkes_batch_membership(
        labellings, "control", treatment_time, hawkes_use_filtration_history
      )
      treat_mem <- .hawkes_batch_membership(
        labellings, "treated", treatment_time, hawkes_use_filtration_history
      )
      total_area <- spatstat.geom::area(as.owin(statespace))
      area_ctrl <- total_area
      area_treat <- total_area
      W_ctrl <- rep(1.0, length(ctrl_mem$t))
      W_treat <- rep(1.0, length(treat_mem$t))
      if (!is.null(hawkes_bg_var) && hawkes_bg_var %in% names(ctrl_mem$geom)) {
        W_ctrl <- as.numeric(ctrl_mem$geom[[hawkes_bg_var]])
        W_treat <- as.numeric(treat_mem$geom[[hawkes_bg_var]])
        W_ctrl[!is.finite(W_ctrl)] <- 0
        W_treat[!is.finite(W_treat)] <- 0
      }
      if (!is.null(treated_state_space)) {
        zero_area <- spatstat.geom::area(as.owin(treated_state_space))
        area_ctrl <- max(1e-12, total_area - zero_area)
        W_ctrl[inside.owin(ctrl_mem$x, ctrl_mem$y, w = treated_state_space)] <- 0
      }
      if (!is.null(control_state_space)) {
        zero_area <- spatstat.geom::area(as.owin(control_state_space))
        area_treat <- max(1e-12, total_area - zero_area)
        W_treat[inside.owin(treat_mem$x, treat_mem$y, w = control_state_space)] <- 0
      }
      ctrl_liks <- loglik_hawk_filtration_batch(
        params = control_par,
        t = ctrl_mem$t, x = ctrl_mem$x, y = ctrl_mem$y,
        is_observed = ctrl_mem$is_observed, member = ctrl_mem$member,
        W_val = W_ctrl, areaS = area_ctrl,
        t_start = treatment_time, t_end = max_data_t,
        t_trunc = t_trunc,
        kernel = hawkes_kernel, spatial_kernel = hawkes_spatial_kernel,
        spatial_q = hawkes_spatial_q, spatial_d = hawkes_spatial_d
      )
      treat_liks <- loglik_hawk_filtration_batch(
        params = treat_par,
        t = treat_mem$t, x = treat_mem$x, y = treat_mem$y,
        is_observed = treat_mem$is_observed, member = treat_mem$member,
        W_val = W_treat, areaS = area_treat,
        t_start = treatment_time, t_end = max_data_t,
        t_trunc = t_trunc,
        kernel = hawkes_kernel, spatial_kernel = hawkes_spatial_kernel,
        spatial_q = hawkes_spatial_q, spatial_d = hawkes_spatial_d
      )
      return(ctrl_liks + treat_liks)
    }
    sapply(labellings, function(y) {
      post_idx <- y$t >= treatment_time
      realiz <- y[post_idx, , drop = FALSE]
      history_control <- if (etas_use_filtration_history) {
        y[!post_idx & y$inferred_process == "control",
          c("x", "y", "t", "mag"), drop = FALSE]
      } else {
        NULL
      }
      history_treated <- if (etas_use_filtration_history) {
        y[!post_idx & y$inferred_process == "treated",
          c("x", "y", "t", "mag"), drop = FALSE]
      } else {
        NULL
      }
      include <- which(realiz$inferred_process == "control")
      if (length(include) == 0) return(-Inf)
      control_lik <- loglik_fn(
        params = control_par, realiz = realiz[include, ],
        windowT = sem_windowT,
        windowS = statespace, zero_background_region = treated_state_space,
        t_trunc = t_trunc, history = history_control, ...
      )
      include <- which(realiz$inferred_process == "treated")
      if (length(include) == 0) return(-Inf)
      treat_lik <- loglik_fn(
        params = treat_par, realiz = realiz[include, ],
        windowT = sem_windowT,
        windowS = statespace, zero_background_region = control_state_space,
        t_trunc = t_trunc, history = history_treated, ...
      )
      return(control_lik + treat_lik)
    })
  }

  t_main_sem_start <- proc.time()[3]
  while (counter < N_iter) {
    if (verbose) {
      cat(sprintf("\n--- SEM Outer Iteration %d / %d ---\n", counter + 1, N_iter))
    }
    if (adaptive_counter == 0) {
      if (verbose) cat("[SEM] Initial adaptive step (runs once)...\n")
      t_adapt_start <- proc.time()[3]
      adapt <- adaptive_step(starting_data = starting_data)
      if (verbose) cat(sprintf("[SEM] Adaptive step complete (took %.1fs)\n", proc.time()[3] - t_adapt_start))
      baseline_adaptive_labelling <- adapt$adaptive_labelling
      c_params <- adapt$control_par
      t_params <- adapt$treated_par
      if (is_biv_etas && !is.null(adapt$etas_bivariate_params)) {
        dots$etas_bivariate_params <- adapt$etas_bivariate_params
        last_inner <- if (is.list(adapt$fits) && length(adapt$fits) > 0L) {
          adapt$fits[[length(adapt$fits)]]
        } else {
          NULL
        }
        last_biv_optim <- list(
          convergence = if (!is.null(adapt$etas_bivariate_convergence)) {
            adapt$etas_bivariate_convergence
          } else if (!is.null(last_inner$convergence)) {
            last_inner$convergence
          } else {
            NA_integer_
          },
          value = if (!is.null(adapt$etas_bivariate_value)) {
            adapt$etas_bivariate_value
          } else if (!is.null(last_inner$value)) {
            last_inner$value
          } else {
            NA_real_
          },
          branching_radius = if (!is.null(last_inner$branching_radius)) {
            last_inner$branching_radius
          } else {
            NA_real_
          }
        )
      }
      adaptive_counter <- 1
      adaptive_history[[1]] <- list(
        max_metric_flips = adapt$max_metric_flips,
        average_flips = adapt$average_flips,
        metrics = adapt$metrics
      )
    }
    if (reset && adaptive_counter > 0) {
      if (verbose) cat("[SEM] Resetting outer counter to 0\n")
      counter <- 0
    }
    if (!is.null(baseline_adaptive_labelling)) {
      if (verbose) cat(sprintf("[SEM] Generating %d labellings from baseline...\n", N_labellings))
      t_gen_start <- proc.time()[3]
      # em_style_labelling returns rbind(pre, post). Propose on post only
      # (matching the adaptive inner loop), then reattach pre. Reuse the
      # proposal simulation cache across the N_labellings draws.
      baseline_post <- baseline_adaptive_labelling[
        baseline_adaptive_labelling$t >= treatment_time, , drop = FALSE
      ]
      baseline_post <- baseline_post[order(baseline_post$t), , drop = FALSE]
      post_inds <- as.numeric(tileindex(
        baseline_post$x, baseline_post$y, partition, close.gaps = FALSE
      ))
      pre_for_proposals <- if (
        is_biv_etas ||
        (is_etas && etas_use_filtration_history) ||
        (!is_etas && hawkes_use_filtration_history)
      ) {
        pre
      } else {
        pre[0, , drop = FALSE]
      }
      proposal_filt_labels <- if (!is.null(pre_for_proposals$inferred_process)) {
        pre_for_proposals$inferred_process
      } else if (!is.null(pre_for_proposals$process)) {
        pre_for_proposals$process
      } else {
        pre_for_proposals$location_process
      }
      proposal_filt_by_proc <- if (!is.null(proposal_filt_labels)) {
        split(pre_for_proposals, proposal_filt_labels)
      } else {
        NULL
      }
      proposal_sim_cache <- NULL
      labellings <- vector("list", N_labellings)
      for (i in seq_len(N_labellings)) {
        prop_result <- simulation_labeling_hawkes_hawkes_fast(
          baseline_post,
          partition = partition, partition_process = partition_processes,
          statespace = statespace,
          state_spaces = adaptive_control$state_spaces,
          windowT = c(treatment_time, max(starting_data$t)),
          hawkes_params_control = c_params[[length(c_params)]],
          hawkes_params_treated = t_params[[length(t_params)]],
          change_factor = adaptive_control$change_factor,
          filtration = pre_for_proposals, proximity_weight = 0,
          filt_by_proc = proposal_filt_by_proc,
          verbose = isTRUE(adaptive_control$verbose),
          points_tile_index = post_inds,
          model_type = model_type,
          proposal_sim_cache = proposal_sim_cache,
          return_proposal_sim_cache = TRUE,
          ...
        )
        prop_post <- if (is.list(prop_result) && !is.null(prop_result$data)) {
          if (is.null(proposal_sim_cache)) {
            proposal_sim_cache <- prop_result$proposal_sim_cache
          }
          prop_result$data
        } else {
          prop_result
        }
        labellings[[i]] <- rbind(pre, prop_post)
      }
      # baseline_adaptive_labelling already includes pre; do not rbind again.
      labellings[[length(labellings) + 1]] <- baseline_adaptive_labelling
      if (verbose) cat(sprintf("[SEM] Labelling generation complete (%d proposals + baseline, took %.1fs)\n",
                               N_labellings, proc.time()[3] - t_gen_start))
    }
    
    if (verbose) cat("[SEM] Calculating importance weights...\n")
    t_weights_start <- proc.time()[3]
    raw_weights <- do.call(calculate_weights, c(list(
      labellings = labellings,
      treat_par = unlist(t_params[[length(t_params)]]),
      control_par = unlist(c_params[[length(c_params)]])
    ), dots_no_trunc))
    keepers <- which(raw_weights != 0 & !is.na(raw_weights) & is.finite(raw_weights))
    if (verbose) {
      cat(sprintf("[SEM] Weight calculation complete (took %.1fs)\n", proc.time()[3] - t_weights_start))
      cat(sprintf("  raw weights range: [%.2g, %.2g]\n", 
                  if(length(raw_weights)>0) min(raw_weights, na.rm=T) else NA, 
                  if(length(raw_weights)>0) max(raw_weights, na.rm=T) else NA))
      cat(sprintf("  keepers: %d / %d labellings\n", length(keepers), length(labellings)))
    }
    if (length(keepers) == 0) {
      warning("All weights are zero/NA/Inf — skipping outer optim")
      counter <- counter + 1
      next
    }
    weights <- normalize_weights(raw_weights[keepers])
    if (verbose) {
      cat("  normalized weights: ", paste(signif(weights, 4), collapse = ", "), "\n")
    }

    # Cache process/time splits once per labelling for non-bivariate outer
    # optim. Bivariate path uses the shared-geometry batch kernel instead.
    prepared_labellings <- if (is_biv_etas) {
      NULL
    } else {
      lapply(labellings[keepers], function(y) {
        post_idx <- y$t >= treatment_time
        pre_idx <- !post_idx
        list(
          post_control = y[post_idx & y$inferred_process == "control", , drop = FALSE],
          post_treated = y[post_idx & y$inferred_process == "treated", , drop = FALSE],
          filt_control = y[pre_idx & y$inferred_process == "control", , drop = FALSE],
          filt_treated = y[pre_idx & y$inferred_process == "treated", , drop = FALSE]
        )
      })
    }

    fp <- if (!is.null(adaptive_control$fixed_params)) adaptive_control$fixed_params else NULL
    outer_maxit <- adaptive_control$outer_maxit
    outer_maxit_biv <- if (!is.null(adaptive_control$outer_maxit_biv)) {
      adaptive_control$outer_maxit_biv
    } else {
      outer_maxit
    }
    outer_optim_trace <- if (!is.null(adaptive_control$outer_optim_trace)) {
      suppressWarnings(as.integer(adaptive_control$outer_optim_trace))
    } else if (isTRUE(verbose)) {
      1L
    } else {
      0L
    }
    if (!is.finite(outer_optim_trace) || is.na(outer_optim_trace) || outer_optim_trace < 0L) {
      outer_optim_trace <- 0L
    }
    outer_optim_report <- if (!is.null(adaptive_control$outer_optim_report)) {
      suppressWarnings(as.integer(adaptive_control$outer_optim_report))
    } else {
      10L
    }
    if (!is.finite(outer_optim_report) || is.na(outer_optim_report) || outer_optim_report < 1L) {
      outer_optim_report <- 10L
    }
    param_refit_cadence <- suppressWarnings(as.integer(adaptive_control$param_refit_cadence))
    if (!is.finite(param_refit_cadence) || is.na(param_refit_cadence) || param_refit_cadence < 1L) {
      param_refit_cadence <- 1L
    }
    do_param_refit <- param_refit_cadence <= 1L ||
      counter == 0L ||
      (counter + 1L) >= N_iter ||
      (counter %% param_refit_cadence) == 0L
    if (!do_param_refit) {
      t_params[[length(t_params) + 1L]] <- t_params[[length(t_params)]]
      if (isTRUE(adaptive_control$update_control_params)) {
        c_params[[length(c_params) + 1L]] <- c_params[[length(c_params)]]
      }
      counter <- counter + 1L
      next
    }

    if (is_biv_etas) {
      # --- Bivariate ETAS: joint optimization over 15-parameter vector ---
      biv_names <- .etas_bivariate_par_names
      biv_fp_idx <- if (!is.null(fp)) match(names(fp), biv_names) else integer(0)
      biv_fp_idx <- biv_fp_idx[!is.na(biv_fp_idx)]
      biv_fr_idx <- setdiff(seq_along(biv_names), biv_fp_idx)

      biv_par <- if (!is.null(dots$etas_bivariate_params)) {
        unlist(dots$etas_bivariate_params)
      } else {
        init_bivariate_from_independent(
          c_params[[length(c_params)]], t_params[[length(t_params)]])
      }
      if (is.null(names(biv_par))) names(biv_par) <- biv_names
      if (!is.null(fp)) {
        for (nm in intersect(names(fp), biv_names)) biv_par[[nm]] <- fp[[nm]]
      }
      biv_par <- .etas_project_subcritical_biv(
        biv_par,
        beta_gr = etas_beta_eff,
        gap_min = etas_gap_min,
        rho_max = etas_init_rho_cap
      )
      biv_free_names <- biv_names[biv_fr_idx]
      biv_alpha_free_pos <- if (hard_subcritical) {
        which(biv_free_names %in% .etas_biv_alpha_names)
      } else {
        integer(0)
      }
      biv_log_names <- c("mu_0", "mu_1", .etas_biv_A_names, "c", "D")
      biv_log_free_pos <- if (etas_log_transform) {
        which(biv_free_names %in% biv_log_names)
      } else {
        integer(0)
      }
      biv_to_optim <- function(natural_free) {
        optim_free <- natural_free
        if (length(biv_log_free_pos)) {
          optim_free[biv_log_free_pos] <- log(pmax(
            as.numeric(natural_free[biv_log_free_pos]), 1e-12
          ))
        }
        if (length(biv_alpha_free_pos)) {
          slack <- pmax(
            etas_beta_eff - etas_gap_min - natural_free[biv_alpha_free_pos],
            1e-12
          )
          optim_free[biv_alpha_free_pos] <- log(slack)
        }
        optim_free
      }
      biv_from_optim <- function(optim_free) {
        natural_free <- optim_free
        if (length(biv_log_free_pos)) {
          natural_free[biv_log_free_pos] <- exp(optim_free[biv_log_free_pos])
        }
        if (length(biv_alpha_free_pos)) {
          natural_free[biv_alpha_free_pos] <-
            etas_beta_eff - etas_gap_min - exp(optim_free[biv_alpha_free_pos])
        }
        natural_free
      }

      biv_wT <- if (use_pre_history_for_biv) {
        c(min(starting_data$t), max_data_t)
      } else {
        sem_windowT
      }
      area_control <- spatstat.geom::area(control_state_space)
      area_treated <- spatstat.geom::area(treated_state_space)
      # Shared geometry across SEM labellings (only labels differ). Build W
      # masks once, then a process_id matrix for the batched kernel.
      keeper_labs <- labellings[keepers]
      geom0 <- if (use_pre_history_for_biv) {
        keeper_labs[[1]]
      } else {
        keeper_labs[[1]][keeper_labs[[1]]$t >= treatment_time, , drop = FALSE]
      }
      geom0 <- geom0[order(geom0$t), , drop = FALSE]
      nn <- nrow(geom0)
      aS0 <- area_control; aS1 <- area_treated
      if (aS0 <= 0) aS0 <- 1; if (aS1 <= 0) aS1 <- 1
      W0_shared <- rep(1.0, nn)
      W1_shared <- rep(1.0, nn)
      in_treated_m <- inside.owin(geom0$x, geom0$y, treated_state_space)
      if (is.null(control_background_everywhere_before)) {
        W0_shared[in_treated_m] <- 0
      } else {
        # Control background covers the whole domain before the policy cutoff.
        W0_shared[in_treated_m & geom0$t >= control_background_everywhere_before] <- 0
      }
      W1_shared[inside.owin(geom0$x, geom0$y, control_state_space)] <- 0
      if (!is.null(treated_background_zero_before)) {
        W1_shared[geom0$t < treated_background_zero_before] <- 0
      }
      if (!is.null(background_rate_var) && background_rate_var %in% names(geom0)) {
        W_cov <- as.numeric(geom0[[background_rate_var]])
        if (length(W_cov) != nn) stop("background_rate_var length mismatch in SEM bivariate precomp.")
        W_cov[!is.finite(W_cov)] <- 0
        min_pos <- suppressWarnings(min(W_cov[W_cov > 0], na.rm = TRUE))
        if (!is.finite(min_pos)) min_pos <- 1e-12
        W_cov[W_cov <= 0] <- min_pos
        W0_shared <- W0_shared * W_cov
        W1_shared <- W1_shared * W_cov
      }
      K_keep <- length(keeper_labs)
      pid_mat <- matrix(0L, nn, K_keep)
      for (kk in seq_len(K_keep)) {
        r <- if (use_pre_history_for_biv) {
          keeper_labs[[kk]]
        } else {
          keeper_labs[[kk]][keeper_labs[[kk]]$t >= treatment_time, , drop = FALSE]
        }
        r <- r[order(r$t), , drop = FALSE]
        if (nrow(r) != nn) {
          stop("SEM bivariate labellings must share event geometry for batch likelihood.")
        }
        proc_col <- if ("inferred_process" %in% names(r)) r$inferred_process else r$location_process
        pid_mat[, kk] <- if (is.character(proc_col)) {
          as.integer(proc_col == "treated")
        } else {
          as.integer(proc_col)
        }
      }
      W0_mat <- matrix(W0_shared, nn, K_keep)
      W1_mat <- matrix(W1_shared, nn, K_keep)
      tt_shift <- geom0$t - biv_wT[1]
      t_max_biv <- biv_wT[2] - biv_wT[1]
      biv_ll_extra_names <- intersect(
        names(dots),
        c(
          "m0", "beta_gr", "enforce_finite_trigger_moments",
          "p_lower_bound", "q_lower_bound", "finite_moment_soft_width",
          "finite_moment_soft_weight", "finite_moment_soft_power",
          "enforce_alpha_subcritical", "alpha_beta_gap_min",
          "alpha_m_lower_bound",
          "alpha_beta_soft_gap", "alpha_beta_soft_weight",
          "alpha_beta_soft_power", "max_branching_radius",
          "stability_barrier_start", "stability_barrier_weight",
          "stability_barrier_power"
        )
      )
      biv_ll_extra <- dots[biv_ll_extra_names]
      if (etas_soft_barrier) {
        biv_ll_extra$max_branching_radius <- Inf
        biv_ll_extra$stability_barrier_start <- biv_barrier_ctrl$stability_barrier_start
        biv_ll_extra$stability_barrier_weight <- biv_barrier_ctrl$stability_barrier_weight
        biv_ll_extra$stability_barrier_power <- biv_barrier_ctrl$stability_barrier_power
      }
      biv_n_threads <- if (!is.null(adaptive_control$biv_n_threads)) {
        as.integer(adaptive_control$biv_n_threads)
      } else {
        1L
      }
      if (!is.finite(biv_n_threads) || biv_n_threads < 1L) biv_n_threads <- 1L
      biv_obj <- function(par15) {
        liks <- do.call(
          loglik_etas_bivariate_batch,
          c(
            list(
              params = par15,
              t = tt_shift, x = geom0$x, y = geom0$y, mag = geom0$mag,
              process_ids = pid_mat, W0s = W0_mat, W1s = W1_mat,
              areaS_0 = aS0, areaS_1 = aS1, t_max = t_max_biv,
              windowT = biv_wT,
              treated_background_zero_before = treated_background_zero_before,
              control_background_everywhere_before = control_background_everywhere_before,
              control_background_pre_mass_ratio = control_background_pre_mass_ratio,
              t_trunc = t_trunc, t_already_shifted = TRUE,
              n_threads = biv_n_threads
            ),
            biv_ll_extra
          )
        )
        sum(liks * weights)
      }

      if (verbose) {
        cat("\n--- Outer optim: bivariate ETAS (joint) ---\n")
        cat(sprintf("  starting par: %s\n",
                    paste(biv_names, signif(biv_par, 5), sep = "=", collapse = "  ")))
        cat(sprintf("  optim control: trace=%d report=%d maxit=%d\n",
                    outer_optim_trace, outer_optim_report, outer_maxit_biv))
      }

      t0 <- proc.time()[3]
      run_biv_nm <- function(start_par) {
        par_local <- start_par
        names(par_local) <- biv_names
        if (length(biv_fr_idx) == 0) {
          return(list(
            res = list(par = numeric(0), convergence = 0),
            par = par_local
          ))
        }
        biv_wrap <- function(optim_free) {
          p15 <- par_local
          p15[biv_fr_idx] <- biv_from_optim(optim_free)
          biv_obj(p15)
        }
        biv_opt_start <- biv_to_optim(par_local[biv_fr_idx])
        res_local <- tryCatch(
          optim(par = biv_opt_start, fn = biv_wrap, method = "Nelder-Mead",
                control = list(
                  fnscale = -1,
                  trace = outer_optim_trace,
                  REPORT = outer_optim_report,
                  maxit = outer_maxit_biv
                )),
          error = function(e) {
            cat(sprintf("  [bivariate] OPTIM ERROR: %s\n", e$message))
            list(par = biv_opt_start, convergence = -99)
          }
        )
        par_local[biv_fr_idx] <- biv_from_optim(res_local$par)
        list(res = res_local, par = par_local)
      }
      nm_biv <- run_biv_nm(biv_par)
      biv_res <- nm_biv$res
      biv_par <- nm_biv$par
      names(biv_par) <- biv_names
      if (hard_subcritical) {
        biv_par <- .etas_project_subcritical_biv(
          biv_par,
          beta_gr = etas_beta_eff,
          gap_min = etas_gap_min,
          rho_max = etas_rho_max
        )
      }
      rho1 <- .etas_biv_spectral_radius(biv_par, etas_beta_eff)
      near_cap_biv <- is.finite(rho1) && rho1 >= (etas_rho_max * etas_near_cap_frac)
      distrust_biv <- !is.null(biv_res$convergence) &&
        (biv_res$convergence != 0L || near_cap_biv)
      if (near_cap_biv || distrust_biv) {
        restart_rho <- min(0.70, etas_init_rho_cap)
        restart_start <- .etas_scale_to_target_rho(
          biv_par, etas_beta_eff, restart_rho
        )
        restart_start <- .etas_project_subcritical_biv(
          restart_start,
          beta_gr = etas_beta_eff,
          gap_min = etas_gap_min,
          rho_max = etas_init_rho_cap
        )
        nm_biv2 <- run_biv_nm(restart_start)
        out2 <- nm_biv2$par
        if (hard_subcritical) {
          out2 <- .etas_project_subcritical_biv(
            out2,
            beta_gr = etas_beta_eff,
            gap_min = etas_gap_min,
            rho_max = etas_rho_max
          )
        }
        if (is.finite(biv_obj(out2)) && biv_obj(out2) > biv_obj(biv_par) + 1e-8) {
          biv_res <- nm_biv2$res
          biv_par <- out2
          if (verbose) {
            cat(sprintf(
              "  [bivariate] interior restart improved objective (rho %.3f -> %.3f)\n",
              rho1, .etas_biv_spectral_radius(biv_par, etas_beta_eff)
            ))
          }
        }
      }
      polish_hi <- etas_rho_max * (1 - 1e-6)
      polish_obj <- function(rho) {
        q <- .etas_scale_to_target_rho(biv_par, etas_beta_eff, rho)
        biv_obj(q)
      }
      polish <- tryCatch(
        stats::optimize(polish_obj, interval = c(0.05, polish_hi), maximum = TRUE),
        error = function(e) NULL
      )
      if (!is.null(polish) && is.finite(polish$objective) &&
          polish$objective > biv_obj(biv_par) + 1e-8) {
        biv_par <- .etas_scale_to_target_rho(biv_par, etas_beta_eff, polish$maximum)
        if (hard_subcritical) {
          biv_par <- .etas_project_subcritical_biv(
            biv_par,
            beta_gr = etas_beta_eff,
            gap_min = etas_gap_min,
            rho_max = etas_rho_max
          )
        }
      }
      biv_res$par <- biv_par
      biv_res$branching_radius <- .etas_biv_spectral_radius(
        biv_par, etas_beta_eff
      )
      dots$etas_bivariate_params <- biv_par
      last_biv_optim <- list(
        convergence = if (!is.null(biv_res$convergence)) {
          as.integer(biv_res$convergence)[1L]
        } else {
          NA_integer_
        },
        value = suppressWarnings(as.numeric(biv_res$value)[1L]),
        branching_radius = biv_res$branching_radius
      )

      # Extract marginal params for downstream compatibility
      t_params[[length(t_params) + 1]] <- as.list(c(
        mu = biv_par[["mu_1"]], A = biv_par[["A_11"]],
        alpha_m = biv_par[["alpha_m_11"]],
        c = biv_par[["c"]], p = biv_par[["p"]],
        D = biv_par[["D"]], gamma = biv_par[["gamma"]], q = biv_par[["q"]]))
      c_params[[length(c_params) + 1]] <- as.list(c(
        mu = biv_par[["mu_0"]], A = biv_par[["A_00"]],
        alpha_m = biv_par[["alpha_m_00"]],
        c = biv_par[["c"]], p = biv_par[["p"]],
        D = biv_par[["D"]], gamma = biv_par[["gamma"]], q = biv_par[["q"]]))

      if (verbose) {
        cat(sprintf("  [bivariate] convergence: %s\n", biv_res$convergence))
        cat(sprintf("  [bivariate] final par: %s\n",
                    paste(biv_names, signif(biv_par, 5), sep = "=", collapse = "  ")))
        cat(sprintf("  [bivariate] took %.1fs\n", proc.time()[3] - t0))
      }
    } else {
    # --- Independent models: separate control/treated optimization ---
    all_names <- if (is_etas) .etas_par_names else c("mu", "alpha", "beta", "K")
    fp_idx <- if (!is.null(fp)) match(names(fp), all_names) else integer(0)
    fr_idx <- setdiff(seq_along(all_names), fp_idx)

    run_outer_optim <- function(process_label, zero_bg_region, par_list) {
      split_key <- if (process_label == "control") "post_control" else "post_treated"
      hawkes_batch <- NULL
      prepared_for_process <- NULL
      if (!is_etas) {
        keeper_labs <- labellings[keepers]
        mem <- .hawkes_batch_membership(
          keeper_labs, process_label, treatment_time, hawkes_use_filtration_history
        )
        total_area <- spatstat.geom::area(as.owin(statespace))
        active_area <- total_area
        W_val <- rep(1.0, length(mem$t))
        if (!is.null(hawkes_bg_var) && hawkes_bg_var %in% names(mem$geom)) {
          W_val <- as.numeric(mem$geom[[hawkes_bg_var]])
          W_val[!is.finite(W_val)] <- 0
        }
        if (!is.null(zero_bg_region)) {
          zero_area <- spatstat.geom::area(as.owin(zero_bg_region))
          active_area <- max(1e-12, total_area - zero_area)
          W_val[inside.owin(mem$x, mem$y, w = zero_bg_region)] <- 0
        }
        hawkes_n_threads <- if (!is.null(adaptive_control$hawkes_n_threads)) {
          as.integer(adaptive_control$hawkes_n_threads)
        } else {
          1L
        }
        if (!is.finite(hawkes_n_threads) || hawkes_n_threads < 1L) hawkes_n_threads <- 1L
        hawkes_batch <- list(
          mem = mem, W_val = W_val, active_area = active_area,
          n_threads = hawkes_n_threads
        )
      } else {
        prepared_for_process <- lapply(prepared_labellings, function(parts) {
          post_part <- parts[[split_key]]
          history_key <- if (process_label == "control") {
            "filt_control"
          } else {
            "filt_treated"
          }
          pc <- precompute_loglik_args(post_part, statespace, zero_bg_region)
          list(
            post_part = post_part,
            history = if (etas_use_filtration_history) {
              parts[[history_key]][, c("x", "y", "t", "mag"), drop = FALSE]
            } else {
              NULL
            },
            precomp = list(
              active_area = pc$active_area,
              in_zero_bg = pc$in_zero_bg_all
            )
          )
        })
      }
      etas_ll_extra_names <- if (is_etas) {
        intersect(
          names(dots),
          c(
            "m0", "beta_gr", "enforce_finite_trigger_moments",
            "p_lower_bound", "q_lower_bound", "finite_moment_soft_width",
            "finite_moment_soft_weight", "finite_moment_soft_power",
            "enforce_alpha_subcritical", "alpha_beta_gap_min",
            "alpha_m_lower_bound",
            "alpha_beta_soft_gap", "alpha_beta_soft_weight",
            "alpha_beta_soft_power", "max_branching_ratio",
            "stability_barrier_start", "stability_barrier_weight",
            "stability_barrier_power"
          )
        )
      } else {
        character(0)
      }
      etas_ll_extra <- dots[etas_ll_extra_names]
      if (is_etas && etas_soft_barrier) {
        etas_ll_extra$max_branching_ratio <- Inf
        etas_ll_extra$stability_barrier_start <- univ_barrier_ctrl$stability_barrier_start
        etas_ll_extra$stability_barrier_weight <- univ_barrier_ctrl$stability_barrier_weight
        etas_ll_extra$stability_barrier_power <- univ_barrier_ctrl$stability_barrier_power
      }
      obj_fn <- function(params) {
        if (!is_etas) {
          liks <- loglik_hawk_filtration_batch(
            params = params,
            t = hawkes_batch$mem$t,
            x = hawkes_batch$mem$x,
            y = hawkes_batch$mem$y,
            is_observed = hawkes_batch$mem$is_observed,
            member = hawkes_batch$mem$member,
            W_val = hawkes_batch$W_val,
            areaS = hawkes_batch$active_area,
            t_start = treatment_time,
            t_end = max_data_t,
            t_trunc = t_trunc,
            kernel = hawkes_kernel,
            spatial_kernel = hawkes_spatial_kernel,
            spatial_q = hawkes_spatial_q,
            spatial_d = hawkes_spatial_d,
            n_threads = hawkes_batch$n_threads
          )
          return(sum(liks * weights))
        }
        liks <- sapply(prepared_for_process, function(parts) {
          if (is.null(parts)) return(-Inf)
          do.call(
            loglik_fn,
            c(
              list(
                params = params,
                realiz = parts$post_part,
                windowT = sem_windowT,
                windowS = statespace,
                zero_background_region = zero_bg_region,
                background_rate_var = "W",
                precomp = parts$precomp,
                t_trunc = t_trunc,
                history = parts$history
              ),
              etas_ll_extra
            )
          )
        })
        return(sum(liks * weights))
      }

      full_vec <- unlist(par_list)
      if (!is.null(fp)) {
        for (nm in intersect(names(fp), all_names)) full_vec[[nm]] <- fp[[nm]]
      }
      if (is_etas) {
        full_vec <- .etas_project_subcritical(
          full_vec,
          beta_gr = etas_beta_eff,
          gap_min = etas_gap_min,
          eta_max = etas_init_eta_cap
        )
      }
      free_names <- all_names[fr_idx]
      alpha_free_pos <- if (is_etas && hard_subcritical) {
        which(free_names == "alpha_m")
      } else {
        integer(0)
      }
      log_free_pos <- if (is_etas && etas_log_transform) {
        which(free_names %in% c("mu", "A", "c", "D"))
      } else {
        integer(0)
      }
      to_optim <- function(natural_free) {
        optim_free <- natural_free
        if (length(log_free_pos)) {
          optim_free[log_free_pos] <- log(pmax(
            as.numeric(natural_free[log_free_pos]), 1e-12
          ))
        }
        if (length(alpha_free_pos)) {
          slack <- pmax(
            etas_beta_eff - etas_gap_min - natural_free[alpha_free_pos],
            1e-12
          )
          optim_free[alpha_free_pos] <- log(slack)
        }
        optim_free
      }
      from_optim <- function(optim_free) {
        natural_free <- optim_free
        if (length(log_free_pos)) {
          natural_free[log_free_pos] <- exp(optim_free[log_free_pos])
        }
        if (length(alpha_free_pos)) {
          natural_free[alpha_free_pos] <-
            etas_beta_eff - etas_gap_min - exp(optim_free[alpha_free_pos])
        }
        natural_free
      }
      if (verbose) {
        cat(sprintf("  [%s] starting par: %s\n", process_label,
                    paste(all_names, signif(full_vec, 5), sep = "=", collapse = "  ")))
        cat(sprintf("  [%s] starting lik: %s\n", process_label, signif(obj_fn(full_vec), 6)))
        cat(sprintf("  [%s] optim control: trace=%d report=%d maxit=%d\n",
                    process_label, outer_optim_trace, outer_optim_report, outer_maxit))
      }

      t0 <- proc.time()[3]
      run_outer_nm <- function(start_natural) {
        if (length(fr_idx) == 0) {
          return(list(
            res = list(par = numeric(0), convergence = 0),
            out_par = start_natural
          ))
        }
        opt_start <- to_optim(start_natural[fr_idx])
        wrap_fn <- function(optim_free) {
          p4 <- start_natural
          p4[fr_idx] <- from_optim(optim_free)
          obj_fn(p4)
        }
        res_local <- tryCatch(
          optim(par = opt_start, fn = wrap_fn, method = "Nelder-Mead",
                control = list(
                  fnscale = -1,
                  trace = outer_optim_trace,
                  REPORT = outer_optim_report,
                  maxit = outer_maxit
                )),
          error = function(e) {
            cat(sprintf("  [%s] OPTIM ERROR: %s\n", process_label, e$message))
            list(par = opt_start, convergence = -99)
          }
        )
        out_local <- start_natural
        out_local[fr_idx] <- from_optim(res_local$par)
        list(res = res_local, out_par = out_local)
      }
      nm1 <- run_outer_nm(full_vec)
      res <- nm1$res
      out_par <- nm1$out_par
      if (is_etas && hard_subcritical) {
        out_par <- .etas_project_subcritical(
          out_par,
          beta_gr = etas_beta_eff,
          gap_min = etas_gap_min,
          eta_max = etas_eta_max
        )
      }
      if (is_etas && "A" %in% free_names) {
        br1 <- .etas_univ_branching_ratio(out_par, etas_beta_eff)
        near_cap <- is.finite(br1) && br1 >= (etas_eta_max * etas_near_cap_frac)
        distrust <- !is.null(res$convergence) &&
          (res$convergence != 0L || near_cap)
        if (near_cap || distrust) {
          restart_eta <- min(0.70, etas_init_eta_cap)
          restart_start <- .etas_scale_to_target_eta(
            out_par, etas_beta_eff, restart_eta
          )
          restart_start <- .etas_project_subcritical(
            restart_start,
            beta_gr = etas_beta_eff,
            gap_min = etas_gap_min,
            eta_max = etas_init_eta_cap
          )
          nm2 <- run_outer_nm(restart_start)
          out2 <- nm2$out_par
          if (hard_subcritical) {
            out2 <- .etas_project_subcritical(
              out2,
              beta_gr = etas_beta_eff,
              gap_min = etas_gap_min,
              eta_max = etas_eta_max
            )
          }
          if (is.finite(obj_fn(out2)) && obj_fn(out2) > obj_fn(out_par) + 1e-8) {
            res <- nm2$res
            out_par <- out2
            if (verbose) {
              cat(sprintf(
                "  [%s] interior restart improved outer objective (eta %.3f -> %.3f)\n",
                process_label, br1,
                .etas_univ_branching_ratio(out_par, etas_beta_eff)
              ))
            }
          }
        }
        # 1D productivity polish on the weighted SEM objective.
        polish_hi <- etas_eta_max * (1 - 1e-6)
        polish_obj <- function(eta) {
          q <- .etas_scale_to_target_eta(out_par, etas_beta_eff, eta)
          obj_fn(q)
        }
        polish <- tryCatch(
          stats::optimize(polish_obj, interval = c(0.05, polish_hi), maximum = TRUE),
          error = function(e) NULL
        )
        if (!is.null(polish) && is.finite(polish$objective) &&
            polish$objective > obj_fn(out_par) + 1e-8) {
          out_par <- .etas_scale_to_target_eta(out_par, etas_beta_eff, polish$maximum)
          if (hard_subcritical) {
            out_par <- .etas_project_subcritical(
              out_par,
              beta_gr = etas_beta_eff,
              gap_min = etas_gap_min,
              eta_max = etas_eta_max
            )
          }
        }
      }
      res$par <- out_par
      if (is_etas) {
        res$branching_ratio <- .etas_univ_branching_ratio(
          res$par, etas_beta_eff
        )
      }

      if (verbose) {
        cat(sprintf("  [%s] convergence: %s\n", process_label,
                    if (!is.null(res$convergence)) res$convergence else "NULL"))
        cat(sprintf("  [%s] final par:   %s\n", process_label,
                    paste(all_names, signif(as.numeric(res$par), 5), sep = "=", collapse = "  ")))
        cat(sprintf("  [%s] final lik: %s  (took %.1fs)\n", process_label,
                    signif(obj_fn(res$par), 6), proc.time()[3] - t0))
      }
      return(res$par)
    }

    if (verbose) cat("\n--- Outer optim: treated ---\n")
    fit_t <- as.list(run_outer_optim("treated", control_state_space, t_params[[length(t_params)]]))
    t_params[[length(t_params) + 1]] <- fit_t

    if (isTRUE(adaptive_control$update_control_params)) {
      if (verbose) cat("\n--- Outer optim: control ---\n")
      fit_c <- as.list(run_outer_optim("control", treated_state_space, c_params[[length(c_params)]]))
      c_params[[length(c_params) + 1]] <- fit_c
    }
    } # end non-bivariate branch

    counter <- counter + 1
  }
  t_main_sem_end <- proc.time()[3]

  print(paste0("Total time taken for SEM: ", signif(proc.time()[3] - t_global, 2)))
  adaptive_summary <- list(
    accuracies = adapt$accuracies,
    average_flips = adapt$average_flips,
    max_metric_flips = adapt$max_metric_flips,
    metrics = adapt$metrics,
    all_metrics = adapt$all_metrics,
    class_results = adapt$class_results,
    adaptive_labelling = if (!is.null(baseline_adaptive_labelling)) baseline_adaptive_labelling else adapt$adaptive_labelling
  )
  result <- list(
    hawkes_params_control = c_params[[length(c_params)]],
    hawkes_params_treated = t_params[[length(t_params)]],
    c_params = c_params,
    t_params = t_params,
    labellings = labellings,
    keepers = keepers,
    weights = weights,
    adaptive = adaptive_summary,
    adaptive_history = adaptive_history,
    time = proc.time()[3] - t_global,
    time_main_sem = t_main_sem_end - t_main_sem_start
  )
  if (is_biv_etas && !is.null(dots$etas_bivariate_params)) {
    result$etas_bivariate_params <- dots$etas_bivariate_params
    result$etas_bivariate_convergence <- if (!is.null(last_biv_optim)) {
      suppressWarnings(as.integer(last_biv_optim$convergence)[1L])
    } else {
      NA_integer_
    }
    result$etas_bivariate_value <- if (!is.null(last_biv_optim)) {
      suppressWarnings(as.numeric(last_biv_optim$value)[1L])
    } else {
      NA_real_
    }
    result$stability <- list(
      law = "bivariate",
      beta_gr = etas_beta_eff,
      threshold = etas_rho_max,
      branching_metric = .etas_biv_spectral_radius(
        dots$etas_bivariate_params, etas_beta_eff
      )
    )
  } else if (is_etas) {
    result$stability <- list(
      law = "univariate",
      beta_gr = etas_beta_eff,
      threshold = etas_eta_max,
      eta_control = .etas_univ_branching_ratio(
        result$hawkes_params_control, etas_beta_eff
      ),
      eta_treated = .etas_univ_branching_ratio(
        result$hawkes_params_treated, etas_beta_eff
      )
    )
  }
  return(result)
}

#' Check convergence of the SEM algorithm
#'
#' Plots the proportion of points in each labeling category over iterations.
#'
#' @param result Output from adaptive_SEM
#' @return Data frame of proportions per iteration
#' @export
convergence_check <- function(result) {
  check <- do.call(rbind, lapply(result$adaptive$class_results, as.data.frame))
  check$total <- rowSums(check)
  for (i in names(check)) {
    check[[i]] <- check[[i]] / check$total
  }
  check$iter <- 1:nrow(check)
  check <- check[, which(names(check) != "total")]
  check_long <- check %>% reshape2::melt(id.vars = "iter")
  p <- ggplot2::ggplot(check_long, ggplot2::aes(x = .data$iter, y = .data$value, color = .data$variable)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      title = "Proportion of Points in Each Category Over SEM Iterations",
      x = "Iteration", y = "Proportion"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  print(p)
  return(check)
}
