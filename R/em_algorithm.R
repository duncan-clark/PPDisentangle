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
  treated_background_zero_before <- if ("treated_background_zero_before" %in% names(dots)) {
    as.numeric(dots$treated_background_zero_before)
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
  hard_subcritical <- !isFALSE(dots$hard_subcritical)
  if ((is_etas || is_biv_etas) && !is.finite(etas_beta_eff)) {
    stop("adaptive_SEM requires a finite positive beta_gr for stable ETAS fitting.")
  }
  if (is_biv_etas && !is.null(dots$etas_bivariate_params)) {
    dots$etas_bivariate_params <- .etas_project_subcritical_biv(
      dots$etas_bivariate_params,
      beta_gr = etas_beta_eff,
      gap_min = etas_gap_min,
      rho_max = etas_rho_max
    )
  }
  if (is_etas) {
    hawkes_params_control <- as.list(.etas_project_subcritical(
      hawkes_params_control,
      beta_gr = etas_beta_eff,
      gap_min = etas_gap_min,
      eta_max = etas_eta_max
    ))
    hawkes_params_treated <- as.list(.etas_project_subcritical(
      hawkes_params_treated,
      beta_gr = etas_beta_eff,
      gap_min = etas_gap_min,
      eta_max = etas_eta_max
    ))
  }
  dots$beta_gr <- etas_beta_eff
  dots$alpha_beta_gap_min <- etas_gap_min
  dots$max_branching_ratio <- etas_eta_max
  dots$max_branching_radius <- etas_rho_max
  dots$hard_subcritical <- hard_subcritical
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
    outer_maxit = 500, outer_maxit_biv = NULL, param_refit_cadence = 1
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
      time_taken = t_adapt_end - t_adapt_start
    ))
  }

  weights <- rep(1, N_labellings)
  em_iter <- 0
  counter <- 0
  adaptive_counter <- 0
  adaptive_history <- list()
  baseline_adaptive_labelling <- NULL

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
    sapply(labellings, function(y) {
      post_idx <- y$t >= treatment_time
      pre_idx <- !post_idx
      realiz <- if (is_biv_etas && use_pre_history_for_biv) {
        y
      } else if (is_etas) {
        y[post_idx, , drop = FALSE]
      } else {
        y
      }
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
        return(loglik_etas_bivariate(
          params = biv_par, realiz = realiz,
          windowT = biv_wT,
          windowS = statespace,
          control_state_space = control_state_space,
          treated_state_space = treated_state_space,
          background_rate_var = background_rate_var,
          treated_background_zero_before = treated_background_zero_before,
          t_trunc = t_trunc,
          m0 = dots$m0,
          beta_gr = dots$beta_gr,
          max_branching_radius = if (!is.null(dots$max_branching_radius)) dots$max_branching_radius else 0.98,
          alpha_beta_gap_min = if (!is.null(dots$alpha_beta_gap_min)) dots$alpha_beta_gap_min else 1e-4,
          enforce_alpha_subcritical = if (!is.null(dots$enforce_alpha_subcritical)) dots$enforce_alpha_subcritical else TRUE
        ))
      }
      include <- if (is_etas) {
        which(realiz$inferred_process == "control")
      } else {
        which(post_idx & y$inferred_process == "control")
      }
      if (length(include) == 0) return(-Inf)
      if (!is_etas) {
        filt_control <- if (hawkes_use_filtration_history) {
          y[pre_idx & y$inferred_process == "control", , drop = FALSE]
        } else {
          y[0, c("x", "y", "t"), drop = FALSE]
        }
        control_lik <- hawkes_loglik_with_filtration(
          params = control_par,
          post_realiz = y[include, , drop = FALSE],
          filt_realiz = filt_control,
          zero_bg_region = treated_state_space
        )
      } else {
        control_lik <- loglik_fn(
          params = control_par, realiz = realiz[include, ],
          windowT = sem_windowT,
          windowS = statespace, zero_background_region = treated_state_space,
          t_trunc = t_trunc, ...
        )
      }
      include <- if (is_etas) {
        which(realiz$inferred_process == "treated")
      } else {
        which(post_idx & y$inferred_process == "treated")
      }
      if (length(include) == 0) return(-Inf)
      if (!is_etas) {
        filt_treated <- if (hawkes_use_filtration_history) {
          y[pre_idx & y$inferred_process == "treated", , drop = FALSE]
        } else {
          y[0, c("x", "y", "t"), drop = FALSE]
        }
        treat_lik <- hawkes_loglik_with_filtration(
          params = treat_par,
          post_realiz = y[include, , drop = FALSE],
          filt_realiz = filt_treated,
          zero_bg_region = control_state_space
        )
      } else {
        treat_lik <- loglik_fn(
          params = treat_par, realiz = realiz[include, ],
          windowT = sem_windowT,
          windowS = statespace, zero_background_region = control_state_space,
          t_trunc = t_trunc, ...
        )
      }
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
      labellings <- lapply(seq_len(N_labellings), function(i) {
        simulation_labeling_hawkes_hawkes_fast(
          baseline_adaptive_labelling,
          partition = partition, partition_process = partition_processes,
          statespace = statespace,
          state_spaces = adaptive_control$state_spaces,
          windowT = c(treatment_time, max(starting_data$t)),
          hawkes_params_control = hawkes_params_control,
          hawkes_params_treated = t_params[[length(t_params)]],
          change_factor = adaptive_control$change_factor,
          filtration = pre, proximity_weight = 0,
          verbose = isTRUE(adaptive_control$verbose),
          model_type = model_type, ...
        )
      })
      baseline_with_pre <- rbind(pre, baseline_adaptive_labelling)
      labellings[[length(labellings) + 1]] <- baseline_with_pre
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

    # Cache process/time splits once per labelling to avoid repeated filtering
    # inside objective function evaluations.
    prepared_labellings <- lapply(labellings[keepers], function(y) {
      post_idx <- y$t >= treatment_time
      pre_idx <- !post_idx
      list(
        post_control = y[post_idx & y$inferred_process == "control", , drop = FALSE],
        post_treated = y[post_idx & y$inferred_process == "treated", , drop = FALSE],
        filt_control = y[pre_idx & y$inferred_process == "control", , drop = FALSE],
        filt_treated = y[pre_idx & y$inferred_process == "treated", , drop = FALSE]
      )
    })

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
        rho_max = etas_rho_max
      )
      biv_free_names <- biv_names[biv_fr_idx]
      biv_alpha_free_pos <- if (hard_subcritical) {
        which(biv_free_names %in% .etas_biv_alpha_names)
      } else {
        integer(0)
      }
      biv_to_optim <- function(natural_free) {
        optim_free <- natural_free
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
      biv_precomps <- lapply(labellings[keepers], function(y) {
        r <- if (use_pre_history_for_biv) y else y[y$t >= treatment_time, ]
        r <- r[order(r$t), ]
        nn <- nrow(r)
        W0 <- rep(1.0, nn); W1 <- rep(1.0, nn)
        aS0 <- area_control
        aS1 <- area_treated
        W0[inside.owin(r$x, r$y, treated_state_space)] <- 0
        W1[inside.owin(r$x, r$y, control_state_space)] <- 0
        if (!is.null(treated_background_zero_before)) {
          W1[r$t < treated_background_zero_before] <- 0
        }
        if (!is.null(background_rate_var) && background_rate_var %in% names(r)) {
          W_cov <- r[[background_rate_var]]
          if (length(W_cov) != nn) stop("background_rate_var length mismatch in SEM bivariate precomp.")
          W_cov <- as.numeric(W_cov)
          W_cov[!is.finite(W_cov)] <- 0
          min_pos <- suppressWarnings(min(W_cov[W_cov > 0], na.rm = TRUE))
          if (!is.finite(min_pos)) min_pos <- 1e-12
          W_cov[W_cov <= 0] <- min_pos
          W0 <- W0 * W_cov
          W1 <- W1 * W_cov
        }
        if (aS0 <= 0) aS0 <- 1; if (aS1 <= 0) aS1 <- 1
        list(realiz = r, precomp = list(W_0 = W0, W_1 = W1,
                                        areaS_0 = aS0, areaS_1 = aS1))
      })
      biv_ll_extra_names <- intersect(
        names(dots),
        c(
          "m0", "beta_gr", "enforce_finite_trigger_moments",
          "p_lower_bound", "q_lower_bound", "finite_moment_soft_width",
          "finite_moment_soft_weight", "finite_moment_soft_power",
          "enforce_alpha_subcritical", "alpha_beta_gap_min",
          "alpha_beta_soft_gap", "alpha_beta_soft_weight",
          "alpha_beta_soft_power", "max_branching_radius",
          "stability_barrier_start", "stability_barrier_weight",
          "stability_barrier_power"
        )
      )
      biv_ll_extra <- dots[biv_ll_extra_names]
      biv_obj <- function(par15) {
        liks <- sapply(biv_precomps, function(pc) {
          do.call(
            loglik_etas_bivariate,
            c(
              list(
                params = par15, realiz = pc$realiz,
                windowT = biv_wT, windowS = statespace,
                control_state_space = control_state_space,
                treated_state_space = treated_state_space,
                background_rate_var = background_rate_var,
                treated_background_zero_before = treated_background_zero_before,
                t_trunc = t_trunc, precomp = pc$precomp
              ),
              biv_ll_extra
            )
          )
        })
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
      if (length(biv_fr_idx) > 0) {
        biv_wrap <- function(optim_free) {
          p15 <- biv_par
          p15[biv_fr_idx] <- biv_from_optim(optim_free)
          biv_obj(p15)
        }
        biv_opt_start <- biv_to_optim(biv_par[biv_fr_idx])
        biv_res <- tryCatch(
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
        biv_par[biv_fr_idx] <- biv_from_optim(biv_res$par)
      } else {
        biv_res <- list(par = numeric(0), convergence = 0)
      }
      names(biv_par) <- biv_names
      if (hard_subcritical) {
        biv_par <- .etas_project_subcritical_biv(
          biv_par,
          beta_gr = etas_beta_eff,
          gap_min = etas_gap_min,
          rho_max = etas_rho_max
        )
      }
      biv_res$par <- biv_par
      biv_res$branching_radius <- .etas_biv_spectral_radius(
        biv_par, etas_beta_eff
      )
      dots$etas_bivariate_params <- biv_par

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
      filt_key <- if (process_label == "control") "filt_control" else "filt_treated"
      prepared_for_process <- lapply(prepared_labellings, function(parts) {
        post_part <- parts[[split_key]]
        if (!is_etas) {
          filt_part <- if (hawkes_use_filtration_history) {
            parts[[filt_key]]
          } else {
            post_part[0, c("x", "y", "t"), drop = FALSE]
          }
          if (nrow(post_part) < 1L) return(NULL)
          post_part <- post_part[order(post_part$t), , drop = FALSE]
          filt_part <- filt_part[filt_part$t < treatment_time, , drop = FALSE]
          filt_part <- filt_part[order(filt_part$t), , drop = FALSE]

          W_post <- if (!is.null(hawkes_bg_var) && hawkes_bg_var %in% names(post_part)) {
            as.numeric(post_part[[hawkes_bg_var]])
          } else {
            rep(1, nrow(post_part))
          }
          W_post[!is.finite(W_post)] <- 0
          total_area <- spatstat.geom::area(statespace)
          active_area <- total_area
          if (!is.null(zero_bg_region)) {
            zero_area <- spatstat.geom::area(zero_bg_region)
            active_area <- max(1e-12, total_area - zero_area)
            W_post[inside.owin(post_part$x, post_part$y, w = zero_bg_region)] <- 0
          }

          list(
            post_t = as.numeric(post_part$t),
            post_x = as.numeric(post_part$x),
            post_y = as.numeric(post_part$y),
            W_post = as.numeric(W_post),
            parent_t = as.numeric(c(filt_part$t, post_part$t)),
            parent_x = as.numeric(c(filt_part$x, post_part$x)),
            parent_y = as.numeric(c(filt_part$y, post_part$y)),
            active_area = active_area
          )
        } else {
          pc <- precompute_loglik_args(post_part, statespace, zero_bg_region)
          list(
            post_part = post_part,
            precomp = list(
              active_area = pc$active_area,
              in_zero_bg = pc$in_zero_bg_all
            )
          )
        }
      })
      etas_ll_extra_names <- if (is_etas) {
        intersect(
          names(dots),
          c(
            "m0", "beta_gr", "enforce_finite_trigger_moments",
            "p_lower_bound", "q_lower_bound", "finite_moment_soft_width",
            "finite_moment_soft_weight", "finite_moment_soft_power",
            "enforce_alpha_subcritical", "alpha_beta_gap_min",
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
      obj_fn <- function(params) {
        liks <- sapply(prepared_for_process, function(parts) {
          if (is.null(parts)) return(-Inf)
          if (!is_etas) {
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
            loglik <- hawkes_loglik_inhom_filtration_cpp(
              post_t = parts$post_t,
              post_x = parts$post_x,
              post_y = parts$post_y,
              W_val = parts$W_post,
              parent_t = parts$parent_t,
              parent_x = parts$parent_x,
              parent_y = parts$parent_y,
              mu = mu,
              alpha = alpha,
              beta = beta,
              K = K,
              areaS = parts$active_area,
              t_start = treatment_time,
              t_end = max_data_t,
              adjust_factor = 1.0,
              t_trunc = if (!is.null(t_trunc)) t_trunc else -1.0,
              kernel_type = hawkes_kernel_type(hawkes_kernel),
              cc = if (is.null(cc)) 1.0 else as.numeric(cc),
              p = if (is.null(p)) 2.0 else as.numeric(p),
              spatial_kernel_type = hawkes_spatial_kernel_type(hawkes_spatial_kernel),
              spatial_q = q_spatial,
              spatial_d = d_spatial
            )
            if (!is.finite(loglik)) return(-Inf)
            loglik
          } else {
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
                  t_trunc = t_trunc
                ),
                etas_ll_extra
              )
            )
          }
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
          eta_max = etas_eta_max
        )
      }
      free_names <- all_names[fr_idx]
      alpha_free_pos <- if (is_etas && hard_subcritical) {
        which(free_names == "alpha_m")
      } else {
        integer(0)
      }
      to_optim <- function(natural_free) {
        optim_free <- natural_free
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
      if (length(fr_idx) > 0) {
        opt_start <- to_optim(full_vec[fr_idx])
        wrap_fn <- function(optim_free) {
          p4 <- full_vec
          p4[fr_idx] <- from_optim(optim_free)
          obj_fn(p4)
        }
        res <- tryCatch(
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
        out_par <- full_vec
        out_par[fr_idx] <- from_optim(res$par)
      } else {
        res <- list(par = numeric(0), convergence = 0)
        out_par <- full_vec
      }
      if (is_etas && hard_subcritical) {
        out_par <- .etas_project_subcritical(
          out_par,
          beta_gr = etas_beta_eff,
          gap_min = etas_gap_min,
          eta_max = etas_eta_max
        )
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
