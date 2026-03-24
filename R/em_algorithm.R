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

  dots <- list(...)
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
    include_starting_data = FALSE, iter = 100, n_props = 10,
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
    result <- em_style_labelling(
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
      metric_name = "post_likelihood",
      optim_method = "max",
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
      model_type = model_type,
      ...
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
    p <- unlist(params)
    mu <- p[1]; alpha <- p[2]; beta <- p[3]; K <- p[4]
    if (!is.finite(mu) || !is.finite(alpha) || !is.finite(beta) || !is.finite(K)) return(-Inf)
    if (mu < 0 || alpha < 0 || beta <= 0 || K < 0 || K >= 1) return(-Inf)
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
      t_trunc = if (!is.null(t_trunc)) t_trunc else -1
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
          treated_background_zero_before = treated_background_zero_before,
          t_trunc = t_trunc
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
      labellings <- lapply(1:N_labellings, function(i) {
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
      if (verbose) cat(sprintf("[SEM] Labelling generation complete (%d + baseline, took %.1fs)\n",
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
      biv_obj <- function(par15) {
        liks <- sapply(biv_precomps, function(pc) {
          loglik_etas_bivariate(
            params = par15, realiz = pc$realiz,
            windowT = biv_wT, windowS = statespace,
            control_state_space = control_state_space,
            treated_state_space = treated_state_space,
            background_rate_var = background_rate_var,
            treated_background_zero_before = treated_background_zero_before,
            t_trunc = t_trunc, precomp = pc$precomp
          )
        })
        sum(liks * weights)
      }

      if (verbose) {
        cat("\n--- Outer optim: bivariate ETAS (joint) ---\n")
        cat(sprintf("  starting par: %s\n",
                    paste(biv_names, signif(biv_par, 5), sep = "=", collapse = "  ")))
      }

      t0 <- proc.time()[3]
      if (length(biv_fp_idx) > 0) {
        biv_wrap <- function(free_par) {
          p15 <- biv_par; p15[biv_fr_idx] <- free_par
          biv_obj(p15)
        }
        biv_res <- tryCatch(
          optim(par = biv_par[biv_fr_idx], fn = biv_wrap, method = "Nelder-Mead",
                control = list(fnscale = -1, trace = 0, maxit = outer_maxit_biv)),
          error = function(e) {
            cat(sprintf("  [bivariate] OPTIM ERROR: %s\n", e$message))
            list(par = biv_par[biv_fr_idx], convergence = -99)
          }
        )
        biv_par[biv_fr_idx] <- biv_res$par
      } else {
        biv_res <- tryCatch(
          optim(par = biv_par, fn = biv_obj, method = "Nelder-Mead",
                control = list(fnscale = -1, trace = 0, maxit = outer_maxit_biv)),
          error = function(e) {
            cat(sprintf("  [bivariate] OPTIM ERROR: %s\n", e$message))
            list(par = biv_par, convergence = -99)
          }
        )
        biv_par <- biv_res$par
      }
      names(biv_par) <- biv_names
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
      obj_fn <- function(params) {
        liks <- sapply(prepared_labellings, function(parts) {
          post_part <- parts[[split_key]]
          if (!is_etas) {
            filt_part <- if (hawkes_use_filtration_history) {
              parts[[filt_key]]
            } else {
              post_part[0, c("x", "y", "t"), drop = FALSE]
            }
            hawkes_loglik_with_filtration(
              params = params,
              post_realiz = as.data.frame(post_part),
              filt_realiz = as.data.frame(filt_part),
              zero_bg_region = zero_bg_region
            )
          } else {
            loglik_fn(
              params = params,
              realiz = post_part,
              windowT = sem_windowT,
              windowS = statespace, zero_background_region = zero_bg_region,
              background_rate_var = "W",
              t_trunc = t_trunc
            )
          }
        })
        return(sum(liks * weights))
      }

      full_vec <- unlist(par_list)
      if (verbose) {
        cat(sprintf("  [%s] starting par: %s\n", process_label,
                    paste(all_names, signif(full_vec, 5), sep = "=", collapse = "  ")))
        cat(sprintf("  [%s] starting lik: %s\n", process_label, signif(obj_fn(full_vec), 6)))
      }

      t0 <- proc.time()[3]
      if (length(fp_idx) > 0) {
        wrap_fn <- function(free_par) {
          p4 <- full_vec; p4[fr_idx] <- free_par
          obj_fn(p4)
        }
        res <- tryCatch(
          optim(par = full_vec[fr_idx], fn = wrap_fn, method = "Nelder-Mead",
                control = list(fnscale = -1, trace = 0, maxit = outer_maxit)),
          error = function(e) { cat(sprintf("  [%s] OPTIM ERROR: %s\n", process_label, e$message)); list(par = full_vec[fr_idx], convergence = -99) }
        )
        out_par <- full_vec; out_par[fr_idx] <- res$par
        res$par <- out_par
      } else {
        res <- tryCatch(
          optim(par = full_vec, fn = obj_fn, method = "Nelder-Mead",
                control = list(fnscale = -1, trace = 0, maxit = outer_maxit)),
          error = function(e) { cat(sprintf("  [%s] OPTIM ERROR: %s\n", process_label, e$message)); list(par = full_vec, convergence = -99) }
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
