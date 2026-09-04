# ============================================================================
# Bivariate ETAS with cross-excitation
#
# Two-component ETAS where events in each process can trigger offspring in
# both processes via a 2x2 kernel matrix:
#
#   g_{kl}: excitation from process l (parent) to process k (child)
#     g_00 = control self-excitation
#     g_01 = treated -> control cross-excitation
#     g_10 = control -> treated cross-excitation
#     g_11 = treated self-excitation
#
# Each kernel component has its own (A, alpha_m) productivity pair.
# Structural parameters (c, p, D, gamma, q) are shared across all kernels.
# Setting A_01 = A_10 = 0 recovers two independent ETAS processes.
# ============================================================================

.etas_bivariate_par_names <- c(
  "mu_0", "mu_1",
  "A_00", "alpha_m_00",
  "A_11", "alpha_m_11",
  "A_01", "alpha_m_01",
  "A_10", "alpha_m_10",
  "c", "p", "D", "gamma", "q"
)

.etas_biv_alpha_names <- c("alpha_m_00", "alpha_m_11", "alpha_m_01", "alpha_m_10")
.etas_biv_A_names <- c("A_00", "A_11", "A_01", "A_10")

#' Gutenberg-Richter branching spectral radius for bivariate ETAS
#'
#' Channel mean offspring is \code{eta_ij = A_ij * beta / (beta - alpha_ij)}
#' (Inf if \code{alpha_ij >= beta}). Returns the spectral radius of the
#' 2x2 matrix with entries \code{eta_00, eta_01; eta_10, eta_11}.
#' @keywords internal
.etas_biv_spectral_radius <- function(params, beta_gr) {
  pv <- if (is.list(params)) unlist(params) else params
  beta_eff <- .etas_resolve_beta_gr(beta_gr)
  if (!is.finite(beta_eff)) return(NA_real_)
  eta_ij <- function(A, alpha) {
    A <- suppressWarnings(as.numeric(A))
    alpha <- suppressWarnings(as.numeric(alpha))
    if (!is.finite(A) || !is.finite(alpha)) return(NA_real_)
    gap <- beta_eff - alpha
    if (!is.finite(gap) || gap <= 1e-8) return(Inf)
    A * beta_eff / gap
  }
  a <- eta_ij(pv[["A_00"]], pv[["alpha_m_00"]])
  b <- eta_ij(pv[["A_01"]], pv[["alpha_m_01"]])
  c <- eta_ij(pv[["A_10"]], pv[["alpha_m_10"]])
  d <- eta_ij(pv[["A_11"]], pv[["alpha_m_11"]])
  if (!all(is.finite(c(a, b, c, d)))) return(Inf)
  # Closed-form Perron root of the nonnegative 2x2 branching matrix.
  disc <- (a - d) * (a - d) + 4 * b * c
  if (disc < 0) disc <- 0
  0.5 * (a + d + sqrt(disc))
}

#' Whether a bivariate ETAS parameter vector is complete, finite, and
#' strictly inside the configured spectral-radius margin.
#' @keywords internal
.etas_biv_params_ok <- function(params, beta_gr, rho_max = 0.98,
                                required_names = .etas_bivariate_par_names) {
  if (is.null(params)) return(FALSE)
  pv <- if (is.list(params)) unlist(params) else params
  if (is.null(names(pv)) || !all(required_names %in% names(pv))) return(FALSE)
  vals <- suppressWarnings(as.numeric(pv[required_names]))
  if (length(vals) != length(required_names) || !all(is.finite(vals))) {
    return(FALSE)
  }
  rho_max <- suppressWarnings(as.numeric(rho_max)[1L])
  if (!is.finite(rho_max) || rho_max <= 0) rho_max <- 0.98
  rho <- .etas_biv_spectral_radius(pv, beta_gr)
  is.finite(rho) && rho < rho_max
}

#' Whether a naive \code{fit_etas_bivariate} result is usable.
#'
#' Requires optimizer \code{convergence == 0}, a finite objective, a complete
#' finite parameter vector, and spectral radius strictly below \code{rho_max}.
#' Finite parameters after a failed Nelder-Mead run are not a successful fit.
#' @keywords internal
.etas_biv_fit_ok <- function(fit, beta_gr, rho_max = 0.98,
                             required_names = .etas_bivariate_par_names) {
  if (is.null(fit) || is.null(fit$par)) return(FALSE)
  if (is.null(fit$convergence) || as.integer(fit$convergence)[1L] != 0L) {
    return(FALSE)
  }
  val <- suppressWarnings(as.numeric(fit$value)[1L])
  if (!is.finite(val)) return(FALSE)
  .etas_biv_params_ok(fit$par, beta_gr, rho_max, required_names)
}

#' Whether a bivariate SEM result is usable.
#'
#' Requires a complete finite subcritical parameter vector and last inner /
#' outer Nelder-Mead \code{convergence == 0}. Finite parameters after a failed
#' M-step are not a successful SEM fit. Missing convergence fails closed.
#' @keywords internal
.etas_biv_sem_ok <- function(sem, beta_gr, rho_max = 0.98,
                             required_names = .etas_bivariate_par_names) {
  if (is.null(sem)) return(FALSE)
  if (!.etas_biv_params_ok(
    sem$etas_bivariate_params, beta_gr, rho_max, required_names
  )) {
    return(FALSE)
  }
  conv <- suppressWarnings(as.integer(sem$etas_bivariate_convergence)[1L])
  if (length(conv) != 1L || is.na(conv) || conv != 0L) return(FALSE)
  val <- suppressWarnings(as.numeric(sem$etas_bivariate_value)[1L])
  if (length(val) == 1L && !is.na(val) && !is.finite(val)) return(FALSE)
  TRUE
}

#' Project bivariate ETAS params into a hard-subcritical set
#'
#' Clamps each \code{alpha_m_*} into
#' \code{(alpha_m_lower_bound, beta - gap_min)} and, if needed,
#' scales all \code{A_*} so the GR spectral radius is at most \code{rho_max}.
#' @keywords internal
.etas_project_subcritical_biv <- function(params,
                                          beta_gr,
                                          gap_min = 1e-4,
                                          rho_max = 0.98,
                                          alpha_m_lower_bound = 0) {
  pv <- if (is.list(params)) unlist(params) else params
  if (is.null(names(pv))) names(pv) <- .etas_bivariate_par_names
  bounds <- .etas_alpha_bounds(
    beta_gr, gap_min = gap_min, alpha_m_lower_bound = alpha_m_lower_bound
  )
  beta_eff <- bounds$beta
  if (!is.finite(beta_eff)) return(pv)
  rho_max <- suppressWarnings(as.numeric(rho_max))
  if (!is.finite(rho_max) || rho_max <= 0) rho_max <- 0.98
  for (nm in .etas_biv_alpha_names) {
    if (is.finite(pv[[nm]])) {
      a <- pv[[nm]]
      if (is.finite(bounds$lo)) a <- max(a, bounds$lo + 1e-8)
      if (is.finite(bounds$hi)) a <- min(a, bounds$hi - 1e-4)
      pv[[nm]] <- a
    }
  }
  rho <- .etas_biv_spectral_radius(pv, beta_eff)
  if (is.finite(rho) && rho > rho_max && rho > 0) {
    scale <- (rho_max / rho) * (1 - 1e-6)
    for (nm in .etas_biv_A_names) {
      if (is.finite(pv[[nm]]) && pv[[nm]] > 0) pv[[nm]] <- pv[[nm]] * scale
    }
  }
  pv
}

.etas_scale_to_target_rho <- function(params, beta_gr, target_rho) {
  pv <- if (is.list(params)) unlist(params) else params
  if (is.null(names(pv))) names(pv) <- .etas_bivariate_par_names
  target_rho <- suppressWarnings(as.numeric(target_rho))
  if (length(target_rho) != 1L || !is.finite(target_rho) || target_rho <= 0) {
    return(pv)
  }
  rho <- .etas_biv_spectral_radius(pv, beta_gr)
  if (!is.finite(rho) || rho <= 0) return(pv)
  scale <- target_rho / rho
  for (nm in .etas_biv_A_names) {
    if (is.finite(pv[[nm]]) && pv[[nm]] > 0) pv[[nm]] <- pv[[nm]] * scale
  }
  pv
}

.etas_raw_loglik_biv <- function(params, ll_args) {
  args <- ll_args
  args$params <- params
  args$max_branching_radius <- Inf
  args$stability_barrier_weight <- 0
  tryCatch(
    as.numeric(do.call(loglik_etas_bivariate, args)),
    error = function(e) -1e15
  )
}

.etas_polish_productivity_biv <- function(params,
                                          ll_args,
                                          beta_gr,
                                          rho_max = 0.98,
                                          rho_lo = 0.05) {
  pv <- if (is.list(params)) unlist(params) else params
  if (is.null(names(pv))) names(pv) <- .etas_bivariate_par_names
  rho_max <- suppressWarnings(as.numeric(rho_max))
  if (!is.finite(rho_max) || rho_max <= 0) rho_max <- 0.98
  rho_lo <- suppressWarnings(as.numeric(rho_lo))
  if (!is.finite(rho_lo) || rho_lo <= 0) rho_lo <- 0.05
  hi <- rho_max * (1 - 1e-6)
  if (!(hi > rho_lo)) {
    return(list(par = pv, value = .etas_raw_loglik_biv(pv, ll_args), polished = FALSE))
  }
  obj <- function(rho) {
    q <- .etas_scale_to_target_rho(pv, beta_gr, rho)
    .etas_raw_loglik_biv(q, ll_args)
  }
  opt <- tryCatch(
    stats::optimize(obj, interval = c(rho_lo, hi), maximum = TRUE),
    error = function(e) NULL
  )
  if (is.null(opt) || !is.finite(opt$objective)) {
    return(list(par = pv, value = .etas_raw_loglik_biv(pv, ll_args), polished = FALSE))
  }
  best <- .etas_scale_to_target_rho(pv, beta_gr, opt$maximum)
  list(par = best, value = as.numeric(opt$objective), polished = TRUE)
}

#' Background time exposure under an activation cutoff
#'
#' Length of the part of \code{windowT} on which a background component is
#' switched on when it only activates at \code{cutoff} (e.g. the treated
#' background under \code{treated_background_zero_before}). \code{NULL} or
#' non-finite \code{cutoff} means the background is on for the whole window.
#' @keywords internal
.etas_bg_exposure <- function(windowT, cutoff = NULL) {
  tval <- windowT[2] - windowT[1]
  if (is.null(cutoff)) return(tval)
  cut <- suppressWarnings(as.numeric(cutoff))
  if (length(cut) != 1L || !is.finite(cut)) return(tval)
  max(0, windowT[2] - max(windowT[1], cut))
}

#' Parse a background policy cutoff time
#'
#' Returns a single numeric cutoff (possibly infinite) or \code{NULL} when the
#' input is \code{NULL}/\code{NA}/malformed, meaning "no policy cutoff".
#' @keywords internal
.etas_bg_cutoff <- function(x) {
  if (is.null(x)) return(NULL)
  v <- suppressWarnings(as.numeric(x))
  if (length(v) != 1L || is.na(v)) return(NULL)
  v
}

#' Effective control-background exposure when the pre-cutoff support is larger
#'
#' Under \code{control_background_everywhere_before = cutoff} the control
#' background covers the whole domain before \code{cutoff} and only its own
#' region afterwards. With the spatial density held fixed (density-continuous
#' convention), each pre-cutoff day charges \code{mu_0 * mass_ratio} instead of
#' \code{mu_0}, where \code{mass_ratio} is the full-domain to control-region
#' background mass ratio (area ratio for a flat background, KDE mass ratio for
#' a weighted one). The returned value is the effective time exposure
#' \code{mass_ratio * pre_len + post_len} to multiply \code{mu_0} by in the
#' compensator.
#' @keywords internal
.etas_bg_exposure_control <- function(windowT, cutoff = NULL, mass_ratio = NULL) {
  tval <- windowT[2] - windowT[1]
  cut <- .etas_bg_cutoff(cutoff)
  if (is.null(cut)) return(tval)
  mr <- suppressWarnings(as.numeric(mass_ratio))
  if (length(mr) != 1L || !is.finite(mr) || is.na(mr) || mr <= 0) mr <- 1
  cut <- min(max(cut, windowT[1]), windowT[2])
  mr * (cut - windowT[1]) + (windowT[2] - cut)
}

#' Log-likelihood for bivariate ETAS with cross-excitation
#'
#' @param params Named list or vector with the 15 bivariate ETAS parameters.
#' @param realiz Data frame with columns x, y, t, mag, and inferred_process
#'   (character "control"/"treated" or integer 0/1).
#' @param windowT Numeric vector c(start, end).
#' @param windowS An owin object for the full observation window.
#' @param m0 Reference magnitude. NULL = min(mag).
#' @param control_state_space owin for the control region.
#' @param treated_state_space owin for the treated region.
#' @param background_rate_var Optional column of \code{realiz} with a
#'   mean-one spatial background weight (e.g. KDE), multiplied into both
#'   background weights. Skipped when a \code{precomp} is supplied
#'   (precomputed weights must already include it).
#' @param treated_background_zero_before Optional activation time of the
#'   treated background. Event-side weights before this time are zeroed
#'   (baked into \code{precomp} weights when one is supplied) and the
#'   \code{mu_1} compensator is charged only from this time onward, so this
#'   argument must be passed even alongside a \code{precomp}.
#' @param control_background_everywhere_before Optional policy cutoff time
#'   (same scale as \code{realiz$t}/\code{windowT}) before which the control
#'   background covers the whole domain. Events before this time keep their
#'   control background weight even inside \code{treated_state_space}; from
#'   this time onward the control background is zero there (the pre-fix
#'   behaviour of zeroing at all times corresponds to \code{NULL}). The
#'   \code{mu_0} compensator is charged \code{mu_0 *
#'   control_background_pre_mass_ratio} per pre-cutoff unit time
#'   (density-continuous convention). Must be passed even alongside a
#'   \code{precomp} (which must already bake the event-side mask).
#' @param control_background_pre_mass_ratio Ratio of full-domain to
#'   control-region background mass used for the pre-cutoff \code{mu_0}
#'   compensator charge. Defaults to \code{(areaS_0 + areaS_1) / areaS_0},
#'   which is exact for a flat background; callers using
#'   \code{background_rate_var} (e.g. a KDE) must pass the KDE mass ratio
#'   \code{integral(full) / integral(control region)}.
#' @param t_trunc Temporal truncation (NULL = none).
#' @param precomp Optional list with \code{W_0}, \code{W_1}, \code{areaS_0},
#'   \code{areaS_1}, and optionally \code{process_id}. If \code{t_shifted}
#'   is supplied and matches \code{nrow(realiz)}, geometry checks (sort,
#'   time-window subset, \code{inside.owin}) are skipped and \code{t_shifted}
#'   is passed to the C++ kernel. Optional \code{x}, \code{y}, \code{mag},
#'   \code{t_expo_0}, \code{t_expo_1}, \code{t_max} hoist further per-eval work.
#' @param ... Ignored (kept so older callers can still pass \code{precomp}
#'   through \code{...}).
#' @return Scalar log-likelihood.
#' @export
loglik_etas_bivariate <- function(params,
                                  realiz,
                                  windowT,
                                  windowS,
                                  m0 = NULL,
                                  control_state_space = NULL,
                                  treated_state_space = NULL,
                                  background_rate_var = NULL,
                                  treated_background_zero_before = NULL,
                                  control_background_everywhere_before = NULL,
                                  control_background_pre_mass_ratio = NULL,
                                  beta_gr = NULL,
                                  enforce_finite_trigger_moments = TRUE,
                                  p_lower_bound = 2.001,
                                  q_lower_bound = 1.501,
                                  finite_moment_soft_width = 0.05,
                                  finite_moment_soft_weight = 2000,
                                  finite_moment_soft_power = 2,
                                  enforce_alpha_subcritical = TRUE,
                                  alpha_beta_gap_min = 1e-4,
                                  # Default: larger magnitudes are more productive (alpha_m > 0).
                                  # Set to -Inf to disable the lower bound.
                                  alpha_m_lower_bound = 0,
                                  alpha_beta_soft_gap = 0.05,
                                  alpha_beta_soft_weight = 2000,
                                  alpha_beta_soft_power = 2,
                                  # Hard upper bound on GR branching spectral radius.
                                  # NULL/Inf disables the hard reject.
                                  max_branching_radius = Inf,
                                  # Soft stability barrier is off by default
                                  # (weight 0); the hard cap is the guardrail.
                                  stability_barrier_start = 0.95,
                                  stability_barrier_weight = 0,
                                  stability_barrier_power = 4,
                                  t_trunc = NULL,
                                  precomp = NULL,
                                  ...) {
  if (is.list(params) && !is.null(names(params))) {
    pv <- unlist(params)
  } else {
    pv <- as.numeric(params)
    if (is.null(names(pv))) names(pv) <- .etas_bivariate_par_names
  }

  mu_0       <- pv[["mu_0"]]
  mu_1       <- pv[["mu_1"]]
  A_00       <- pv[["A_00"]]
  alpha_m_00 <- pv[["alpha_m_00"]]
  A_11       <- pv[["A_11"]]
  alpha_m_11 <- pv[["alpha_m_11"]]
  A_01       <- pv[["A_01"]]
  alpha_m_01 <- pv[["alpha_m_01"]]
  A_10       <- pv[["A_10"]]
  alpha_m_10 <- pv[["alpha_m_10"]]
  cc         <- pv[["c"]]
  p          <- pv[["p"]]
  D          <- pv[["D"]]
  gamma_p    <- pv[["gamma"]]
  q          <- pv[["q"]]

  core_par <- c(mu_0, mu_1, A_00, A_11, A_01, A_10, cc, p, D, gamma_p, q)
  if (any(!is.finite(core_par))) return(-1e15)
  if (min(mu_0, mu_1, A_00, A_11, cc, D) < 0 ||
      A_01 < 0 || A_10 < 0 || q <= 1 || gamma_p < 0) {
    return(-1e15)
  }
  if (!.etas_omori_p_admissible(p, t_trunc)) return(-1e15)
  p_min <- suppressWarnings(as.numeric(p_lower_bound))
  if (length(p_min) != 1L || !is.finite(p_min) || is.na(p_min)) p_min <- 2.001
  q_min <- suppressWarnings(as.numeric(q_lower_bound))
  if (length(q_min) != 1L || !is.finite(q_min) || is.na(q_min)) q_min <- 1.501
  if (isTRUE(enforce_finite_trigger_moments) && (p <= p_min || q <= q_min)) return(-1e15)

  if (is.null(precomp)) {
    dots_pre <- list(...)
    if (!is.null(dots_pre$precomp)) precomp <- dots_pre$precomp
  }
  t_shifted <- if (!is.null(precomp)) precomp$t_shifted else NULL
  n_realiz <- nrow(realiz)
  frozen_geom <- !is.null(t_shifted) && length(t_shifted) == n_realiz && n_realiz > 0L

  if (!frozen_geom) {
    if (!all(is.finite(realiz$t)) || !all(is.finite(realiz$x)) ||
        !all(is.finite(realiz$y)) || !all(is.finite(realiz$mag))) {
      return(-1e15)
    }
    if (is.unsorted(realiz$t)) realiz <- realiz[order(realiz$t), ]
    t_idx <- realiz$t >= windowT[1] & realiz$t <= windowT[2]
    if (!all(t_idx)) realiz <- realiz[t_idx, ]
    if (!all(spatstat.geom::inside.owin(
      realiz$x, realiz$y, spatstat.geom::as.owin(windowS)
    ))) {
      return(-1e15)
    }
    n <- nrow(realiz)
    if (n == 0) return(-1e15)
    tt <- realiz$t - windowT[1]
    xx <- realiz$x
    yy <- realiz$y
    mm <- realiz$mag
  } else {
    n <- n_realiz
    tt <- t_shifted
    xx <- if (!is.null(precomp$x) && length(precomp$x) == n) precomp$x else realiz$x
    yy <- if (!is.null(precomp$y) && length(precomp$y) == n) precomp$y else realiz$y
    mm <- if (!is.null(precomp$mag) && length(precomp$mag) == n) precomp$mag else realiz$mag
  }

  if (is.null(m0)) m0 <- min(realiz$mag)
  beta_eff <- .etas_resolve_beta_gr(beta_gr, realiz = realiz, m0 = m0)
  alpha_vals <- c(alpha_m_00, alpha_m_11, alpha_m_01, alpha_m_10)
  alpha_lo <- suppressWarnings(as.numeric(alpha_m_lower_bound))
  if (length(alpha_lo) != 1L || is.na(alpha_lo)) alpha_lo <- 0
  if (is.finite(alpha_lo) && any(!is.finite(alpha_vals) | alpha_vals <= alpha_lo)) {
    return(-1e15)
  }
  if (isTRUE(enforce_alpha_subcritical)) {
    gap_min <- suppressWarnings(as.numeric(alpha_beta_gap_min))
    if (length(gap_min) != 1L || !is.finite(gap_min) || is.na(gap_min) || gap_min < 0) gap_min <- 1e-4
    if (!is.finite(beta_eff) || is.na(beta_eff) || any(alpha_vals >= (beta_eff - gap_min))) return(-1e15)
  }
  rho_max <- suppressWarnings(as.numeric(max_branching_radius))
  if (length(rho_max) == 1L && is.finite(rho_max) && rho_max > 0 && is.finite(beta_eff)) {
    rho_hat <- .etas_biv_spectral_radius(pv, beta_eff)
    if (!is.finite(rho_hat) || rho_hat >= rho_max) return(-1e15)
  }

  # Build process_id vector (0 = control, 1 = treated)
  if (!is.null(precomp) && !is.null(precomp$process_id) && length(precomp$process_id) == n) {
    process_id <- as.integer(precomp$process_id)
  } else if ("inferred_process" %in% names(realiz)) {
    proc_col <- realiz$inferred_process
    if (is.character(proc_col)) {
      process_id <- as.integer(proc_col == "treated")
    } else {
      process_id <- as.integer(proc_col)
    }
  } else if ("process" %in% names(realiz)) {
    proc_col <- realiz$process
    if (is.character(proc_col)) {
      process_id <- as.integer(proc_col == "treated")
    } else {
      process_id <- as.integer(proc_col)
    }
  } else if ("location_process" %in% names(realiz)) {
    proc_col <- realiz$location_process
    if (is.character(proc_col)) {
      process_id <- as.integer(proc_col == "treated")
    } else {
      process_id <- as.integer(proc_col)
    }
  } else {
    stop("realiz must have an inferred_process, process, or location_process column")
  }

  # Background weights: W_0 = 0 in treated region, W_1 = 0 in control region
  if (!is.null(precomp)) {
    W_0     <- precomp$W_0
    W_1     <- precomp$W_1
    areaS_0 <- precomp$areaS_0
    areaS_1 <- precomp$areaS_1
  } else {
    W_0 <- rep(1.0, n)
    W_1 <- rep(1.0, n)

    total_area <- spatstat.geom::area(as.owin(windowS))

    if (!is.null(treated_state_space)) {
      areaS_1 <- spatstat.geom::area(as.owin(treated_state_space))
      in_treated <- inside.owin(realiz$x, realiz$y, treated_state_space)
      cb_cut <- .etas_bg_cutoff(control_background_everywhere_before)
      if (is.null(cb_cut)) {
        W_0[in_treated] <- 0
      } else {
        # Pre-cutoff the control background covers the whole domain; it is
        # only excluded from the treated region once treatment switches on.
        W_0[in_treated & realiz$t >= cb_cut] <- 0
      }
    } else {
      areaS_1 <- total_area / 2
    }

    if (!is.null(control_state_space)) {
      areaS_0 <- spatstat.geom::area(as.owin(control_state_space))
      in_control <- inside.owin(realiz$x, realiz$y, control_state_space)
      W_1[in_control] <- 0
    } else {
      areaS_0 <- total_area - areaS_1
    }

    if (!is.finite(areaS_0) || areaS_0 <= 0) areaS_0 <- 1
    if (!is.finite(areaS_1) || areaS_1 <= 0) areaS_1 <- 1
  }
  if (!is.finite(areaS_0) || areaS_0 <= 0) areaS_0 <- 1
  if (!is.finite(areaS_1) || areaS_1 <= 0) areaS_1 <- 1

  # Optional inhomogeneous background covariate / treated-background time mask.
  # When precomp is supplied these are assumed already baked into W_0/W_1
  # (SEM outer/inner precomp does that); re-applying would double-weight.
  if (is.null(precomp)) {
    if (!is.null(background_rate_var) && background_rate_var %in% names(realiz)) {
      W_cov <- realiz[[background_rate_var]]
      if (length(W_cov) != n) stop("background_rate_var length mismatch in realiz.")
      W_cov <- as.numeric(W_cov)
      W_cov[!is.finite(W_cov)] <- 0
      min_pos <- suppressWarnings(min(W_cov[W_cov > 0], na.rm = TRUE))
      if (!is.finite(min_pos)) min_pos <- 1e-12
      W_cov[W_cov <= 0] <- min_pos
      W_0 <- W_0 * W_cov
      W_1 <- W_1 * W_cov
    }
    if (!is.null(treated_background_zero_before)) {
      W_1[realiz$t < as.numeric(treated_background_zero_before)] <- 0
    }
  }

  tval <- if (!is.null(precomp) && !is.null(precomp$t_max) &&
              is.finite(precomp$t_max)) {
    precomp$t_max
  } else {
    windowT[2] - windowT[1]
  }
  # Compensator exposure: the treated background is only on from the cutoff
  # onward, so mu_1 must be charged over that span only (never over the full
  # window). Callers may hoist the two exposures onto precomp when the
  # window and policy cutoffs are frozen (SEM M-step / MLE).
  if (!is.null(precomp) && !is.null(precomp$t_expo_0) && !is.null(precomp$t_expo_1) &&
      is.finite(precomp$t_expo_0) && is.finite(precomp$t_expo_1)) {
    t_expo_0 <- precomp$t_expo_0
    t_expo_1 <- precomp$t_expo_1
  } else {
    t_expo_1 <- .etas_bg_exposure(windowT, treated_background_zero_before)
    cb_mass <- suppressWarnings(as.numeric(control_background_pre_mass_ratio))
    if (length(cb_mass) != 1L || !is.finite(cb_mass) || is.na(cb_mass) || cb_mass <= 0) {
      cb_mass <- (areaS_0 + areaS_1) / areaS_0
    }
    t_expo_0 <- .etas_bg_exposure_control(
      windowT, control_background_everywhere_before, cb_mass
    )
  }

  loglik <- etas_bivariate_loglik_cpp(
    t         = tt,
    x         = xx,
    y         = yy,
    mag       = mm,
    process_id = as.integer(process_id),
    W_val_0   = W_0,
    W_val_1   = W_1,
    mu_0 = mu_0, mu_1 = mu_1,
    A_00 = A_00, alpha_m_00 = alpha_m_00,
    A_11 = A_11, alpha_m_11 = alpha_m_11,
    A_01 = A_01, alpha_m_01 = alpha_m_01,
    A_10 = A_10, alpha_m_10 = alpha_m_10,
    cc = cc, p = p, D = D, gamma_par = gamma_p, q = q,
    m0 = m0,
    areaS_0 = areaS_0, areaS_1 = areaS_1,
    t_max = tval,
    t_trunc = if (!is.null(t_trunc)) t_trunc else -1.0,
    bg_exposure_0 = t_expo_0,
    bg_exposure_1 = t_expo_1
  )
  if (!is.finite(loglik)) return(-1e15)

  # Hybrid boundary handling: hard admissibility plus smooth interior penalties.
  if (isTRUE(enforce_finite_trigger_moments)) {
    soft_width <- suppressWarnings(as.numeric(finite_moment_soft_width))
    soft_weight <- suppressWarnings(as.numeric(finite_moment_soft_weight))
    soft_power <- suppressWarnings(as.numeric(finite_moment_soft_power))
    if (length(soft_width) != 1L || !is.finite(soft_width) || is.na(soft_width) || soft_width <= 0) soft_width <- 0
    if (length(soft_weight) != 1L || !is.finite(soft_weight) || is.na(soft_weight) || soft_weight <= 0) soft_weight <- 0
    if (length(soft_power) != 1L || !is.finite(soft_power) || is.na(soft_power) || soft_power <= 0) soft_power <- 2
    if (soft_width > 0 && soft_weight > 0) {
      p_excess <- max(0, (p_min + soft_width) - p)
      q_excess <- max(0, (q_min + soft_width) - q)
      if (p_excess > 0) loglik <- loglik - soft_weight * (p_excess ^ soft_power)
      if (q_excess > 0) loglik <- loglik - soft_weight * (q_excess ^ soft_power)
    }
  }
  if (isTRUE(enforce_alpha_subcritical)) {
    beta_soft <- suppressWarnings(as.numeric(beta_gr))
    if (length(beta_soft) != 1L || !is.finite(beta_soft) || is.na(beta_soft) || beta_soft <= 0) {
      mag_delta <- as.numeric(realiz$mag) - as.numeric(m0)
      mag_delta <- mag_delta[is.finite(mag_delta) & mag_delta > 0]
      beta_soft <- if (length(mag_delta) > 0L) 1 / mean(mag_delta) else NA_real_
    }
    gap_min <- suppressWarnings(as.numeric(alpha_beta_gap_min))
    if (length(gap_min) != 1L || !is.finite(gap_min) || is.na(gap_min) || gap_min < 0) gap_min <- 1e-4
    soft_gap <- suppressWarnings(as.numeric(alpha_beta_soft_gap))
    if (length(soft_gap) != 1L || !is.finite(soft_gap) || is.na(soft_gap) || soft_gap <= 0) soft_gap <- 0
    soft_weight <- suppressWarnings(as.numeric(alpha_beta_soft_weight))
    if (length(soft_weight) != 1L || !is.finite(soft_weight) || is.na(soft_weight) || soft_weight <= 0) soft_weight <- 0
    soft_power <- suppressWarnings(as.numeric(alpha_beta_soft_power))
    if (length(soft_power) != 1L || !is.finite(soft_power) || is.na(soft_power) || soft_power <= 0) soft_power <- 2
    if (is.finite(beta_soft) && soft_gap > 0 && soft_weight > 0) {
      alpha_gaps <- beta_soft - c(alpha_m_00, alpha_m_11, alpha_m_01, alpha_m_10)
      min_gap <- min(alpha_gaps, na.rm = TRUE)
      gap_excess <- max(0, (gap_min + soft_gap) - min_gap)
      if (gap_excess > 0) loglik <- loglik - soft_weight * (gap_excess ^ soft_power)
    }
  }
  # Soft stability barrier (supplementary to hard max_branching_radius reject).
  barrier_weight <- suppressWarnings(as.numeric(stability_barrier_weight))
  if (length(barrier_weight) != 1L || !is.finite(barrier_weight) || is.na(barrier_weight) || barrier_weight <= 0) barrier_weight <- 0
  if (barrier_weight > 0 && is.finite(beta_eff)) {
    rho <- .etas_biv_spectral_radius(pv, beta_eff)
    if (!is.finite(rho)) return(-1e15)
    barrier_start <- suppressWarnings(as.numeric(stability_barrier_start))
    if (length(barrier_start) != 1L || !is.finite(barrier_start) || is.na(barrier_start)) barrier_start <- 0.95
    barrier_power <- suppressWarnings(as.numeric(stability_barrier_power))
    if (length(barrier_power) != 1L || !is.finite(barrier_power) || is.na(barrier_power) || barrier_power <= 0) barrier_power <- 2
    excess <- max(0, rho - barrier_start)
    if (excess > 0) {
      loglik <- loglik - barrier_weight * (excess ^ barrier_power)
    }
  }
  loglik
}

#' Shared soft-penalty term for bivariate ETAS (params only; labelling-invariant)
#' @keywords internal
.etas_bivariate_param_penalty <- function(pv,
                                          beta_eff,
                                          enforce_finite_trigger_moments = TRUE,
                                          p_lower_bound = 2.001,
                                          q_lower_bound = 1.501,
                                          finite_moment_soft_width = 0.05,
                                          finite_moment_soft_weight = 2000,
                                          finite_moment_soft_power = 2,
                                          enforce_alpha_subcritical = TRUE,
                                          alpha_beta_gap_min = 1e-4,
                                          alpha_beta_soft_gap = 0.05,
                                          alpha_beta_soft_weight = 2000,
                                          alpha_beta_soft_power = 2,
                                          stability_barrier_start = 0.95,
                                          stability_barrier_weight = 0,
                                          stability_barrier_power = 4) {
  pen <- 0
  p <- pv[["p"]]; q <- pv[["q"]]
  if (isTRUE(enforce_finite_trigger_moments)) {
    p_min <- suppressWarnings(as.numeric(p_lower_bound))
    if (length(p_min) != 1L || !is.finite(p_min) || is.na(p_min)) p_min <- 2.001
    q_min <- suppressWarnings(as.numeric(q_lower_bound))
    if (length(q_min) != 1L || !is.finite(q_min) || is.na(q_min)) q_min <- 1.501
    soft_width <- suppressWarnings(as.numeric(finite_moment_soft_width))
    soft_weight <- suppressWarnings(as.numeric(finite_moment_soft_weight))
    soft_power <- suppressWarnings(as.numeric(finite_moment_soft_power))
    if (length(soft_width) != 1L || !is.finite(soft_width) || soft_width <= 0) soft_width <- 0
    if (length(soft_weight) != 1L || !is.finite(soft_weight) || soft_weight <= 0) soft_weight <- 0
    if (length(soft_power) != 1L || !is.finite(soft_power) || soft_power <= 0) soft_power <- 2
    if (soft_width > 0 && soft_weight > 0) {
      p_excess <- max(0, (p_min + soft_width) - p)
      q_excess <- max(0, (q_min + soft_width) - q)
      if (p_excess > 0) pen <- pen + soft_weight * (p_excess ^ soft_power)
      if (q_excess > 0) pen <- pen + soft_weight * (q_excess ^ soft_power)
    }
  }
  if (isTRUE(enforce_alpha_subcritical) && is.finite(beta_eff)) {
    gap_min <- suppressWarnings(as.numeric(alpha_beta_gap_min))
    if (length(gap_min) != 1L || !is.finite(gap_min) || gap_min < 0) gap_min <- 1e-4
    soft_gap <- suppressWarnings(as.numeric(alpha_beta_soft_gap))
    if (length(soft_gap) != 1L || !is.finite(soft_gap) || soft_gap <= 0) soft_gap <- 0
    soft_weight <- suppressWarnings(as.numeric(alpha_beta_soft_weight))
    if (length(soft_weight) != 1L || !is.finite(soft_weight) || soft_weight <= 0) soft_weight <- 0
    soft_power <- suppressWarnings(as.numeric(alpha_beta_soft_power))
    if (length(soft_power) != 1L || !is.finite(soft_power) || soft_power <= 0) soft_power <- 2
    if (soft_gap > 0 && soft_weight > 0) {
      alpha_gaps <- beta_eff - pv[.etas_biv_alpha_names]
      min_gap <- min(alpha_gaps, na.rm = TRUE)
      gap_excess <- max(0, (gap_min + soft_gap) - min_gap)
      if (gap_excess > 0) pen <- pen + soft_weight * (gap_excess ^ soft_power)
    }
  }
  barrier_weight <- suppressWarnings(as.numeric(stability_barrier_weight))
  if (length(barrier_weight) == 1L && is.finite(barrier_weight) && barrier_weight > 0 &&
      is.finite(beta_eff)) {
    rho <- .etas_biv_spectral_radius(pv, beta_eff)
    if (!is.finite(rho)) return(Inf)
    barrier_start <- suppressWarnings(as.numeric(stability_barrier_start))
    if (length(barrier_start) != 1L || !is.finite(barrier_start)) barrier_start <- 0.95
    barrier_power <- suppressWarnings(as.numeric(stability_barrier_power))
    if (length(barrier_power) != 1L || !is.finite(barrier_power) || barrier_power <= 0) {
      barrier_power <- 2
    }
    excess <- max(0, rho - barrier_start)
    if (excess > 0) pen <- pen + barrier_weight * (excess ^ barrier_power)
  }
  pen
}

#' Batched bivariate ETAS log-likelihoods for shared-geometry labellings
#'
#' Evaluates the joint log-likelihood under many labellings that share the
#' same event geometry \code{(t,x,y,mag)} and differ only in process labels
#' (and optionally background masks). Pairwise kernel terms are computed once.
#' Soft penalties that depend only on parameters are applied once to every
#' entry (matching \code{loglik_etas_bivariate}).
#'
#' @param params Named length-15 parameter vector/list.
#' @param t,x,y,mag Shared event geometry (t already shifted so window starts at 0,
#'   or pass \code{windowT} and unshifted times).
#' @param process_ids Integer matrix \code{n x K} (0/1).
#' @param W0s,W1s Numeric matrices \code{n x K} of background weights, or
#'   \code{n x 1} to recycle the same mask across labellings.
#' @param areaS_0,areaS_1 Active areas.
#' @param t_max Window length. If NULL, uses \code{diff(windowT)}.
#' @param windowT Optional \code{c(start,end)}; used when \code{t_max} is NULL,
#'   to shift \code{t} when \code{t_already_shifted=FALSE}, and to anchor
#'   \code{treated_background_zero_before}.
#' @param treated_background_zero_before Optional activation time (raw scale)
#'   of the treated background; requires \code{windowT}. Callers must still
#'   zero the pre-cutoff entries of \code{W1s} themselves; this argument only
#'   fixes the \code{mu_1} compensator exposure.
#' @param control_background_everywhere_before Optional policy cutoff (raw
#'   scale) before which the control background covers the whole domain;
#'   requires \code{windowT}. Callers must themselves keep the pre-cutoff
#'   treated-region entries of \code{W0s} unmasked; this argument only fixes
#'   the \code{mu_0} compensator exposure.
#' @param control_background_pre_mass_ratio Full-domain to control-region
#'   background mass ratio for the pre-cutoff \code{mu_0} charge. Defaults to
#'   \code{(areaS_0 + areaS_1) / areaS_0} (flat background); pass the KDE
#'   mass ratio when \code{W0s} carries KDE weights.
#' @param t_already_shifted If FALSE, subtracts \code{windowT[1]} from \code{t}.
#' @param n_threads Worker threads for the C++ kernel (default 1).
#' @param ... Same barrier/penalty arguments as \code{loglik_etas_bivariate}.
#' @return Length-K numeric vector of log-likelihoods.
#' @export
loglik_etas_bivariate_batch <- function(params,
                                        t, x, y, mag,
                                        process_ids,
                                        W0s, W1s,
                                        areaS_0, areaS_1,
                                        t_max = NULL,
                                        windowT = NULL,
                                        treated_background_zero_before = NULL,
                                        control_background_everywhere_before = NULL,
                                        control_background_pre_mass_ratio = NULL,
                                        m0 = NULL,
                                        beta_gr = NULL,
                                        enforce_finite_trigger_moments = TRUE,
                                        p_lower_bound = 2.001,
                                        q_lower_bound = 1.501,
                                        finite_moment_soft_width = 0.05,
                                        finite_moment_soft_weight = 2000,
                                        finite_moment_soft_power = 2,
                                        enforce_alpha_subcritical = TRUE,
                                        alpha_beta_gap_min = 1e-4,
                                        alpha_m_lower_bound = 0,
                                        alpha_beta_soft_gap = 0.05,
                                        alpha_beta_soft_weight = 2000,
                                        alpha_beta_soft_power = 2,
                                        max_branching_radius = Inf,
                                        stability_barrier_start = 0.95,
                                        stability_barrier_weight = 0,
                                        stability_barrier_power = 4,
                                        t_trunc = NULL,
                                        t_already_shifted = TRUE,
                                        n_threads = 1L,
                                        ...) {
  if (is.list(params) && !is.null(names(params))) {
    pv <- unlist(params)
  } else {
    pv <- as.numeric(params)
    if (is.null(names(pv))) names(pv) <- .etas_bivariate_par_names
  }
  K <- ncol(process_ids)
  if (is.null(K) || K < 1L) return(numeric(0))
  reject <- rep(-1e15, K)

  mu_0 <- pv[["mu_0"]]; mu_1 <- pv[["mu_1"]]
  A_00 <- pv[["A_00"]]; alpha_m_00 <- pv[["alpha_m_00"]]
  A_11 <- pv[["A_11"]]; alpha_m_11 <- pv[["alpha_m_11"]]
  A_01 <- pv[["A_01"]]; alpha_m_01 <- pv[["alpha_m_01"]]
  A_10 <- pv[["A_10"]]; alpha_m_10 <- pv[["alpha_m_10"]]
  cc <- pv[["c"]]; p <- pv[["p"]]; D <- pv[["D"]]
  gamma_p <- pv[["gamma"]]; q <- pv[["q"]]

  core_par <- c(mu_0, mu_1, A_00, A_11, A_01, A_10, cc, p, D, gamma_p, q)
  if (any(!is.finite(core_par))) return(reject)
  if (min(mu_0, mu_1, A_00, A_11, cc, D) < 0 ||
      A_01 < 0 || A_10 < 0 || q <= 1 || gamma_p < 0) {
    return(reject)
  }
  if (!.etas_omori_p_admissible(p, t_trunc)) return(reject)
  p_min <- suppressWarnings(as.numeric(p_lower_bound))
  if (length(p_min) != 1L || !is.finite(p_min) || is.na(p_min)) p_min <- 2.001
  q_min <- suppressWarnings(as.numeric(q_lower_bound))
  if (length(q_min) != 1L || !is.finite(q_min) || is.na(q_min)) q_min <- 1.501
  if (isTRUE(enforce_finite_trigger_moments) && (p <= p_min || q <= q_min)) return(reject)

  n <- length(t)
  if (n == 0L || length(x) != n || length(y) != n || length(mag) != n) return(reject)
  if (nrow(process_ids) != n || nrow(W0s) != n || nrow(W1s) != n) {
    stop("process_ids/W0s/W1s must have n = length(t) rows.")
  }
  w0_ok <- ncol(W0s) == K || ncol(W0s) == 1L
  w1_ok <- ncol(W1s) == K || ncol(W1s) == 1L
  if (!isTRUE(w0_ok) || !isTRUE(w1_ok)) {
    stop("W0s/W1s must be n x K, or n x 1 to recycle shared masks across labellings.")
  }
  if (!all(is.finite(t)) || !all(is.finite(x)) || !all(is.finite(y)) || !all(is.finite(mag))) {
    return(reject)
  }

  if (is.null(m0)) m0 <- min(mag)
  beta_eff <- .etas_resolve_beta_gr(beta_gr, m0 = m0)
  alpha_vals <- c(alpha_m_00, alpha_m_11, alpha_m_01, alpha_m_10)
  alpha_lo <- suppressWarnings(as.numeric(alpha_m_lower_bound))
  if (length(alpha_lo) != 1L || is.na(alpha_lo)) alpha_lo <- 0
  if (is.finite(alpha_lo) && any(!is.finite(alpha_vals) | alpha_vals <= alpha_lo)) {
    return(reject)
  }
  if (isTRUE(enforce_alpha_subcritical)) {
    gap_min <- suppressWarnings(as.numeric(alpha_beta_gap_min))
    if (length(gap_min) != 1L || !is.finite(gap_min) || gap_min < 0) gap_min <- 1e-4
    if (!is.finite(beta_eff) || any(alpha_vals >= (beta_eff - gap_min))) return(reject)
  }
  rho_max <- suppressWarnings(as.numeric(max_branching_radius))
  if (length(rho_max) == 1L && is.finite(rho_max) && rho_max > 0 && is.finite(beta_eff)) {
    rho_hat <- .etas_biv_spectral_radius(pv, beta_eff)
    if (!is.finite(rho_hat) || rho_hat >= rho_max) return(reject)
  }

  if (is.null(t_max)) {
    if (is.null(windowT)) stop("Provide t_max or windowT.")
    t_max <- windowT[2] - windowT[1]
  }
  tt <- if (isTRUE(t_already_shifted)) {
    as.numeric(t)
  } else {
    if (is.null(windowT)) stop("t_already_shifted=FALSE requires windowT.")
    as.numeric(t) - windowT[1]
  }
  # mu_1 compensator exposure under the treated-background activation cutoff.
  bg_expo_1 <- t_max
  if (!is.null(treated_background_zero_before)) {
    if (is.null(windowT)) {
      stop("treated_background_zero_before requires windowT in loglik_etas_bivariate_batch.")
    }
    bg_expo_1 <- .etas_bg_exposure(windowT, treated_background_zero_before)
  }
  # mu_0 compensator exposure when the control background covers the whole
  # domain before its policy cutoff (density-continuous convention).
  bg_expo_0 <- t_max
  if (!is.null(control_background_everywhere_before)) {
    if (is.null(windowT)) {
      stop("control_background_everywhere_before requires windowT in loglik_etas_bivariate_batch.")
    }
    cb_mass <- suppressWarnings(as.numeric(control_background_pre_mass_ratio))
    if (length(cb_mass) != 1L || !is.finite(cb_mass) || is.na(cb_mass) || cb_mass <= 0) {
      cb_mass <- (areaS_0 + areaS_1) / areaS_0
    }
    bg_expo_0 <- .etas_bg_exposure_control(
      windowT, control_background_everywhere_before, cb_mass
    )
  }

  liks <- etas_bivariate_loglik_batch_cpp(
    t = tt, x = as.numeric(x), y = as.numeric(y), mag = as.numeric(mag),
    process_ids = as.matrix(process_ids),
    W0s = as.matrix(W0s), W1s = as.matrix(W1s),
    mu_0 = mu_0, mu_1 = mu_1,
    A_00 = A_00, alpha_m_00 = alpha_m_00,
    A_11 = A_11, alpha_m_11 = alpha_m_11,
    A_01 = A_01, alpha_m_01 = alpha_m_01,
    A_10 = A_10, alpha_m_10 = alpha_m_10,
    cc = cc, p = p, D = D, gamma_par = gamma_p, q = q,
    m0 = m0, areaS_0 = areaS_0, areaS_1 = areaS_1,
    t_max = t_max,
    t_trunc = if (!is.null(t_trunc)) as.numeric(t_trunc) else -1.0,
    bg_exposure_0 = bg_expo_0,
    bg_exposure_1 = bg_expo_1,
    n_threads = as.integer(n_threads)
  )
  liks[!is.finite(liks)] <- -1e15
  pen <- .etas_bivariate_param_penalty(
    pv, beta_eff,
    enforce_finite_trigger_moments = enforce_finite_trigger_moments,
    p_lower_bound = p_lower_bound, q_lower_bound = q_lower_bound,
    finite_moment_soft_width = finite_moment_soft_width,
    finite_moment_soft_weight = finite_moment_soft_weight,
    finite_moment_soft_power = finite_moment_soft_power,
    enforce_alpha_subcritical = enforce_alpha_subcritical,
    alpha_beta_gap_min = alpha_beta_gap_min,
    alpha_beta_soft_gap = alpha_beta_soft_gap,
    alpha_beta_soft_weight = alpha_beta_soft_weight,
    alpha_beta_soft_power = alpha_beta_soft_power,
    stability_barrier_start = stability_barrier_start,
    stability_barrier_weight = stability_barrier_weight,
    stability_barrier_power = stability_barrier_power
  )
  if (!is.finite(pen)) return(reject)
  liks - pen
}


#' Sequential MAP process labels under bivariate ETAS
#'
#' Time-ordered assignment: non-assignable events (typically pre-treatment)
#' keep \code{process_id_init}; each assignable event is labelled
#' \code{argmax_k \lambda_k} given already-assigned parents. Matches the
#' intensity used by \code{loglik_etas_bivariate}.
#'
#' @param params Named bivariate ETAS parameter vector.
#' @param t,x,y,mag Event geometry (\code{t} on the same scale as
#'   \code{t_trunc}; not shifted unless the caller shifted it).
#' @param assignable Integer 0/1; 1 means this row may be relabelled.
#' @param process_id_init Integer 0/1 labels used for non-assignable rows.
#' @param W0,W1 Background weight vectors.
#' @param areaS_0,areaS_1 Active areas.
#' @param m0 Reference magnitude.
#' @param t_trunc Temporal truncation; \code{NULL} disables.
#' @return Integer vector of process ids (0 = control, 1 = treated).
#' @keywords internal
.sequential_label_etas_bivariate <- function(params, t, x, y, mag,
                                             assignable, process_id_init,
                                             W0, W1, areaS_0, areaS_1,
                                             m0, t_trunc = NULL,
                                             sample_bernoulli = FALSE) {
  pv <- if (is.list(params) && !is.null(names(params))) unlist(params) else {
    v <- as.numeric(params)
    if (is.null(names(v))) names(v) <- .etas_bivariate_par_names
    v
  }
  tt <- suppressWarnings(as.numeric(t_trunc))
  if (length(tt) != 1L || !is.finite(tt) || tt <= 0) tt <- -1.0
  etas_bivariate_sequential_map_cpp(
    t = as.numeric(t), x = as.numeric(x), y = as.numeric(y), mag = as.numeric(mag),
    assignable = as.integer(assignable),
    process_id_init = as.integer(process_id_init),
    W_val_0 = as.numeric(W0), W_val_1 = as.numeric(W1),
    mu_0 = pv[["mu_0"]], mu_1 = pv[["mu_1"]],
    A_00 = pv[["A_00"]], alpha_m_00 = pv[["alpha_m_00"]],
    A_11 = pv[["A_11"]], alpha_m_11 = pv[["alpha_m_11"]],
    A_01 = pv[["A_01"]], alpha_m_01 = pv[["alpha_m_01"]],
    A_10 = pv[["A_10"]], alpha_m_10 = pv[["alpha_m_10"]],
    cc = pv[["c"]], p = pv[["p"]], D = pv[["D"]],
    gamma_par = if (is.finite(pv[["gamma"]])) pv[["gamma"]] else 0,
    q = pv[["q"]],
    m0 = m0, areaS_0 = areaS_0, areaS_1 = areaS_1, t_trunc = tt,
    sample_bernoulli = isTRUE(sample_bernoulli)
  )
}

sequential_map_etas_bivariate <- function(params, t, x, y, mag,
                                          assignable, process_id_init,
                                          W0, W1, areaS_0, areaS_1,
                                          m0, t_trunc = NULL) {
  .sequential_label_etas_bivariate(
    params, t, x, y, mag, assignable, process_id_init,
    W0, W1, areaS_0, areaS_1, m0, t_trunc = t_trunc,
    sample_bernoulli = FALSE
  )
}

#' Sequential Bernoulli process labels for bivariate ETAS
#'
#' Same walk as \code{sequential_map_etas_bivariate}, but each assignable
#' event is drawn \code{Z_i ~ Bern(\lambda_1 / (\lambda_0 + \lambda_1))}
#' using already-sampled parents.
#' @keywords internal
sequential_bernoulli_etas_bivariate <- function(params, t, x, y, mag,
                                               assignable, process_id_init,
                                               W0, W1, areaS_0, areaS_1,
                                               m0, t_trunc = NULL) {
  .sequential_label_etas_bivariate(
    params, t, x, y, mag, assignable, process_id_init,
    W0, W1, areaS_0, areaS_1, m0, t_trunc = t_trunc,
    sample_bernoulli = TRUE
  )
}

#' Fit bivariate ETAS via MLE
#'
#' @param params_init Initial parameter values (named list or vector of 15).
#' @param realiz Data frame with x, y, t, mag, and process labels.
#' @param windowT Numeric c(start, end).
#' @param windowS owin observation window.
#' @param m0 Reference magnitude (NULL = auto).
#' @param control_state_space owin for control region.
#' @param treated_state_space owin for treated region.
#' @param treated_background_zero_before Optional treated-background
#'   activation time (see \code{loglik_etas_bivariate}).
#' @param control_background_everywhere_before Optional policy cutoff before
#'   which the control background covers the whole domain (see
#'   \code{loglik_etas_bivariate}).
#' @param control_background_pre_mass_ratio Full-domain to control-region
#'   background mass ratio for the pre-cutoff \code{mu_0} charge; pass the
#'   KDE mass ratio when using \code{background_rate_var}.
#' @param maxit Maximum optim iterations.
#' @param fixed_params Named list of parameters to hold fixed.
#' @param symmetric If TRUE, constrain A_01=A_10 and alpha_m_01=alpha_m_10.
#' @param trace Trace level for optim.
#' @param t_trunc Temporal truncation.
#' @param hard_subcritical If TRUE (default), reparameterize free
#'   \code{alpha_m_*} into
#'   \code{(alpha_m_lower_bound, beta - gap)} (default lower bound 0),
#'   project the start into \code{rho < max_branching_radius}, and project
#'   the returned fit.
#' @param ... Passed to loglik_etas_bivariate (include \code{beta_gr}).
#' @return An optim result with par as a named length-15 vector.
#' @export
fit_etas_bivariate <- function(params_init,
                               realiz,
                               windowT,
                               windowS,
                               m0 = NULL,
                               control_state_space = NULL,
                               treated_state_space = NULL,
                               background_rate_var = NULL,
                               treated_background_zero_before = NULL,
                               control_background_everywhere_before = NULL,
                               control_background_pre_mass_ratio = NULL,
                               maxit = 5000,
                               fixed_params = NULL,
                               symmetric = FALSE,
                               trace = 0,
                               t_trunc = NULL,
                               hard_subcritical = TRUE,
                               log_transform = TRUE,
                               soft_branching_barrier = TRUE,
                               polish_productivity = TRUE,
                               interior_restart = TRUE,
                               ...) {
  dots <- list(...)
  all_names <- .etas_bivariate_par_names

  if (inherits(params_init, "list")) params_init <- unlist(params_init)
  if (!is.null(names(params_init))) {
    full_init <- params_init[all_names]
  } else {
    full_init <- params_init
    names(full_init) <- all_names
  }
  if (!is.null(fixed_params)) {
    for (nm in names(fixed_params)) full_init[nm] <- fixed_params[[nm]]
  }

  if (is.null(m0)) m0 <- min(realiz$mag)

  p_min <- suppressWarnings(as.numeric(dots$p_lower_bound))
  if (length(p_min) != 1L || !is.finite(p_min) || is.na(p_min)) p_min <- 2.001
  q_min <- suppressWarnings(as.numeric(dots$q_lower_bound))
  if (length(q_min) != 1L || !is.finite(q_min) || is.na(q_min)) q_min <- 1.501
  if (!isFALSE(dots$enforce_finite_trigger_moments)) {
    if (is.finite(full_init["p"])) full_init["p"] <- max(full_init["p"], p_min + 1e-3)
    if (is.finite(full_init["q"])) full_init["q"] <- max(full_init["q"], q_min + 1e-3)
  }
  gap_min <- suppressWarnings(as.numeric(dots$alpha_beta_gap_min))
  if (length(gap_min) != 1L || !is.finite(gap_min) || is.na(gap_min) || gap_min < 0) gap_min <- 1e-4
  rho_max <- suppressWarnings(as.numeric(dots$max_branching_radius))
  if (length(rho_max) != 1L || !is.finite(rho_max) || rho_max <= 0) rho_max <- 0.98
  alpha_lo <- suppressWarnings(as.numeric(dots$alpha_m_lower_bound))
  if (length(alpha_lo) != 1L || is.na(alpha_lo)) alpha_lo <- 0
  init_margin <- suppressWarnings(as.numeric(dots$init_branching_margin))
  if (length(init_margin) != 1L || !is.finite(init_margin) ||
      init_margin <= 0 || init_margin > 1) {
    init_margin <- 0.9
  }
  near_cap_frac <- suppressWarnings(as.numeric(dots$near_cap_frac))
  if (length(near_cap_frac) != 1L || !is.finite(near_cap_frac) ||
      near_cap_frac <= 0 || near_cap_frac > 1) {
    near_cap_frac <- 0.99
  }
  restart_rho <- suppressWarnings(as.numeric(dots$interior_restart_rho))
  if (length(restart_rho) != 1L || !is.finite(restart_rho) || restart_rho <= 0) {
    restart_rho <- min(0.70, rho_max * init_margin)
  }
  dots$alpha_beta_gap_min <- gap_min
  dots$alpha_m_lower_bound <- alpha_lo
  beta_eff <- .etas_resolve_beta_gr(dots$beta_gr, realiz = realiz, m0 = m0)
  alpha_bounds <- .etas_alpha_bounds(
    beta_eff, gap_min = gap_min, alpha_m_lower_bound = alpha_lo
  )
  init_rho_cap <- rho_max * init_margin
  if (isTRUE(hard_subcritical) || !isFALSE(dots$enforce_alpha_subcritical) ||
      is.finite(alpha_lo)) {
    if (is.finite(beta_eff)) {
      full_init <- .etas_project_subcritical_biv(
        full_init, beta_eff, gap_min = gap_min, rho_max = init_rho_cap,
        alpha_m_lower_bound = alpha_lo
      )
    }
  }

  use_soft_barrier <- isTRUE(soft_branching_barrier) && is.finite(rho_max)
  barrier_ctrl <- .etas_soft_barrier_controls(
    rho_max,
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

  # Handle symmetric constraint: fix A_10 = A_01 and alpha_m_10 = alpha_m_01
  sym_fixed <- NULL
  if (symmetric) {
    sym_fixed <- c("A_10", "alpha_m_10")
    full_init["A_10"] <- full_init["A_01"]
    full_init["alpha_m_10"] <- full_init["alpha_m_01"]
  }

  all_fixed_names <- unique(c(names(fixed_params), sym_fixed))
  fixed_idx <- if (length(all_fixed_names) > 0) {
    match(all_fixed_names, all_names)
  } else {
    integer(0)
  }
  free_idx <- setdiff(seq_along(all_names), fixed_idx)
  free_names <- all_names[free_idx]
  A_free <- intersect(.etas_biv_A_names, free_names)
  # Hard alpha constraint via reparameterization onto (alpha_lo, alpha_hi):
  #   alpha = alpha_lo + (alpha_hi - alpha_lo) * plogis(u)
  # with alpha_hi = beta - gap_min (default alpha_lo = 0).
  use_alpha_reparam <- (isTRUE(hard_subcritical) || is.finite(alpha_lo)) &&
    is.finite(beta_eff)
  alpha_free_pos <- if (use_alpha_reparam) {
    which(free_names %in% .etas_biv_alpha_names)
  } else {
    integer(0)
  }
  log_names <- c("mu_0", "mu_1", .etas_biv_A_names, "c", "D")
  log_free_pos <- if (isTRUE(log_transform)) {
    which(free_names %in% log_names)
  } else {
    integer(0)
  }
  free_to_natural <- function(free_par) {
    nat <- free_par
    if (length(log_free_pos)) {
      nat[log_free_pos] <- exp(free_par[log_free_pos])
    }
    if (length(alpha_free_pos)) {
      nat[alpha_free_pos] <- .etas_opt_to_alpha(
        free_par[alpha_free_pos], alpha_bounds$lo, alpha_bounds$hi
      )
    }
    nat
  }
  natural_to_free <- function(nat_par) {
    free_par <- nat_par
    if (length(log_free_pos)) {
      pos <- pmax(as.numeric(nat_par[log_free_pos]), 1e-12)
      free_par[log_free_pos] <- log(pos)
    }
    if (length(alpha_free_pos)) {
      free_par[alpha_free_pos] <- .etas_alpha_to_opt(
        nat_par[alpha_free_pos], alpha_bounds$lo, alpha_bounds$hi
      )
    }
    free_par
  }

  finite_rows <- is.finite(realiz$t) & is.finite(realiz$x) & is.finite(realiz$y) & is.finite(realiz$mag)
  if (!all(finite_rows)) realiz <- realiz[finite_rows, , drop = FALSE]
  realiz <- realiz[order(realiz$t), , drop = FALSE]
  t_idx <- is.finite(realiz$t) & realiz$t >= windowT[1] & realiz$t <= windowT[2]
  if (!all(t_idx)) realiz <- realiz[t_idx, , drop = FALSE]
  n <- nrow(realiz)

  W_0 <- rep(1.0, n); W_1 <- rep(1.0, n)
  total_area <- spatstat.geom::area(as.owin(windowS))
  if (!is.null(treated_state_space)) {
    areaS_1 <- spatstat.geom::area(as.owin(treated_state_space))
    in_treated_fit <- inside.owin(realiz$x, realiz$y, treated_state_space)
    cb_cut_fit <- .etas_bg_cutoff(control_background_everywhere_before)
    if (is.null(cb_cut_fit)) {
      W_0[in_treated_fit] <- 0
    } else {
      # Control background covers the whole domain before the policy cutoff.
      W_0[in_treated_fit & realiz$t >= cb_cut_fit] <- 0
    }
  } else { areaS_1 <- total_area / 2 }
  if (!is.null(control_state_space)) {
    areaS_0 <- spatstat.geom::area(as.owin(control_state_space))
    W_1[inside.owin(realiz$x, realiz$y, control_state_space)] <- 0
  } else { areaS_0 <- total_area - areaS_1 }
  if (!is.finite(areaS_0) || areaS_0 <= 0) areaS_0 <- max(1, total_area / 2)
  if (!is.finite(areaS_1) || areaS_1 <= 0) areaS_1 <- max(1, total_area / 2)

  if (!is.null(background_rate_var) && background_rate_var %in% names(realiz)) {
    W_cov <- realiz[[background_rate_var]]
    if (length(W_cov) != n) stop("background_rate_var length mismatch in realiz.")
    W_cov <- as.numeric(W_cov)
    W_cov[!is.finite(W_cov)] <- 0
    min_pos <- suppressWarnings(min(W_cov[W_cov > 0], na.rm = TRUE))
    if (!is.finite(min_pos)) min_pos <- 1e-12
    W_cov[W_cov <= 0] <- min_pos
    W_0 <- W_0 * W_cov
    W_1 <- W_1 * W_cov
  }
  # Policy mask must be baked in here: loglik_etas_bivariate skips event-side
  # masking when precomp is supplied (it still uses the cutoff, forwarded via
  # ll_args, for the mu_1 compensator exposure).
  if (!is.null(treated_background_zero_before)) {
    W_1[realiz$t < as.numeric(treated_background_zero_before)] <- 0
  }
  cb_mass_fit <- suppressWarnings(as.numeric(control_background_pre_mass_ratio))
  if (length(cb_mass_fit) != 1L || !is.finite(cb_mass_fit) ||
      is.na(cb_mass_fit) || cb_mass_fit <= 0) {
    cb_mass_fit <- (areaS_0 + areaS_1) / areaS_0
  }
  precomp <- list(
    W_0 = W_0, W_1 = W_1, areaS_0 = areaS_0, areaS_1 = areaS_1,
    t_expo_0 = .etas_bg_exposure_control(
      windowT, control_background_everywhere_before, cb_mass_fit
    ),
    t_expo_1 = .etas_bg_exposure(windowT, treated_background_zero_before),
    t_max = windowT[2] - windowT[1]
  )
  # Freeze geometry for Nelder-Mead only when every point is already inside
  # windowS; otherwise loglik must keep the inside.owin reject.
  if (n > 0L && all(spatstat.geom::inside.owin(
    realiz$x, realiz$y, spatstat.geom::as.owin(windowS)
  ))) {
    precomp$t_shifted <- realiz$t - windowT[1]
    precomp$x <- realiz$x
    precomp$y <- realiz$y
    precomp$mag <- realiz$mag
    proc_col <- if ("inferred_process" %in% names(realiz)) {
      realiz$inferred_process
    } else if ("process" %in% names(realiz)) {
      realiz$process
    } else if ("location_process" %in% names(realiz)) {
      realiz$location_process
    } else {
      NULL
    }
    if (!is.null(proc_col)) {
      precomp$process_id <- if (is.character(proc_col)) {
        as.integer(proc_col == "treated")
      } else {
        as.integer(proc_col)
      }
    }
  }

  build_ll_args <- function(par15, for_optim = TRUE) {
    ll_args <- list(
      params = par15, realiz = realiz,
      windowT = windowT, windowS = windowS,
      m0 = m0,
      control_state_space = control_state_space,
      treated_state_space = treated_state_space,
      background_rate_var = background_rate_var,
      treated_background_zero_before = treated_background_zero_before,
      control_background_everywhere_before = control_background_everywhere_before,
      control_background_pre_mass_ratio = control_background_pre_mass_ratio,
      t_trunc = t_trunc,
      precomp = precomp
    )
    if (!is.null(dots$beta_gr) || is.finite(beta_eff)) {
      ll_args$beta_gr <- if (!is.null(dots$beta_gr)) dots$beta_gr else beta_eff
    }
    ll_args$alpha_beta_gap_min <- gap_min
    ll_args$alpha_m_lower_bound <- alpha_lo
    if (!is.null(dots$enforce_finite_trigger_moments)) ll_args$enforce_finite_trigger_moments <- dots$enforce_finite_trigger_moments
    if (!is.null(dots$p_lower_bound)) ll_args$p_lower_bound <- dots$p_lower_bound
    if (!is.null(dots$q_lower_bound)) ll_args$q_lower_bound <- dots$q_lower_bound
    if (!is.null(dots$enforce_alpha_subcritical)) ll_args$enforce_alpha_subcritical <- dots$enforce_alpha_subcritical
    if (!is.null(dots$finite_moment_soft_width)) ll_args$finite_moment_soft_width <- dots$finite_moment_soft_width
    if (!is.null(dots$finite_moment_soft_weight)) ll_args$finite_moment_soft_weight <- dots$finite_moment_soft_weight
    if (!is.null(dots$finite_moment_soft_power)) ll_args$finite_moment_soft_power <- dots$finite_moment_soft_power
    if (!is.null(dots$alpha_beta_soft_gap)) ll_args$alpha_beta_soft_gap <- dots$alpha_beta_soft_gap
    if (!is.null(dots$alpha_beta_soft_weight)) ll_args$alpha_beta_soft_weight <- dots$alpha_beta_soft_weight
    if (!is.null(dots$alpha_beta_soft_power)) ll_args$alpha_beta_soft_power <- dots$alpha_beta_soft_power
    if (isTRUE(for_optim) && use_soft_barrier) {
      ll_args$max_branching_radius <- Inf
      ll_args$stability_barrier_start <- barrier_ctrl$stability_barrier_start
      ll_args$stability_barrier_weight <- barrier_ctrl$stability_barrier_weight
      ll_args$stability_barrier_power <- barrier_ctrl$stability_barrier_power
    } else if (isTRUE(for_optim)) {
      ll_args$max_branching_radius <- rho_max
      if (!is.null(dots$stability_barrier_start)) ll_args$stability_barrier_start <- dots$stability_barrier_start
      if (!is.null(dots$stability_barrier_weight)) ll_args$stability_barrier_weight <- dots$stability_barrier_weight
      if (!is.null(dots$stability_barrier_power)) ll_args$stability_barrier_power <- dots$stability_barrier_power
    } else {
      ll_args$max_branching_radius <- Inf
      ll_args$stability_barrier_weight <- 0
    }
    ll_args
  }

  profile_fn <- function(free_par) {
    par15 <- full_init
    par15[free_idx] <- free_to_natural(free_par)

    if (symmetric) {
      par15["A_10"] <- par15["A_01"]
      par15["alpha_m_10"] <- par15["alpha_m_01"]
    }

    ll_val <- tryCatch(
      do.call(loglik_etas_bivariate, build_ll_args(par15, for_optim = TRUE)),
      error = function(e) -1e15
    )
    if (!is.finite(ll_val) || is.na(ll_val)) return(-1e15)
    as.numeric(ll_val)
  }

  run_nm <- function(start_natural) {
    fit_local <- stats::optim(
      par     = natural_to_free(start_natural[free_idx]),
      fn      = profile_fn,
      method  = "Nelder-Mead",
      control = list(fnscale = -1, trace = trace, maxit = maxit)
    )
    par_local <- full_init
    par_local[free_idx] <- free_to_natural(fit_local$par)
    if (symmetric) {
      par_local["A_10"] <- par_local["A_01"]
      par_local["alpha_m_10"] <- par_local["alpha_m_01"]
    }
    if ((isTRUE(hard_subcritical) || is.finite(alpha_lo)) && is.finite(beta_eff)) {
      par_local <- .etas_project_subcritical_biv(
        par_local, beta_eff, gap_min = gap_min, rho_max = rho_max,
        alpha_m_lower_bound = alpha_lo
      )
    }
    list(fit = fit_local, par = par_local)
  }

  finalize_candidate <- function(par15, fit_meta = NULL) {
    if ((isTRUE(hard_subcritical) || is.finite(alpha_lo)) && is.finite(beta_eff)) {
      par15 <- .etas_project_subcritical_biv(
        par15, beta_eff, gap_min = gap_min, rho_max = rho_max,
        alpha_m_lower_bound = alpha_lo
      )
    }
    if (symmetric) {
      par15["A_10"] <- par15["A_01"]
      par15["alpha_m_10"] <- par15["alpha_m_01"]
    }
    raw_ll <- .etas_raw_loglik_biv(par15, build_ll_args(par15, for_optim = FALSE))
    rho <- if (is.finite(beta_eff)) .etas_biv_spectral_radius(par15, beta_eff) else NA_real_
    list(par = par15, value = raw_ll, branching_radius = rho, fit = fit_meta)
  }

  primary <- run_nm(full_init)
  best <- finalize_candidate(primary$par, primary$fit)
  near_cap <- is.finite(best$branching_radius) &&
    best$branching_radius >= (rho_max * near_cap_frac)
  distrust_conv <- !is.null(best$fit$convergence) &&
    (best$fit$convergence != 0L || near_cap)

  if (isTRUE(interior_restart) && length(A_free) > 0L && is.finite(beta_eff) &&
      (near_cap || distrust_conv ||
       !is.finite(best$value) || best$value <= -1e10)) {
    restart_start <- .etas_scale_to_target_rho(best$par, beta_eff, restart_rho)
    restart_start <- .etas_project_subcritical_biv(
      restart_start, beta_eff, gap_min = gap_min, rho_max = init_rho_cap,
      alpha_m_lower_bound = alpha_lo
    )
    restarted <- run_nm(restart_start)
    cand <- finalize_candidate(restarted$par, restarted$fit)
    if (is.finite(cand$value) &&
        (!is.finite(best$value) || cand$value > best$value + 1e-8)) {
      best <- cand
      best$interior_restarted <- TRUE
    }
  }

  if (isTRUE(polish_productivity) && length(A_free) > 0L && is.finite(beta_eff)) {
    polished <- .etas_polish_productivity_biv(
      best$par, build_ll_args(best$par, for_optim = FALSE),
      beta_eff, rho_max = rho_max
    )
    if (isTRUE(polished$polished) && is.finite(polished$value) &&
        (!is.finite(best$value) || polished$value > best$value + 1e-8)) {
      best <- finalize_candidate(polished$par, best$fit)
      best$polished_productivity <- TRUE
    }
  }

  fit <- if (!is.null(best$fit)) best$fit else list(convergence = 0L, counts = NA, message = NULL)
  fit$par <- best$par
  fit$value <- best$value
  fit$m0 <- m0
  fit$hard_subcritical <- isTRUE(hard_subcritical)
  fit$alpha_m_lower_bound <- alpha_lo
  fit$branching_radius <- best$branching_radius
  fit$soft_branching_barrier <- use_soft_barrier
  fit$log_transform <- isTRUE(log_transform)
  fit$polished_productivity <- isTRUE(best$polished_productivity)
  fit$interior_restarted <- isTRUE(best$interior_restarted)
  fit$near_cap <- is.finite(best$branching_radius) &&
    best$branching_radius >= (rho_max * near_cap_frac)
  return(fit)
}

.etas_biv_background_lmax <- function(mu, windowT, state_space, lookup,
                                       ref_area = NULL) {
  if (!is.function(lookup) || is.null(state_space) || !is.finite(mu) ||
      mu < 1e-10) {
    return(NULL)
  }
  area_ss <- spatstat.geom::area(state_space)
  ref_area <- suppressWarnings(as.numeric(ref_area))
  if (length(ref_area) != 1L || !is.finite(ref_area) || ref_area <= 0) {
    ref_area <- area_ss
  }
  duration <- windowT[2] - windowT[1]
  pixel_fun <- function(x, y) {
    w <- suppressWarnings(as.numeric(lookup(x, y)))
    w[!is.finite(w) | w < 0] <- 0
    duration * (mu / ref_area) * w
  }
  intensity_im <- spatstat.geom::as.im(pixel_fun, W = state_space)
  intensity_summary <- summary(intensity_im)
  as.numeric(
    intensity_summary$max + 0.05 * diff(intensity_summary$range)
  )
}

#' Simulate a bivariate ETAS process
#'
#' @param params Named list with 15 bivariate ETAS parameters.
#' @param windowT Numeric c(start, end).
#' @param windowS owin or c(xmin, xmax, ymin, ymax).
#' @param state_spaces Named list with "control" and "treated" owin objects.
#' @param m0 Reference magnitude.
#' @param beta_gr Gutenberg-Richter beta (NULL for resampling).
#' @param mag_pool Magnitude pool.
#' @param filtration Optional history data frame with x, y, t, mag, process_id.
#' @param covariate_lookup Optional function \code{f(x, y)} shared by both
#'   backgrounds, or a named list with \code{control} and \code{treated}
#'   functions. Each function should have spatial mean one on its state space.
#' @param bg_ref_areas Optional named list with \code{control} and/or
#'   \code{treated} reference areas for the background density. Each
#'   background is generated with spatial density \code{mu / ref_area *
#'   w(x,y)} over its state space, so \code{mu} is the expected per-unit-time
#'   count on a region of size \code{ref_area} rather than on the (possibly
#'   larger) simulation support. Defaults to each state space's own area,
#'   which reproduces the total-rate behaviour. Use this to extend a fitted
#'   background density onto a larger support (e.g. a control-everywhere
#'   counterfactual with \code{ref_area = area(fitted control region)}).
#' @param bg_lmax Optional named list of precomputed Poisson rejection bounds
#'   for the control and treated background lookup functions.
#' @param t_trunc Temporal truncation.
#' @param ... Ignored.
#' @return Data frame with x, y, t, mag, process_id, background columns.
#' @export
sim_etas_bivariate <- function(params,
                               windowT,
                               windowS,
                               state_spaces = NULL,
                               m0,
                               beta_gr = NULL,
                               mag_pool = NULL,
                               filtration = NULL,
                               covariate_lookup = NULL,
                               bg_ref_areas = NULL,
                               bg_lmax = NULL,
                               t_trunc = NULL,
                               ...) {
  if (!inherits(windowS, "owin")) {
    windowS <- owin(xrange = c(windowS[1], windowS[2]),
                    yrange = c(windowS[3], windowS[4]))
  }

  pv <- unlist(params)
  mu_0 <- pv[["mu_0"]]; mu_1 <- pv[["mu_1"]]

  x_min <- windowS$xrange[1]; x_max <- windowS$xrange[2]
  y_min <- windowS$yrange[1]; y_max <- windowS$yrange[2]

  use_gr <- !is.null(beta_gr) && beta_gr > 0
  pool <- if (!is.null(mag_pool)) as.numeric(mag_pool) else numeric(0)

  draw_magnitudes <- function(n) {
    if (n == 0) return(numeric(0))
    if (use_gr) m0 + rexp(n, rate = beta_gr)
    else if (length(pool) > 0) sample(pool, n, replace = TRUE)
    else rep(m0, n)
  }

  # Background for each process in its own state space. ref_area is the
  # density reference: intensity is mu / ref_area * w(x,y), so with the
  # default ref_area = area(ss) the total per-unit-time rate is mu, while a
  # smaller ref_area extends the same spatial density onto a larger support.
  gen_bg <- function(mu, ss, proc_id, lookup = NULL, ref_area = NULL,
                     lmax = NULL) {
    if (mu < 1e-10 || is.null(ss)) {
      return(list(x = numeric(0), y = numeric(0), t = numeric(0),
                  mag = numeric(0), proc = integer(0)))
    }
    area_ss <- spatstat.geom::area(ss)
    ref_area <- suppressWarnings(as.numeric(ref_area))
    if (length(ref_area) != 1L || !is.finite(ref_area) || ref_area <= 0) {
      ref_area <- area_ss
    }
    duration <- windowT[2] - windowT[1]
    if (is.function(lookup)) {
      pixel_fun <- function(x, y) {
        w <- suppressWarnings(as.numeric(lookup(x, y)))
        w[!is.finite(w) | w < 0] <- 0
        duration * (mu / ref_area) * w
      }
      lmax <- suppressWarnings(as.numeric(lmax))
      use_lmax <- length(lmax) == 1L && is.finite(lmax) && lmax >= 0
      bg_pp <- if (use_lmax) {
        spatstat.random::rpoispp(pixel_fun, lmax = lmax, win = ss)
      } else {
        spatstat.random::rpoispp(pixel_fun, win = ss)
      }
      bx <- bg_pp$x
      by <- bg_pp$y
      n_ok <- bg_pp$n
    } else {
      n_bg <- rpois(1, mu * duration * (area_ss / ref_area))
      if (n_bg == 0) {
        return(list(x = numeric(0), y = numeric(0), t = numeric(0),
                    mag = numeric(0), proc = integer(0)))
      }
      bg_pp <- spatstat.random::runifpoint(n = n_bg, win = ss)
      bx <- bg_pp$x
      by <- bg_pp$y
      n_ok <- bg_pp$n
    }
    if (n_ok == 0) {
      return(list(x = numeric(0), y = numeric(0), t = numeric(0),
                  mag = numeric(0), proc = integer(0)))
    }
    list(
      x = bx, y = by,
      t = sort(runif(n_ok, windowT[1], windowT[2])),
      mag = draw_magnitudes(n_ok),
      proc = rep(as.integer(proc_id), n_ok)
    )
  }

  ctrl_ss <- if (!is.null(state_spaces)) state_spaces[["control"]] else windowS
  treat_ss <- if (!is.null(state_spaces)) state_spaces[["treated"]] else windowS
  ctrl_lookup <- if (is.list(covariate_lookup)) {
    covariate_lookup[["control"]]
  } else {
    covariate_lookup
  }
  treat_lookup <- if (is.list(covariate_lookup)) {
    covariate_lookup[["treated"]]
  } else {
    covariate_lookup
  }

  ctrl_ref_area <- if (is.list(bg_ref_areas)) bg_ref_areas[["control"]] else NULL
  treat_ref_area <- if (is.list(bg_ref_areas)) bg_ref_areas[["treated"]] else NULL
  ctrl_lmax <- if (is.list(bg_lmax)) bg_lmax[["control"]] else NULL
  treat_lmax <- if (is.list(bg_lmax)) bg_lmax[["treated"]] else NULL

  bg0 <- gen_bg(mu_0, ctrl_ss, 0L, ctrl_lookup, ctrl_ref_area, ctrl_lmax)
  bg1 <- gen_bg(mu_1, treat_ss, 1L, treat_lookup, treat_ref_area, treat_lmax)

  p_x   <- c(bg0$x, bg1$x)
  p_y   <- c(bg0$y, bg1$y)
  p_t   <- c(bg0$t, bg1$t)
  p_mag <- c(bg0$mag, bg1$mag)
  p_proc <- c(bg0$proc, bg1$proc)
  p_bg  <- rep(TRUE, length(p_x))

  # Add filtration. Parents older than windowT[1] - t_trunc cannot place a
  # child in the observation window (out-of-window children never reproduce).
  if (!is.null(filtration)) {
    filtration <- .etas_trim_filtration_to_trunc(filtration, windowT, t_trunc)
    if (is.data.frame(filtration) && nrow(filtration) > 0) {
      f_proc <- if ("process_id" %in% names(filtration)) {
        as.integer(filtration$process_id)
      } else if ("inferred_process" %in% names(filtration)) {
        as.integer(filtration$inferred_process == "treated")
      } else {
        rep(0L, nrow(filtration))
      }
      p_x <- c(p_x, filtration$x)
      p_y <- c(p_y, filtration$y)
      p_t <- c(p_t, filtration$t)
      p_mag <- c(p_mag, if ("mag" %in% names(filtration)) filtration$mag else rep(m0, nrow(filtration)))
      p_proc <- c(p_proc, f_proc)
      p_bg <- c(p_bg, rep(TRUE, nrow(filtration)))
    }
  }

  if (length(p_t) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0), t = numeric(0),
                      mag = numeric(0), process_id = integer(0),
                      background = logical(0)))
  }

  # Sort by time
  ord <- order(p_t)
  p_x <- p_x[ord]; p_y <- p_y[ord]; p_t <- p_t[ord]
  p_mag <- p_mag[ord]; p_proc <- p_proc[ord]; p_bg <- p_bg[ord]

  # windowS is the modeled support, not merely its bounding rectangle. Remove
  # unsupported filtration points before they can generate descendants.
  supported <- spatstat.geom::inside.owin(p_x, p_y, windowS)
  p_x <- p_x[supported]; p_y <- p_y[supported]; p_t <- p_t[supported]
  p_mag <- p_mag[supported]; p_proc <- p_proc[supported]; p_bg <- p_bg[supported]

  any_excitation <- max(pv[["A_00"]], pv[["A_11"]], pv[["A_01"]], pv[["A_10"]]) > 1e-10

  if (any_excitation && length(p_t) > 0) {
    # Generate one generation at a time so children outside a non-rectangular
    # window are killed before they can enter the branching queue.
    generation <- list(
      x = p_x, y = p_y, t = p_t, mag = p_mag,
      process_id = as.integer(p_proc)
    )
    generation_index <- 0L
    repeat {
      generation_index <- generation_index + 1L
      if (generation_index > 10000L) {
        stop("Bivariate ETAS exceeded 10,000 offspring generations.")
      }
      children <- sim_etas_bivariate_children_cpp(
        parent_x = generation$x, parent_y = generation$y,
        parent_t = generation$t, parent_mag = generation$mag,
        parent_process = as.integer(generation$process_id),
        A_00 = pv[["A_00"]], alpha_m_00 = pv[["alpha_m_00"]],
        A_11 = pv[["A_11"]], alpha_m_11 = pv[["alpha_m_11"]],
        A_01 = pv[["A_01"]], alpha_m_01 = pv[["alpha_m_01"]],
        A_10 = pv[["A_10"]], alpha_m_10 = pv[["alpha_m_10"]],
        cc = pv[["c"]], p = pv[["p"]], D = pv[["D"]],
        gamma_par = pv[["gamma"]], q = pv[["q"]],
        m0 = m0, beta_gr = if (use_gr) beta_gr else -1.0,
        t_min = windowT[1], t_max = windowT[2],
        x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max,
        t_trunc = if (!is.null(t_trunc)) t_trunc else -1.0,
        mag_pool = pool,
        max_generations = 1L
      )
      if (length(children$x) < 1L) break

      child_supported <- spatstat.geom::inside.owin(
        children$x, children$y, windowS
      )
      if (!any(child_supported)) break
      generation <- lapply(children, function(x) x[child_supported])
      n_ch <- length(generation$x)
      p_x <- c(p_x, generation$x); p_y <- c(p_y, generation$y)
      p_t <- c(p_t, generation$t); p_mag <- c(p_mag, generation$mag)
      p_proc <- c(p_proc, generation$process_id)
      p_bg <- c(p_bg, rep(FALSE, n_ch))
    }
  }

  # Filter to window
  valid <- p_t >= windowT[1] & p_t <= windowT[2] &
    spatstat.geom::inside.owin(p_x, p_y, windowS)
  data.frame(
    x = p_x[valid], y = p_y[valid], t = p_t[valid],
    mag = p_mag[valid], process_id = p_proc[valid],
    background = p_bg[valid]
  )
}


#' Generate inhomogeneous bivariate ETAS on a partitioned domain
#'
#' Multi-region analogue of generate_inhomogeneous_etas for the bivariate
#' model. Backgrounds are generated per-process in their respective state
#' spaces, then a single joint BFS generates all offspring with cross-excitation.
#'
#' @param Omega Full state space (owin).
#' @param partition spatstat tess object.
#' @param time_window Numeric c(start, end).
#' @param partition_processes Character vector of process names per tile.
#' @param etas_bivariate_params Named list of 15 bivariate ETAS parameters.
#' @param m0 Reference magnitude.
#' @param beta_gr Gutenberg-Richter beta (NULL for resampling).
#' @param mag_pool Magnitude pool.
#' @param state_spaces Named list with "control" and "treated" owin.
#' @param filtration History data frame.
#' @param t_trunc Temporal truncation.
#' @param ... Ignored.
#' @return data.table with x, y, t, mag, process_id, background, location_process.
#' @export
generate_inhomogeneous_etas_bivariate <- function(
    Omega, partition, time_window, partition_processes,
    etas_bivariate_params, m0,
    beta_gr = NULL, mag_pool = NULL,
    state_spaces = NULL, filtration = NULL,
    t_trunc = NULL, ...) {

  if (!inherits(Omega, "owin")) Omega <- as.owin(Omega)

  if (is.null(state_spaces)) {
    state_spaces <- list(
      control = as.owin(partition[partition_processes != "treated"]),
      treated = as.owin(partition[partition_processes == "treated"])
    )
  }
  support_window <- .etas_union_support(state_spaces, Omega)

  result <- sim_etas_bivariate(
    params = etas_bivariate_params,
    windowT = time_window,
    windowS = support_window,
    state_spaces = state_spaces,
    m0 = m0,
    beta_gr = beta_gr,
    mag_pool = mag_pool,
    filtration = filtration,
    t_trunc = t_trunc
  )

  if (nrow(result) == 0) {
    result$location_process <- character(0)
    result$process <- character(0)
    return(data.table::as.data.table(result))
  }

  tile_idx <- as.integer(tileindex(
    result$x, result$y, partition, close.gaps = FALSE
  ))
  supported <- !is.na(tile_idx)
  if (!all(supported)) {
    result <- result[supported, , drop = FALSE]
    tile_idx <- tile_idx[supported]
  }
  if (nrow(result) == 0L) {
    result$tile_index <- integer(0)
    result$location_process <- character(0)
    result$process <- character(0)
    return(data.table::as.data.table(result))
  }
  result$tile_index <- tile_idx
  result$location_process <- partition_processes[
    pmin(pmax(tile_idx, 1), length(partition_processes))
  ]
  result$process <- ifelse(result$process_id == 0, "control", "treated")

  data.table::as.data.table(result)
}


#' Initialize bivariate ETAS parameters from independent fits
#'
#' Maps two independent ETAS parameter sets (control and treated) into a
#' single bivariate parameter vector, with cross-excitation initialized
#' near zero.
#'
#' @param ctrl_params Named list of control ETAS params (mu, A, alpha_m, ...).
#' @param treat_params Named list of treated ETAS params.
#' @param cross_A Initial cross-excitation productivity (default 0.01).
#' @param cross_alpha_m Initial cross-excitation magnitude scaling (NULL = mean of self-terms).
#' @return Named numeric vector of 15 bivariate parameters.
#' @export
init_bivariate_from_independent <- function(ctrl_params,
                                            treat_params,
                                            cross_A = 0.01,
                                            cross_alpha_m = NULL) {
  cp <- unlist(ctrl_params)
  tp <- unlist(treat_params)

  if (is.null(cross_alpha_m)) {
    cross_alpha_m <- mean(c(cp[["alpha_m"]], tp[["alpha_m"]]))
  }

  c(
    mu_0 = unname(cp[["mu"]]),
    mu_1 = unname(tp[["mu"]]),
    A_00 = unname(cp[["A"]]),
    alpha_m_00 = unname(cp[["alpha_m"]]),
    A_11 = unname(tp[["A"]]),
    alpha_m_11 = unname(tp[["alpha_m"]]),
    A_01 = cross_A,
    alpha_m_01 = cross_alpha_m,
    A_10 = cross_A,
    alpha_m_10 = cross_alpha_m,
    c = unname(cp[["c"]]),
    p = unname(cp[["p"]]),
    D = unname(cp[["D"]]),
    gamma = unname(cp[["gamma"]]),
    q = unname(cp[["q"]])
  )
}
