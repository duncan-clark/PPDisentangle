# ============================================================================
# ETAS (Epidemic-Type Aftershock Sequence) model
#
# Implements the spatio-temporal ETAS model of Ogata (1988) and
# Zhuang et al. (2002) with:
#   - Omori-Utsu power-law temporal kernel
#   - Isotropic power-law spatial kernel with magnitude-dependent spread
#   - Magnitude-dependent productivity (Gutenberg-Richter)
#
# The ETAS parameter vector has 8 free components:
#   (mu, A, alpha_m, c, p, D, gamma, q)
# plus a fixed reference magnitude m0 (typically min observed magnitude).
#
# See individual function documentation for the full mathematical
# specification.
# ============================================================================

# Canonical ordering and names for the 8-element ETAS parameter vector.
.etas_par_names <- c("mu", "A", "alpha_m", "c", "p", "D", "gamma", "q")

.etas_param_vector <- function(params, context = "ETAS parameters") {
  pv <- if (is.list(params)) unlist(params, use.names = TRUE) else params
  if (is.null(names(pv))) {
    if (length(pv) != length(.etas_par_names)) {
      stop(context, " must contain exactly ", length(.etas_par_names),
           " values in canonical ETAS order.")
    }
    names(pv) <- .etas_par_names
  }
  missing_names <- setdiff(.etas_par_names, names(pv))
  if (length(missing_names) > 0L) {
    stop(context, " is missing: ", paste(missing_names, collapse = ", "), ".")
  }
  out <- suppressWarnings(as.numeric(pv[.etas_par_names]))
  names(out) <- .etas_par_names
  if (any(!is.finite(out))) {
    stop(context, " must contain finite values for all ETAS parameters.")
  }
  out
}

.etas_resolve_beta_gr <- function(beta_gr, realiz = NULL, m0 = NULL) {
  beta_eff <- suppressWarnings(as.numeric(beta_gr))
  beta_eff <- beta_eff[is.finite(beta_eff) & beta_eff > 0]
  if (length(beta_eff)) return(beta_eff[[1]])
  if (is.null(realiz) || is.null(realiz$mag)) return(NA_real_)
  if (is.null(m0) || !is.finite(m0)) {
    m0 <- suppressWarnings(min(as.numeric(realiz$mag), na.rm = TRUE))
  }
  mag_delta <- as.numeric(realiz$mag) - as.numeric(m0)
  mag_delta <- mag_delta[is.finite(mag_delta) & mag_delta > 0]
  if (!length(mag_delta)) return(NA_real_)
  1 / mean(mag_delta)
}

.etas_univ_branching_ratio <- function(params, beta_gr) {
  pv <- if (is.list(params)) unlist(params) else params
  beta_eff <- .etas_resolve_beta_gr(beta_gr)
  A <- suppressWarnings(as.numeric(pv[["A"]]))
  alpha <- suppressWarnings(as.numeric(pv[["alpha_m"]]))
  gap <- beta_eff - alpha
  if (!is.finite(A) || A < 0 || !is.finite(gap) || gap <= 1e-8) return(Inf)
  A * beta_eff / gap
}

#' Resolve admissible interval for magnitude productivity alpha_m
#'
#' Default lower bound is 0 (larger events more productive). Upper bound is
#' the hard subcritical ceiling \code{beta - gap_min}.
#' @keywords internal
.etas_alpha_bounds <- function(beta_gr,
                               gap_min = 1e-4,
                               alpha_m_lower_bound = 0) {
  beta_eff <- .etas_resolve_beta_gr(beta_gr)
  gap_min <- suppressWarnings(as.numeric(gap_min))
  if (length(gap_min) != 1L || !is.finite(gap_min) || is.na(gap_min) || gap_min < 0) {
    gap_min <- 1e-4
  }
  alpha_lo <- suppressWarnings(as.numeric(alpha_m_lower_bound))
  if (length(alpha_lo) != 1L || is.na(alpha_lo)) alpha_lo <- 0
  if (!is.finite(beta_eff)) {
    return(list(lo = alpha_lo, hi = Inf, beta = beta_eff, gap_min = gap_min))
  }
  alpha_hi <- beta_eff - gap_min
  if (is.finite(alpha_lo) && alpha_lo >= alpha_hi) {
    # Pathological: keep a tiny open interval just below the subcritical cap.
    alpha_lo <- alpha_hi - 1e-3
  }
  list(lo = alpha_lo, hi = alpha_hi, beta = beta_eff, gap_min = gap_min)
}

#' Map alpha_m from natural scale to unconstrained optimiser coordinate
#' @keywords internal
.etas_alpha_to_opt <- function(alpha, lo, hi) {
  alpha <- suppressWarnings(as.numeric(alpha))
  lo <- suppressWarnings(as.numeric(lo))
  hi <- suppressWarnings(as.numeric(hi))
  if (!length(alpha)) return(numeric(0))
  out <- rep(0, length(alpha))
  for (i in seq_along(alpha)) {
    a <- alpha[[i]]
    if (!is.finite(a) || !is.finite(lo) || !is.finite(hi) || hi <= lo) {
      out[[i]] <- 0
      next
    }
    # One-sided (legacy): alpha = hi - exp(u) when lo = -Inf
    if (!is.finite(lo) && is.finite(hi)) {
      out[[i]] <- log(max(hi - a, 1e-12))
      next
    }
    u <- (a - lo) / (hi - lo)
    u <- min(max(u, 1e-12), 1 - 1e-12)
    out[[i]] <- stats::qlogis(u)
  }
  out
}

#' Map unconstrained optimiser coordinate back to alpha_m
#' @keywords internal
.etas_opt_to_alpha <- function(u, lo, hi) {
  u <- suppressWarnings(as.numeric(u))
  lo <- suppressWarnings(as.numeric(lo))
  hi <- suppressWarnings(as.numeric(hi))
  if (!length(u)) return(numeric(0))
  out <- rep(NA_real_, length(u))
  for (i in seq_along(u)) {
    ui <- u[[i]]
    if (!is.finite(ui)) {
      out[[i]] <- if (is.finite(lo)) lo + 1e-8 else NA_real_
      next
    }
    if (!is.finite(lo) && is.finite(hi)) {
      out[[i]] <- hi - exp(ui)
      next
    }
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
      out[[i]] <- ui
      next
    }
    out[[i]] <- lo + (hi - lo) * stats::plogis(ui)
  }
  out
}

.etas_project_subcritical <- function(params,
                                      beta_gr,
                                      gap_min = 1e-4,
                                      eta_max = 0.98,
                                      alpha_m_lower_bound = 0) {
  pv <- if (is.list(params)) unlist(params) else params
  if (is.null(names(pv))) names(pv) <- .etas_par_names
  bounds <- .etas_alpha_bounds(
    beta_gr, gap_min = gap_min, alpha_m_lower_bound = alpha_m_lower_bound
  )
  beta_eff <- bounds$beta
  if (!is.finite(beta_eff)) return(pv)
  eta_max <- suppressWarnings(as.numeric(eta_max))
  if (!is.finite(eta_max) || eta_max <= 0) eta_max <- 0.98
  if (is.finite(pv[["alpha_m"]])) {
    a <- pv[["alpha_m"]]
    if (is.finite(bounds$lo)) a <- max(a, bounds$lo + 1e-8)
    if (is.finite(bounds$hi)) a <- min(a, bounds$hi - 1e-4)
    pv[["alpha_m"]] <- a
  }
  eta <- .etas_univ_branching_ratio(pv, beta_eff)
  if (is.finite(eta) && eta > eta_max && eta > 0 && is.finite(pv[["A"]])) {
    pv[["A"]] <- pv[["A"]] * (eta_max / eta) * (1 - 1e-6)
  }
  pv
}

.etas_scale_to_target_eta <- function(params, beta_gr, target_eta) {
  pv <- if (is.list(params)) unlist(params) else params
  if (is.null(names(pv))) names(pv) <- .etas_par_names
  target_eta <- suppressWarnings(as.numeric(target_eta))
  if (length(target_eta) != 1L || !is.finite(target_eta) || target_eta <= 0) {
    return(pv)
  }
  eta <- .etas_univ_branching_ratio(pv, beta_gr)
  if (!is.finite(eta) || eta <= 0 || !is.finite(pv[["A"]])) return(pv)
  pv[["A"]] <- as.numeric(pv[["A"]]) * (target_eta / eta)
  pv
}

.etas_soft_barrier_controls <- function(cap,
                                        start = NULL,
                                        weight = 5000,
                                        power = 4) {
  cap <- suppressWarnings(as.numeric(cap))
  if (length(cap) != 1L || !is.finite(cap) || cap <= 0) cap <- 0.98
  start <- suppressWarnings(as.numeric(start))
  if (length(start) != 1L || !is.finite(start)) {
    start <- min(0.90, max(0.5, cap * 0.92))
  }
  weight <- suppressWarnings(as.numeric(weight))
  if (length(weight) != 1L || !is.finite(weight) || weight <= 0) weight <- 5000
  power <- suppressWarnings(as.numeric(power))
  if (length(power) != 1L || !is.finite(power) || power <= 0) power <- 4
  list(
    stability_barrier_start = start,
    stability_barrier_weight = weight,
    stability_barrier_power = power
  )
}

.etas_raw_loglik <- function(params, ll_args) {
  args <- ll_args
  args$params <- params
  # Evaluate the statistical objective without the hard branching cliff or
  # soft barrier used only to guide Nelder-Mead.
  args$max_branching_ratio <- Inf
  args$stability_barrier_weight <- 0
  tryCatch(
    as.numeric(do.call(loglik_etas, args)),
    error = function(e) -1e15
  )
}

.etas_polish_productivity_A <- function(params,
                                        ll_args,
                                        beta_gr,
                                        eta_max = 0.98,
                                        eta_lo = 0.05) {
  pv <- if (is.list(params)) unlist(params) else params
  if (is.null(names(pv))) names(pv) <- .etas_par_names
  eta_max <- suppressWarnings(as.numeric(eta_max))
  if (!is.finite(eta_max) || eta_max <= 0) eta_max <- 0.98
  eta_lo <- suppressWarnings(as.numeric(eta_lo))
  if (!is.finite(eta_lo) || eta_lo <= 0) eta_lo <- 0.05
  hi <- eta_max * (1 - 1e-6)
  if (!(hi > eta_lo)) return(list(par = pv, value = .etas_raw_loglik(pv, ll_args), polished = FALSE))
  obj <- function(eta) {
    q <- .etas_scale_to_target_eta(pv, beta_gr, eta)
    .etas_raw_loglik(q, ll_args)
  }
  opt <- tryCatch(
    stats::optimize(obj, interval = c(eta_lo, hi), maximum = TRUE),
    error = function(e) NULL
  )
  if (is.null(opt) || !is.finite(opt$objective)) {
    return(list(par = pv, value = .etas_raw_loglik(pv, ll_args), polished = FALSE))
  }
  best <- .etas_scale_to_target_eta(pv, beta_gr, opt$maximum)
  list(par = best, value = as.numeric(opt$objective), polished = TRUE)
}

#' Log-likelihood for a spatio-temporal ETAS process (C++ accelerated)
#'
#' Evaluates the log-likelihood of the ETAS conditional intensity model
#' on a marked spatio-temporal point pattern.
#'
#' @section Model:
#' The conditional intensity is
#' \deqn{\lambda(t,x,y) = \frac{\mu}{|S|}\,W(x,y) +
#'   \sum_{t_j < t} \kappa(m_j)\,g(t - t_j)\,f(x-x_j, y-y_j \mid m_j)}
#'
#' \describe{
#'   \item{Productivity}{\eqn{\kappa(m) = A\,\exp(\alpha_m\,(m - m_0))}}
#'   \item{Omori-Utsu temporal kernel}{
#'     \eqn{g(\Delta t) = \frac{p-1}{c}\bigl(1 + \Delta t / c\bigr)^{-p}},
#'     \eqn{p > 1,\; c > 0}}
#'   \item{Power-law spatial kernel (Zhuang et al. 2002)}{
#'     \eqn{f(x,y \mid m) = \frac{q-1}{\pi\,d(m)}
#'       \bigl(1 + r^2/d(m)\bigr)^{-q}},
#'     \eqn{d(m) = D\exp(\gamma(m - m_0)),\; q > 1}}
#' }
#'
#' The compensator (integral of \eqn{\lambda}) uses the infinite-plane
#' approximation for the spatial kernel (consistent with the Hawkes
#' implementation in this package).
#'
#' @param params  Named numeric vector or list with elements
#'   \code{mu, A, alpha_m, c, p, D, gamma, q}.
#' @param realiz  Data frame with columns \code{x, y, t, mag}
#'   (and optionally a background-rate column, default \code{"W"}).
#' @param windowT Numeric vector \code{c(start, end)} for the time window.
#' @param windowS An \code{owin} object (or convertible) for the spatial
#'   observation window.
#' @param m0  Reference (cutoff) magnitude.  If \code{NULL} (default),
#'   set to \code{min(realiz$mag)}.
#' @param zero_background_region  Optional \code{owin} where the background
#'   intensity is zero (events may still be triggered there).
#' @param background_rate_var  Column name in \code{realiz} for
#'   inhomogeneous background weights (default \code{"W"}).
#' @param precomp  Optional list from \code{precompute_loglik_args} to skip
#'   redundant area / inside.owin calculations.
#' @param t_trunc  Temporal truncation horizon.  \code{NULL} = no truncation.
#' @param history Optional pre-window event data with columns
#'   \code{x}, \code{y}, \code{t}, and \code{mag}. These events contribute as
#'   triggering parents, but their own intensities are not included in the
#'   likelihood. The compensator is integrated only over \code{windowT}.
#' @param ...  Additional arguments (ignored).
#' @return Scalar log-likelihood value.
#' @export
loglik_etas <- function(params,
                        realiz,
                        windowT,
                        windowS,
                        m0 = NULL,
                        zero_background_region = NULL,
                        background_rate_var = "W",
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
                        max_branching_ratio = Inf,
                        stability_barrier_start = 0.95,
                        stability_barrier_weight = 0,
                        stability_barrier_power = 4,
                        precomp = NULL,
                        t_trunc = NULL,
                        history = NULL,
                        ...) {
  # --- Parse parameters ---
  if (is.list(params)) {
    mu      <- params$mu
    A       <- params$A
    alpha_m <- params$alpha_m
    cc      <- params[["c"]]
    p       <- params$p
    D       <- params$D
    gamma_p <- params$gamma
    q       <- params$q
  } else {
    mu <- params[1]; A <- params[2]; alpha_m <- params[3]; cc <- params[4]
    p  <- params[5]; D <- params[6]; gamma_p <- params[7]; q  <- params[8]
  }

  if (!all(is.finite(c(mu, A, alpha_m, cc, p, D, gamma_p, q)))) return(-1e15)
  if (!all(is.finite(realiz$t)) || !all(is.finite(realiz$x)) || !all(is.finite(realiz$y)) || !all(is.finite(realiz$mag))) {
    return(-1e15)
  }

  # --- Sorting ---
  if (is.unsorted(realiz$t)) realiz <- realiz[order(realiz$t), ]

  # --- Time window filter ---
  t_idx <- realiz$t >= windowT[1] & realiz$t <= windowT[2]
  if (!all(t_idx)) realiz <- realiz[t_idx, ]
  n <- nrow(realiz)
  if (n == 0) return(-1e15)

  # --- Optional pre-window parent history ---
  use_history <- !is.null(history)
  history_data <- realiz[0, c("x", "y", "t", "mag"), drop = FALSE]
  if (use_history) {
    history <- as.data.frame(history)
    required_history_cols <- c("x", "y", "t", "mag")
    if (nrow(history) > 0L &&
        !all(required_history_cols %in% names(history))) {
      return(-1e15)
    }
    if (nrow(history) > 0L) {
      history_data <- history[, required_history_cols, drop = FALSE]
      history_finite <- vapply(history_data, function(z) all(is.finite(z)), logical(1))
      if (!all(history_finite)) return(-1e15)
      history_data <- history_data[
        history_data$t < windowT[1], , drop = FALSE
      ]
      trunc_value <- suppressWarnings(as.numeric(t_trunc))
      if (length(trunc_value) == 1L && is.finite(trunc_value) &&
          !is.na(trunc_value) && trunc_value > 0) {
        history_data <- history_data[
          history_data$t >= windowT[1] - trunc_value, , drop = FALSE
        ]
      }
      history_data <- history_data[order(history_data$t), , drop = FALSE]
    }
  }

  # --- Parameter bounds ---
  if (min(mu, A, cc, D) < 0 || p <= 1 || q <= 1 || gamma_p < 0) return(-1e15)
  p_min <- suppressWarnings(as.numeric(p_lower_bound))
  if (length(p_min) != 1L || !is.finite(p_min) || is.na(p_min)) p_min <- 2.001
  q_min <- suppressWarnings(as.numeric(q_lower_bound))
  if (length(q_min) != 1L || !is.finite(q_min) || is.na(q_min)) q_min <- 1.501
  if (isTRUE(enforce_finite_trigger_moments) && (p <= p_min || q <= q_min)) return(-1e15)

  # --- Reference magnitude ---
  if (is.null(m0)) m0 <- min(realiz$mag)
  beta_eff <- .etas_resolve_beta_gr(beta_gr, realiz = realiz, m0 = m0)
  alpha_lo <- suppressWarnings(as.numeric(alpha_m_lower_bound))
  if (length(alpha_lo) != 1L || is.na(alpha_lo)) alpha_lo <- 0
  if (is.finite(alpha_lo) && is.finite(alpha_m) && alpha_m <= alpha_lo) return(-1e15)
  if (isTRUE(enforce_alpha_subcritical)) {
    gap_min <- suppressWarnings(as.numeric(alpha_beta_gap_min))
    if (length(gap_min) != 1L || !is.finite(gap_min) || is.na(gap_min) || gap_min < 0) gap_min <- 1e-4
    if (!is.finite(beta_eff) || is.na(beta_eff) || alpha_m >= (beta_eff - gap_min)) return(-1e15)
  }
  eta_max <- suppressWarnings(as.numeric(max_branching_ratio))
  if (length(eta_max) == 1L && is.finite(eta_max) && eta_max > 0 &&
      is.finite(beta_eff)) {
    eta_hat <- .etas_univ_branching_ratio(
      c(A = as.numeric(A), alpha_m = as.numeric(alpha_m)),
      beta_eff
    )
    if (!is.finite(eta_hat) || eta_hat >= eta_max) return(-1e15)
  }

  # --- Background weights ---
  W_vec <- if (!is.null(background_rate_var) &&
               background_rate_var %in% names(realiz)) {
    realiz[[background_rate_var]]
  } else {
    rep(1.0, n)
  }

  # --- Spatial areas and zero-background handling ---
  if (!is.null(precomp)) {
    active_area <- precomp$active_area
    if (!is.null(precomp$in_zero_bg)) {
      in_zero <- precomp$in_zero_bg
      if (is.logical(in_zero) && length(in_zero) == n) {
        W_vec[in_zero] <- 0
      }
    }
  } else {
    total_area <- spatstat.geom::area(as.owin(windowS))
    if (!is.null(zero_background_region)) {
      zero_area <- spatstat.geom::area(as.owin(zero_background_region))
      active_area <- total_area - zero_area
      if (active_area <= 0) warning("Zero-background region covers entire window!")
      in_zero_bg <- inside.owin(realiz[, c("x", "y")],
                                w = zero_background_region)
      W_vec[in_zero_bg] <- 0
    } else {
      active_area <- total_area
    }
  }

  tval <- windowT[2] - windowT[1]

  # --- Call C++ ---
  if (use_history) {
    parents <- rbind(
      history_data,
      realiz[, c("x", "y", "t", "mag"), drop = FALSE]
    )
    parents <- parents[order(parents$t), , drop = FALSE]
    loglik <- etas_loglik_inhom_filtration_cpp(
      post_t     = realiz$t,
      post_x     = realiz$x,
      post_y     = realiz$y,
      W_val      = W_vec,
      parent_t   = parents$t,
      parent_x   = parents$x,
      parent_y   = parents$y,
      parent_mag = parents$mag,
      mu         = mu,
      A          = A,
      alpha_m    = alpha_m,
      cc         = cc,
      p          = p,
      D           = D,
      gamma_par  = gamma_p,
      q           = q,
      m0          = m0,
      areaS       = active_area,
      t_start     = windowT[1],
      t_end       = windowT[2],
      t_trunc     = if (!is.null(t_trunc)) t_trunc else -1.0
    )
  } else {
    loglik <- etas_loglik_inhom_cpp(
      t         = realiz$t - windowT[1],
      x         = realiz$x,
      y         = realiz$y,
      mag       = realiz$mag,
      W_val     = W_vec,
      mu        = mu,
      A         = A,
      alpha_m   = alpha_m,
      cc        = cc,
      p         = p,
      D         = D,
      gamma_par = gamma_p,
      q         = q,
      m0        = m0,
      areaS     = active_area,
      t_max     = tval,
      t_trunc   = if (!is.null(t_trunc)) t_trunc else -1.0
    )
  }
  if (!is.finite(loglik)) return(-1e15)

  # Hybrid boundary handling: keep hard admissibility, plus a smooth interior
  # penalty in a narrow band above p_min/q_min and alpha-beta gap_min.
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
      alpha_gap <- beta_soft - as.numeric(alpha_m)
      gap_excess <- max(0, (gap_min + soft_gap) - alpha_gap)
      if (gap_excess > 0) loglik <- loglik - soft_weight * (gap_excess ^ soft_power)
    }
  }

  # Smooth stability barrier on the univariate ETAS branching ratio eta
  # (off by default, weight 0; the hard max_branching_ratio cap is the
  # guardrail).
  barrier_weight <- suppressWarnings(as.numeric(stability_barrier_weight))
  if (length(barrier_weight) != 1L || !is.finite(barrier_weight) || is.na(barrier_weight) || barrier_weight <= 0) barrier_weight <- 0
  if (barrier_weight > 0 && is.finite(beta_eff)) {
    eta <- .etas_univ_branching_ratio(
      c(A = as.numeric(A), alpha_m = as.numeric(alpha_m)),
      beta_eff
    )
    if (!is.finite(eta)) return(-1e15)
    barrier_start <- suppressWarnings(as.numeric(stability_barrier_start))
    if (length(barrier_start) != 1L || !is.finite(barrier_start) || is.na(barrier_start)) barrier_start <- 0.95
    barrier_power <- suppressWarnings(as.numeric(stability_barrier_power))
    if (length(barrier_power) != 1L || !is.finite(barrier_power) || is.na(barrier_power) || barrier_power <= 0) barrier_power <- 2
    excess <- max(0, eta - barrier_start)
    if (excess > 0) {
      loglik <- loglik - barrier_weight * (excess ^ barrier_power)
    }
  }
  return(loglik)
}


#' Fit a spatio-temporal ETAS process via MLE
#'
#' Maximum-likelihood estimation for the 8-parameter ETAS model using
#' \code{optim}.  Supports profile likelihood via \code{fixed_params}
#' (e.g.\ fixing structural parameters \code{c, p, D, gamma, q} while
#' estimating \code{mu, A, alpha_m}).
#'
#' @inheritParams loglik_etas
#' @param params_init  Initial parameter values (named list or numeric vector
#'   of length 8 in canonical order:
#'   \code{mu, A, alpha_m, c, p, D, gamma, q}).
#' @param m0  Reference magnitude (\code{NULL} = auto).
#' @param trace  Trace level for \code{optim} (0 = silent).
#' @param maxit  Maximum iterations.
#' @param method  Optimisation method (default \code{"Nelder-Mead"};
#'   \code{"L-BFGS-B"} supported with explicit bounds).
#' @param lower,upper  Bounds for \code{L-BFGS-B}.
#' @param fixed_params  Named list of parameters to hold fixed during
#'   optimisation.  E.g.\ \code{list(c = 0.01, p = 1.2, q = 1.5)}
#'   optimises only the remaining parameters.
#' @param log_transform  Logical. If \code{TRUE}, optimise in
#'   log-space for strictly positive parameters (\code{mu, A, c, D}).
#'   This flattens ridges in the likelihood surface and enforces
#'   positivity without explicit bounds.  Default \code{FALSE}.
#' @param init_decluster  Logical. If \code{TRUE}, override the
#'   \code{mu} and \code{A} starting values with moment-based estimates
#'   derived from the expected branching ratio.  Uses the supplied
#'   \code{alpha_m} init and the magnitude distribution to decompose
#'   observed event counts into background and triggered components.
#'   Default \code{FALSE}.
#' @param hard_subcritical If TRUE (default), constrain
#'   \code{alpha_m_lower_bound < alpha_m < beta_gr - alpha_beta_gap_min},
#'   initialize inside \code{eta < max_branching_ratio}, and return a projected
#'   stable fit. The default lower bound is 0 (magnitude-increasing productivity).
#' @param ...  Passed through to \code{loglik_etas}.
#' @return An \code{optim} result list.  \code{$par} is always a length-8
#'   named vector (including fixed values) on the *original* scale.
#'   An additional element \code{$m0} stores the reference magnitude used.
#' @export
fit_etas <- function(params_init,
                     realiz,
                     windowT,
                     windowS,
                     m0 = NULL,
                     trace = 0,
                     maxit = 5000,
                     method = "Nelder-Mead",
                     lower = NULL,
                     upper = NULL,
                     fixed_params = NULL,
                     t_trunc = NULL,
                     log_transform = TRUE,
                     init_decluster = FALSE,
                     hard_subcritical = TRUE,
                     soft_branching_barrier = TRUE,
                     polish_productivity = TRUE,
                     interior_restart = TRUE,
                     ...) {
  dots <- list(...)
  dots$precomp <- NULL
  if (inherits(params_init, "list")) params_init <- unlist(params_init)

  all_names <- .etas_par_names  # mu, A, alpha_m, c, p, D, gamma, q
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

  # Moment-based declustering initialization:
  # Given alpha_m and magnitudes, the expected branching ratio is
  #   n = A * mean(exp(alpha_m * (mag - m0)))
  # Total events N in time T satisfy N = mu*T / (1 - n), so
  #   mu = N * (1-n) / T  and  A = n / mean(exp(alpha_m*(mag-m0)))
  if (init_decluster && !is.null(realiz$mag)) {
    N <- nrow(realiz)
    Tlen <- windowT[2] - windowT[1]
    am_init <- full_init["alpha_m"]
    mean_prod <- mean(exp(am_init * (realiz$mag - m0)))
    # Target a branching ratio of ~0.5 (reasonable subcritical default),
    # but cap at what the data can support
    n_hat <- min(0.5, max(0.05, 1 - full_init["mu"] * Tlen / N))
    A_init <- n_hat / mean_prod
    mu_init <- N * (1 - n_hat) / Tlen
    if (is.finite(mu_init) && mu_init > 0 && !"mu" %in% names(fixed_params))
      full_init["mu"] <- mu_init
    if (is.finite(A_init) && A_init > 0 && !"A" %in% names(fixed_params))
      full_init["A"] <- A_init
  }

  # Keep optimiser starts inside constrained region when constraints are active.
  p_min <- suppressWarnings(as.numeric(dots$p_lower_bound))
  if (length(p_min) != 1L || !is.finite(p_min) || is.na(p_min)) p_min <- 2.001
  q_min <- suppressWarnings(as.numeric(dots$q_lower_bound))
  if (length(q_min) != 1L || !is.finite(q_min) || is.na(q_min)) q_min <- 1.501
  if (!isFALSE(dots$enforce_finite_trigger_moments)) {
    if (is.finite(full_init["p"])) full_init["p"] <- max(full_init["p"], p_min + 1e-3)
    if (is.finite(full_init["q"])) full_init["q"] <- max(full_init["q"], q_min + 1e-3)
  }
  beta_eff <- .etas_resolve_beta_gr(dots$beta_gr, realiz = realiz, m0 = m0)
  gap_min <- suppressWarnings(as.numeric(dots$alpha_beta_gap_min))
  if (length(gap_min) != 1L || !is.finite(gap_min) || is.na(gap_min) || gap_min < 0) gap_min <- 1e-4
  eta_max <- suppressWarnings(as.numeric(dots$max_branching_ratio))
  if (length(eta_max) != 1L || !is.finite(eta_max) || eta_max <= 0) eta_max <- 0.98
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
  restart_eta <- suppressWarnings(as.numeric(dots$interior_restart_eta))
  if (length(restart_eta) != 1L || !is.finite(restart_eta) || restart_eta <= 0) {
    restart_eta <- min(0.70, eta_max * init_margin)
  }
  dots$alpha_beta_gap_min <- gap_min
  dots$alpha_m_lower_bound <- alpha_lo
  alpha_bounds <- .etas_alpha_bounds(
    beta_eff, gap_min = gap_min, alpha_m_lower_bound = alpha_lo
  )
  # Start strictly inside the hard cap so Nelder-Mead is not born on the cliff.
  init_eta_cap <- eta_max * init_margin
  if (isTRUE(hard_subcritical) || !isFALSE(dots$enforce_alpha_subcritical) ||
      is.finite(alpha_lo)) {
    if (is.finite(beta_eff)) {
      full_init <- .etas_project_subcritical(
        full_init, beta_eff, gap_min = gap_min, eta_max = init_eta_cap,
        alpha_m_lower_bound = alpha_lo
      )
    }
  }

  # Soft barrier during NM; project back to eta_max after optim.
  use_soft_barrier <- isTRUE(soft_branching_barrier) && is.finite(eta_max)
  barrier_ctrl <- .etas_soft_barrier_controls(
    eta_max,
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
  dots_opt <- dots
  if (use_soft_barrier) {
    dots_opt$max_branching_ratio <- Inf
    dots_opt$stability_barrier_start <- barrier_ctrl$stability_barrier_start
    dots_opt$stability_barrier_weight <- barrier_ctrl$stability_barrier_weight
    dots_opt$stability_barrier_power <- barrier_ctrl$stability_barrier_power
  } else {
    dots_opt$max_branching_ratio <- eta_max
  }

  fixed_idx <- if (!is.null(fixed_params)) match(names(fixed_params), all_names) else integer(0)
  free_idx  <- setdiff(seq_along(all_names), fixed_idx)
  free_names <- all_names[free_idx]
  A_is_free <- "A" %in% free_names

  # Indices (within the free vector) of params that must be positive
  log_idx <- if (log_transform) which(free_names %in% c("mu", "A", "c", "D")) else integer(0)
  alpha_free_pos <- if ((isTRUE(hard_subcritical) || is.finite(alpha_lo)) &&
                          is.finite(beta_eff)) {
    which(free_names == "alpha_m")
  } else {
    integer(0)
  }

  finite_rows <- is.finite(realiz$t) & is.finite(realiz$x) & is.finite(realiz$y) & is.finite(realiz$mag)
  if (!all(finite_rows)) realiz <- realiz[finite_rows, , drop = FALSE]
  realiz <- realiz[order(realiz$t), , drop = FALSE]
  etas_precomp <- precompute_loglik_args(
    realiz,
    windowS,
    if ("zero_background_region" %in% names(dots)) dots$zero_background_region else NULL
  )
  precomp_list <- list(
    active_area = etas_precomp$active_area,
    in_zero_bg = etas_precomp$in_zero_bg_all
  )
  ll_base_args <- c(
    list(
      realiz = realiz,
      windowT = windowT,
      windowS = windowS,
      m0 = m0,
      precomp = precomp_list,
      t_trunc = t_trunc
    ),
    dots
  )
  # Raw LL for selection/polish ignores the NM soft barrier and hard cliff.
  ll_base_args$max_branching_ratio <- Inf
  ll_base_args$stability_barrier_weight <- 0

  opt_to_natural <- function(opt_par) {
    free_par <- opt_par
    if (length(log_idx) > 0) free_par[log_idx] <- exp(opt_par[log_idx])
    if (length(alpha_free_pos) > 0) {
      free_par[alpha_free_pos] <- .etas_opt_to_alpha(
        opt_par[alpha_free_pos], alpha_bounds$lo, alpha_bounds$hi
      )
    }
    free_par
  }
  natural_to_opt <- function(free_par) {
    opt_par <- free_par
    if (length(log_idx) > 0) {
      pos <- pmax(as.numeric(free_par[log_idx]), 1e-12)
      opt_par[log_idx] <- log(pos)
    }
    if (length(alpha_free_pos) > 0) {
      opt_par[alpha_free_pos] <- .etas_alpha_to_opt(
        free_par[alpha_free_pos], alpha_bounds$lo, alpha_bounds$hi
      )
    }
    opt_par
  }

  profile_fn <- function(opt_par, ...) {
    free_par <- opt_to_natural(opt_par)
    par8 <- full_init
    par8[free_idx] <- free_par
    ll_val <- tryCatch(
      loglik_etas(par8, ...),
      error = function(e) -1e15
    )
    if (!is.finite(ll_val) || is.na(ll_val)) return(-1e15)
    as.numeric(ll_val)
  }

  run_nm <- function(start_natural) {
    opt_args <- c(
      list(
        par     = natural_to_opt(start_natural[free_idx]),
        fn      = profile_fn,
        method  = method,
        control = list(fnscale = -1, trace = trace, maxit = maxit),
        realiz  = realiz,
        windowT = windowT,
        windowS = windowS,
        m0      = m0,
        precomp = precomp_list,
        t_trunc = t_trunc
      ),
      dots_opt
    )
    if (method == "L-BFGS-B" && !is.null(lower) && !is.null(upper)) {
      opt_args$lower <- lower[free_idx]
      opt_args$upper <- upper[free_idx]
    }
    fit_local <- do.call(stats::optim, opt_args)
    par_local <- full_init
    par_local[free_idx] <- opt_to_natural(fit_local$par)
    if ((isTRUE(hard_subcritical) || is.finite(alpha_lo)) && is.finite(beta_eff)) {
      par_local <- .etas_project_subcritical(
        par_local, beta_eff, gap_min = gap_min, eta_max = eta_max,
        alpha_m_lower_bound = alpha_lo
      )
    }
    list(fit = fit_local, par = par_local)
  }

  finalize_candidate <- function(par8, fit_meta = NULL) {
    if ((isTRUE(hard_subcritical) || is.finite(alpha_lo)) && is.finite(beta_eff)) {
      par8 <- .etas_project_subcritical(
        par8, beta_eff, gap_min = gap_min, eta_max = eta_max,
        alpha_m_lower_bound = alpha_lo
      )
    }
    raw_ll <- .etas_raw_loglik(par8, ll_base_args)
    eta <- if (is.finite(beta_eff)) {
      .etas_univ_branching_ratio(par8, beta_eff)
    } else {
      NA_real_
    }
    list(par = par8, value = raw_ll, branching_ratio = eta, fit = fit_meta)
  }

  primary <- run_nm(full_init)
  best <- finalize_candidate(primary$par, primary$fit)
  near_cap <- is.finite(best$branching_ratio) &&
    best$branching_ratio >= (eta_max * near_cap_frac)
  distrust_conv <- !is.null(best$fit$convergence) &&
    (best$fit$convergence != 0L || near_cap)

  if (isTRUE(interior_restart) && A_is_free && is.finite(beta_eff) &&
      (near_cap || distrust_conv ||
       !is.finite(best$value) || best$value <= -1e10)) {
    restart_start <- .etas_scale_to_target_eta(best$par, beta_eff, restart_eta)
    restart_start <- .etas_project_subcritical(
      restart_start, beta_eff, gap_min = gap_min, eta_max = init_eta_cap,
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

  if (isTRUE(polish_productivity) && A_is_free && is.finite(beta_eff)) {
    polished <- .etas_polish_productivity_A(
      best$par, ll_base_args, beta_eff, eta_max = eta_max
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
  fit$branching_ratio <- best$branching_ratio
  fit$soft_branching_barrier <- use_soft_barrier
  fit$log_transform <- isTRUE(log_transform)
  fit$polished_productivity <- isTRUE(best$polished_productivity)
  fit$interior_restarted <- isTRUE(best$interior_restarted)
  fit$near_cap <- is.finite(best$branching_ratio) &&
    best$branching_ratio >= (eta_max * near_cap_frac)
  return(fit)
}


#' Simulate a spatio-temporal ETAS process
#'
#' Generates a realisation of the ETAS model via Poisson-cluster (branching)
#' simulation.  Background events are thinned from a homogeneous (or
#' inhomogeneous) Poisson process; offspring are generated in C++ via
#' \code{sim_etas_children_cpp}.
#'
#' @section Magnitude assignment:
#' Offspring magnitudes are drawn from a Gutenberg-Richter distribution
#' \eqn{m = m_0 + \mathrm{Exp}(\beta_{GR})} when \code{beta_gr > 0},
#' or resampled from \code{mag_pool} otherwise.  Background magnitudes
#' use the same mechanism.
#'
#' @param params  Named list with elements
#'   \code{mu, A, alpha_m, c, p, D, gamma, q}.
#' @param windowT Numeric vector \code{c(start, end)}.
#' @param windowS An \code{owin} object or numeric vector
#'   \code{c(xmin, xmax, ymin, ymax)}.
#' @param m0  Reference magnitude.
#' @param beta_gr  Gutenberg-Richter \eqn{\beta = b \ln 10}.  Set to
#'   \code{NULL} or \eqn{\le 0} to use empirical resampling.
#' @param mag_pool  Numeric vector of magnitudes for resampling.
#' @param background_realization  Optional list with \code{x, y, t, mag}
#'   to use as predetermined background events.
#' @param filtration  Optional data frame of history points that can trigger
#'   offspring but are not themselves part of the new realisation.
#' @param covariate_lookup  Function(x, y) returning spatial covariate
#'   values, or \code{NULL} for homogeneous background.
#' @param t_trunc  Temporal truncation (\code{NULL} = none).
#' @param ...  Additional arguments (ignored).
#' @return A list with elements \code{x, y, t, mag, n, background, W}.
#' @export
sim_etas <- function(params,
                     windowT,
                     windowS,
                     m0,
                     beta_gr = NULL,
                     mag_pool = NULL,
                     background_realization = NULL,
                     filtration = NULL,
                     covariate_lookup = NULL,
                     t_trunc = NULL,
                     ...) {
  if (!inherits(windowS, "owin")) {
    windowS <- owin(xrange = c(windowS[1], windowS[2]),
                    yrange = c(windowS[3], windowS[4]))
  }

  if (!is.list(params)) params <- as.list(params)
  mu      <- params$mu
  A       <- params$A
  alpha_m <- params$alpha_m
  cc      <- params[["c"]]
  p       <- params$p
  D       <- params$D
  gamma_p <- params$gamma
  q       <- params$q

  is_rect <- (windowS$type == "rectangle")
  x_min <- windowS$xrange[1]; x_max <- windowS$xrange[2]
  y_min <- windowS$yrange[1]; y_max <- windowS$yrange[2]
  areaS  <- if (is_rect) (x_max - x_min) * (y_max - y_min) else spatstat.geom::area(windowS)

  empty_result <- list(x = numeric(0), y = numeric(0), t = numeric(0),
                       mag = numeric(0), n = integer(0),
                       background = logical(0), W = numeric(0))

  use_gr <- !is.null(beta_gr) && beta_gr > 0
  pool   <- if (!is.null(mag_pool)) as.numeric(mag_pool) else numeric(0)

  draw_magnitudes <- function(n) {
    if (n == 0) return(numeric(0))
    if (use_gr) {
      m0 + rexp(n, rate = beta_gr)
    } else if (length(pool) > 0) {
      sample(pool, n, replace = TRUE)
    } else {
      rep(m0, n)
    }
  }

  get_covariate <- function(x, y) {
    if (is.null(covariate_lookup)) return(rep(1, length(x)))
    if (is.function(covariate_lookup)) {
      res <- covariate_lookup(x, y)
    } else {
      stop("covariate_lookup must be a function or NULL")
    }
    res[is.na(res)] <- 0
    return(res)
  }

  # --- Background events ---
  if (is.null(background_realization)) {
    if (is.null(covariate_lookup)) {
      n_bg <- rpois(1, mu * (windowT[2] - windowT[1]))
      if (n_bg > 0) {
        if (is_rect) {
          bg_x <- runif(n_bg, x_min, x_max)
          bg_y <- runif(n_bg, y_min, y_max)
          n_ok <- n_bg
        } else {
          # Sample directly on irregular windows to preserve the target
          # immigrant rate mu * |windowT| without bbox-rejection thinning.
          bg_pp <- spatstat.random::runifpoint(n = n_bg, win = windowS)
          bg_x <- bg_pp$x
          bg_y <- bg_pp$y
          n_ok <- bg_pp$n
        }
        if (n_ok > 0) {
          p_x <- bg_x; p_y <- bg_y
          p_t <- sort(runif(n_ok, windowT[1], windowT[2]))
          p_mag <- draw_magnitudes(n_ok)
          p_bg <- rep(TRUE, n_ok)
        } else {
          return(empty_result)
        }
      } else {
        if (A < 1e-10 && is.null(filtration)) return(empty_result)
        p_x <- numeric(0); p_y <- numeric(0)
        p_t <- numeric(0); p_mag <- numeric(0); p_bg <- logical(0)
      }
    } else {
      mu_star <- mu / areaS
      pixel_fun <- function(x, y) {
        (windowT[2] - windowT[1]) * mu_star * get_covariate(x, y)
      }
      bg_pp <- spatstat.random::rpoispp(pixel_fun, win = windowS)
      n_bg <- bg_pp$n
      if (n_bg > 0) {
        p_x <- bg_pp$x; p_y <- bg_pp$y
        p_t <- sort(runif(n_bg, windowT[1], windowT[2]))
        p_mag <- draw_magnitudes(n_bg)
        p_bg <- rep(TRUE, n_bg)
      } else {
        p_x <- numeric(0); p_y <- numeric(0)
        p_t <- numeric(0); p_mag <- numeric(0); p_bg <- logical(0)
      }
    }
  } else {
    if (is.data.frame(background_realization)) {
      p_x <- background_realization$x
      p_y <- background_realization$y
      p_t <- background_realization$t
      p_mag <- background_realization$mag
    } else {
      br <- if (is.list(background_realization)) background_realization
            else as.list(background_realization)
      p_x <- br$x; p_y <- br$y; p_t <- br$t; p_mag <- br$mag
    }
    if (is.null(p_mag)) p_mag <- rep(m0, length(p_x))
    p_bg <- rep(TRUE, length(p_x))
    if (length(p_t) > 1L && is.unsorted(p_t)) {
      ord <- order(p_t)
      p_x <- p_x[ord]; p_y <- p_y[ord]; p_t <- p_t[ord]; p_mag <- p_mag[ord]
    }
  }

  # --- Filtration (history) ---
  if (!is.null(filtration)) {
    if (is.data.frame(filtration)) {
      f_x <- filtration$x; f_y <- filtration$y; f_t <- filtration$t
      f_mag <- if ("mag" %in% names(filtration)) filtration$mag else rep(m0, nrow(filtration))
      n_filt <- nrow(filtration)
    } else if (is.list(filtration)) {
      f_x <- filtration$x; f_y <- filtration$y; f_t <- filtration$t
      f_mag <- if (!is.null(filtration$mag)) filtration$mag else rep(m0, length(f_x))
      n_filt <- length(f_x)
    } else {
      f <- as.data.frame(filtration)
      f_x <- f$x; f_y <- f$y; f_t <- f$t
      f_mag <- if ("mag" %in% names(f)) f$mag else rep(m0, nrow(f))
      n_filt <- nrow(f)
    }
    p_x <- c(p_x, f_x); p_y <- c(p_y, f_y)
    p_t <- c(p_t, f_t); p_mag <- c(p_mag, f_mag)
    p_bg <- c(p_bg, rep(TRUE, n_filt))
  }

  # --- Offspring generation ---
  if (A > 1e-10 && length(p_t) > 0) {
    children <- sim_etas_children_cpp(
      parent_x   = p_x,
      parent_y   = p_y,
      parent_t   = p_t,
      parent_mag = p_mag,
      A          = A,
      alpha_m    = alpha_m,
      cc         = cc,
      p          = p,
      D          = D,
      gamma_par  = gamma_p,
      q          = q,
      m0         = m0,
      beta_gr    = if (use_gr) beta_gr else -1.0,
      t_min      = windowT[1],
      t_max      = windowT[2],
      x_min      = x_min,
      x_max      = x_max,
      y_min      = y_min,
      y_max      = y_max,
      t_trunc    = if (!is.null(t_trunc)) t_trunc else -1.0,
      mag_pool   = pool
    )
    n_ch <- length(children$x)
    if (n_ch > 0) {
      combined_x   <- c(p_x, children$x)
      combined_y   <- c(p_y, children$y)
      combined_t   <- c(p_t, children$t)
      combined_mag <- c(p_mag, children$mag)
      combined_bg  <- c(p_bg, rep(FALSE, n_ch))
    } else {
      combined_x <- p_x; combined_y <- p_y
      combined_t <- p_t; combined_mag <- p_mag; combined_bg <- p_bg
    }
  } else {
    combined_x <- p_x; combined_y <- p_y
    combined_t <- p_t; combined_mag <- p_mag; combined_bg <- p_bg
  }

  n_total <- length(combined_t)
  if (n_total == 0) return(empty_result)

  # --- Temporal filter ---
  valid_t <- combined_t >= windowT[1] & combined_t <= windowT[2]
  if (!all(valid_t)) {
    combined_x   <- combined_x[valid_t];   combined_y   <- combined_y[valid_t]
    combined_t   <- combined_t[valid_t];   combined_mag <- combined_mag[valid_t]
    combined_bg  <- combined_bg[valid_t]
  }

  # --- Spatial filter (non-rectangular windows) ---
  if (!is_rect && length(combined_x) > 0) {
    inside <- inside.owin(combined_x, combined_y, windowS)
    if (!all(inside)) {
      combined_x   <- combined_x[inside];   combined_y   <- combined_y[inside]
      combined_t   <- combined_t[inside];   combined_mag <- combined_mag[inside]
      combined_bg  <- combined_bg[inside]
    }
  }

  n_final <- length(combined_x)
  W_vals <- if (!is.null(covariate_lookup)) {
    get_covariate(combined_x, combined_y)
  } else {
    rep(1.0, n_final)
  }

  list(
    x = combined_x, y = combined_y, t = combined_t,
    mag = combined_mag,
    n = rep(n_final, n_final),
    background = combined_bg, W = W_vals
  )
}


#' @rdname sim_etas
#' @export
sim_etas_fast <- sim_etas


#' Generate inhomogeneous ETAS with different parameters per region
#'
#' Multi-region ETAS analogue of \code{generate_inhomogeneous_hawkes}.
#' Each tile in the partition can have its own ETAS parameter set.
#' Offspring can spill across tile boundaries.
#'
#' @param Omega       Full state space (\code{owin}).
#' @param partition   A \code{spatstat::tess} object.
#' @param time_window Numeric vector \code{c(start, end)}.
#' @param partition_processes  Character vector of process names per tile.
#' @param etas_params Named list of ETAS parameter lists
#'   (e.g.\ \code{list(control = ..., treated = ...)}).
#' @param m0  Reference magnitude.
#' @param beta_gr  Gutenberg-Richter \eqn{\beta} (or \code{NULL} for
#'   resampling from \code{mag_pool}).
#' @param mag_pool  Numeric vector for magnitude resampling.
#' @param state_spaces  Optional pre-computed list of per-process \code{owin}.
#' @param filtration  Optional data frame of history points.
#' @param t_trunc Temporal truncation.
#' @param ...  Additional arguments (e.g.\ \code{covariate_lookup}).
#' @return A \code{data.table} with columns
#'   \code{x, y, t, mag, n, background, W, process, location_process,
#'   tile_index}.
#' @export
generate_inhomogeneous_etas <- function(Omega,
                                        partition,
                                        time_window,
                                        partition_processes,
                                        etas_params,
                                        m0,
                                        beta_gr = NULL,
                                        mag_pool = NULL,
                                        state_spaces = NULL,
                                        filtration = NULL,
                                        t_trunc = NULL,
                                        ...) {
  processes <- unique(partition_processes)

  if (!inherits(Omega, "owin")) Omega <- as.owin(Omega)

  if (is.null(state_spaces)) {
    state_spaces <- lapply(processes, function(pr) {
      idx <- which(partition_processes == pr)
      as.owin(partition[idx])
    })
  }

  dots <- list(...)
  has_covariate <- !is.null(dots$covariate_lookup) && is.function(dots$covariate_lookup)
  filt_by_proc_precomputed <- dots$filt_by_proc

  # --- Background per process ---
  all_bg_x <- list(); all_bg_y <- list(); all_bg_t <- list()
  all_bg_mag <- list(); all_bg_w <- list()
  for (k in seq_along(processes)) {
    pr <- processes[k]
    hp <- etas_params[[pr]]
    if (!is.list(hp)) hp <- as.list(hp)
    if (hp$mu < 1e-10) next
    tmp_hp <- hp; tmp_hp$A <- 0
    ev <- sim_etas(params = tmp_hp, windowT = time_window,
                   windowS = state_spaces[[k]],
                   m0 = m0, beta_gr = beta_gr, mag_pool = mag_pool,
                   background_realization = NULL, filtration = NULL,
                   covariate_lookup = dots$covariate_lookup,
                   t_trunc = t_trunc)
    n_ev <- length(ev$x)
    if (n_ev == 0) next
    all_bg_x[[pr]]   <- ev$x
    all_bg_y[[pr]]   <- ev$y
    all_bg_t[[pr]]   <- ev$t
    all_bg_mag[[pr]] <- ev$mag
    all_bg_w[[pr]]   <- if (length(ev$W) == n_ev) ev$W else numeric(n_ev)
  }

  # --- Split filtration by process ---
  filt_by_proc <- filt_by_proc_precomputed
  if (is.null(filt_by_proc) && !is.null(filtration)) {
    if (is.null(filtration$location_process)) {
      filtration$location_process <-
        partition_processes[tileindex(filtration$x, filtration$y, partition)]
    }
    filt_by_proc <- if (is.data.frame(filtration)) {
      split(filtration, filtration$location_process)
    } else {
      filt_df <- as.data.frame(filtration)
      split(filt_df, filt_df$location_process)
    }
  }

  # --- Full simulation per process ---
  out_x <- list(); out_y <- list(); out_t <- list()
  out_mag <- list(); out_bg <- list(); out_w <- list()
  out_proc <- list(); out_loc <- list(); out_tile <- list()
  idx <- 0L

  for (k in seq_along(processes)) {
    pr <- processes[k]
    f  <- if (!is.null(filt_by_proc) && pr %in% names(filt_by_proc)) filt_by_proc[[pr]] else NULL

    if (is.null(all_bg_x[[pr]])) {
      hp <- etas_params[[pr]]
      if (!is.list(hp)) hp <- as.list(hp)
      if (hp$A < 1e-6 && is.null(f)) next
      bg_realization <- list(x = numeric(0), y = numeric(0),
                             t = numeric(0), mag = numeric(0))
    } else {
      bg_realization <- list(x = all_bg_x[[pr]], y = all_bg_y[[pr]],
                             t = all_bg_t[[pr]], mag = all_bg_mag[[pr]])
    }

    new_events <- sim_etas(
      params = etas_params[[pr]], windowT = time_window, windowS = Omega,
      m0 = m0, beta_gr = beta_gr, mag_pool = mag_pool,
      background_realization = bg_realization, filtration = f,
      covariate_lookup = dots$covariate_lookup, t_trunc = t_trunc
    )
    n_new <- length(new_events$t)
    if (n_new == 0) next

    tile_idx <- as.integer(tileindex(new_events$x, new_events$y, partition))
    loc_proc <- partition_processes[tile_idx]

    w_vals <- if (has_covariate) {
      dots$covariate_lookup(new_events$x, new_events$y)
    } else if (length(new_events$W) == n_new) {
      new_events$W
    } else {
      numeric(n_new)
    }

    idx <- idx + 1L
    out_x[[idx]]    <- new_events$x
    out_y[[idx]]    <- new_events$y
    out_t[[idx]]    <- new_events$t
    out_mag[[idx]]  <- new_events$mag
    out_bg[[idx]]   <- new_events$background
    out_w[[idx]]    <- w_vals
    out_proc[[idx]] <- rep(pr, n_new)
    out_loc[[idx]]  <- loc_proc
    out_tile[[idx]] <- tile_idx
  }

  if (idx == 0L) {
    events <- list(x = numeric(0), y = numeric(0), t = numeric(0),
                   mag = numeric(0), n = integer(0),
                   background = logical(0), W = numeric(0),
                   process = character(0), location_process = character(0),
                   tile_index = integer(0))
    setDT(events)
    return(events)
  }

  events <- list(
    x = unlist(out_x, use.names = FALSE),
    y = unlist(out_y, use.names = FALSE),
    t = unlist(out_t, use.names = FALSE),
    mag = unlist(out_mag, use.names = FALSE),
    n = integer(sum(vapply(out_x, length, 0L))),
    background = unlist(out_bg, use.names = FALSE),
    W = unlist(out_w, use.names = FALSE),
    process = unlist(out_proc, use.names = FALSE),
    location_process = unlist(out_loc, use.names = FALSE),
    tile_index = unlist(out_tile, use.names = FALSE)
  )
  setDT(events)
  events$process[events$t < time_window[1]] <- "control"
  return(events)
}


#' Compensator of spatio-temporal ETAS process
#'
#' Computes the cumulative compensator \eqn{\Lambda(t_i) =
#' \int_0^{t_i}\!\int_S \lambda(s,u)\,\mathrm{d}u\,\mathrm{d}s}
#' evaluated at each event time, for use in residual analysis and
#' goodness-of-fit testing (random time-change theorem).
#'
#' @inheritParams loglik_etas
#' @param m0  Reference magnitude.
#' @return Numeric vector of compensator values at each event time.
#' @export
compensator_etas <- function(params,
                             realiz,
                             windowT,
                             windowS,
                             m0 = NULL,
                             zero_background_region = NULL,
                             ...) {
  realiz <- realiz[order(realiz$t), ]
  realiz$t <- realiz$t - windowT[1]
  realiz <- realiz[realiz$t >= 0, ]

  windowS <- as.owin(windowS)

  if (is.list(params)) {
    mu <- params$mu; A <- params$A; alpha_m <- params$alpha_m
    cc <- params[["c"]]; p <- params$p
  } else {
    mu <- params[1]; A <- params[2]; alpha_m <- params[3]
    cc <- params[4]; p <- params[5]
  }

  if (is.null(m0)) m0 <- min(realiz$mag)

  if (is.null(zero_background_region)) {
    adjust_factor <- 1
  } else {
    area_zero <- spatstat.geom::area(as.owin(zero_background_region))
    adjust_factor <- (spatstat.geom::area(windowS) - area_zero) /
                      spatstat.geom::area(windowS)
  }

  kappa_vec <- A * exp(alpha_m * (realiz$mag - m0))

  n <- nrow(realiz)
  incremental <- numeric(n)
  for (i in seq_len(n)) {
    bg_part <- realiz$t[i] * adjust_factor * mu
    trig_part <- 0
    if (i >= 2) {
      for (j in 1:(i - 1)) {
        horizon <- realiz$t[i] - realiz$t[j]
        trig_part <- trig_part +
          kappa_vec[j] * (1 - (1 + horizon / cc)^(-(p - 1)))
      }
    }
    incremental[i] <- bg_part + trig_part
  }
  return(incremental)
}


#' KS test p-value for ETAS goodness-of-fit
#'
#' Applies the random time-change theorem: if the model is correct, the
#' compensator increments are i.i.d.\ Exp(1), so
#' \eqn{1 - \exp(-\Delta\Lambda_i) \sim \mathrm{Uniform}(0,1)}.
#'
#' @inheritParams compensator_etas
#' @param etas_par  Named list of ETAS parameters.
#' @return Numeric p-value from the Kolmogorov-Smirnov test.
#' @export
ks_test_pval_etas <- function(realiz, windowT, windowS, etas_par, m0 = NULL,
                              zero_background_region = NULL, ...) {
  compensators <- compensator_etas(
    params = etas_par,
    realiz = realiz,
    windowT = windowT,
    windowS = windowS,
    m0 = m0,
    zero_background_region = zero_background_region,
    ...
  )
  compensator_incs <- diff(compensators)
  test_dist <- 1 - exp(-compensator_incs)
  test <- ks.test(test_dist, "punif")
  return(test$p.value)
}
