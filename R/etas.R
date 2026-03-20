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
                        precomp = NULL,
                        t_trunc = NULL,
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

  # --- Sorting ---
  if (is.unsorted(realiz$t)) realiz <- realiz[order(realiz$t), ]

  # --- Time window filter ---
  t_idx <- realiz$t >= windowT[1] & realiz$t <= windowT[2]
  if (!all(t_idx)) realiz <- realiz[t_idx, ]
  n <- nrow(realiz)
  if (n == 0) return(-1e15)

  # --- Parameter bounds ---
  if (min(mu, A, cc, D) < 0 || p <= 1 || q <= 1 || gamma_p < 0) return(-1e15)

  # --- Reference magnitude ---
  if (is.null(m0)) m0 <- min(realiz$mag)

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
                     log_transform = FALSE,
                     init_decluster = FALSE,
                     ...) {
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

  fixed_idx <- if (!is.null(fixed_params)) match(names(fixed_params), all_names) else integer(0)
  free_idx  <- setdiff(seq_along(all_names), fixed_idx)
  free_init <- full_init[free_idx]
  free_names <- all_names[free_idx]

  # Indices (within the free vector) of params that must be positive
  log_idx <- if (log_transform) which(free_names %in% c("mu", "A", "c", "D")) else integer(0)

  realiz <- realiz[order(realiz$t), ]

  # Transform initial values to optimisation scale
  opt_init <- free_init
  if (length(log_idx) > 0) opt_init[log_idx] <- log(free_init[log_idx])

  profile_fn <- function(opt_par, ...) {
    free_par <- opt_par
    if (length(log_idx) > 0) free_par[log_idx] <- exp(opt_par[log_idx])
    par8 <- full_init
    par8[free_idx] <- free_par
    loglik_etas(par8, ...)
  }

  opt_args <- c(
    list(
      par     = opt_init,
      fn      = profile_fn,
      method  = method,
      control = list(fnscale = -1, trace = trace, maxit = maxit),
      realiz  = realiz,
      windowT = windowT,
      windowS = windowS,
      m0      = m0,
      t_trunc = t_trunc
    ),
    list(...)
  )
  if (method == "L-BFGS-B" && !is.null(lower) && !is.null(upper)) {
    opt_args$lower <- lower[free_idx]
    opt_args$upper <- upper[free_idx]
  }

  fit <- do.call(stats::optim, opt_args)

  # Back-transform to original scale
  free_result <- fit$par
  if (length(log_idx) > 0) free_result[log_idx] <- exp(fit$par[log_idx])

  par8 <- full_init
  par8[free_idx] <- free_result
  fit$par <- par8
  fit$m0 <- m0
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
