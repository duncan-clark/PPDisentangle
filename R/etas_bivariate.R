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
#' @param t_trunc Temporal truncation (NULL = none).
#' @param ... Ignored.
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
                                  beta_gr = NULL,
                                  stability_barrier_start = 0.95,
                                  stability_barrier_weight = 100,
                                  stability_barrier_power = 2,
                                  t_trunc = NULL,
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

  if (min(mu_0, mu_1, A_00, A_11, cc, D) < 0 ||
      A_01 < 0 || A_10 < 0 || p <= 1 || q <= 1 || gamma_p < 0) {
    return(-1e15)
  }

  if (is.unsorted(realiz$t)) realiz <- realiz[order(realiz$t), ]
  t_idx <- realiz$t >= windowT[1] & realiz$t <= windowT[2]
  if (!all(t_idx)) realiz <- realiz[t_idx, ]
  n <- nrow(realiz)
  if (n == 0) return(-1e15)

  if (is.null(m0)) m0 <- min(realiz$mag)

  # Build process_id vector (0 = control, 1 = treated)
  dots <- list(...)
  precomp <- dots$precomp
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
      W_0[in_treated] <- 0
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

    if (areaS_0 <= 0) areaS_0 <- 1
    if (areaS_1 <= 0) areaS_1 <- 1
  }

  # Optional inhomogeneous background covariate:
  # W scales the baseline intensity for both processes after region masking.
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

  # Optional policy mask: force treated-process background to zero before
  # treatment (or any user-specified cutoff), while keeping control unchanged.
  if (!is.null(treated_background_zero_before)) {
    W_1[realiz$t < as.numeric(treated_background_zero_before)] <- 0
  }

  tval <- windowT[2] - windowT[1]

  loglik <- etas_bivariate_loglik_cpp(
    t         = realiz$t - windowT[1],
    x         = realiz$x,
    y         = realiz$y,
    mag       = realiz$mag,
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
    t_trunc = if (!is.null(t_trunc)) t_trunc else -1.0
  )
  # Smooth stability barrier on the bivariate ETAS offspring matrix spectral
  # radius; activates only when rho exceeds `stability_barrier_start`.
  if (is.finite(stability_barrier_weight) && stability_barrier_weight > 0) {
    beta_eff <- suppressWarnings(as.numeric(beta_gr))
    if (!is.finite(beta_eff) || is.na(beta_eff) || beta_eff <= 0) {
      mag_delta <- as.numeric(realiz$mag) - as.numeric(m0)
      mag_delta <- mag_delta[is.finite(mag_delta) & mag_delta > 0]
      beta_eff <- if (length(mag_delta) > 0L) 1 / mean(mag_delta) else 1
    }
    eta_comp <- function(Ai, alphai) {
      Ai <- suppressWarnings(as.numeric(Ai))
      alphai <- suppressWarnings(as.numeric(alphai))
      if (!is.finite(Ai) || !is.finite(alphai) || !is.finite(beta_eff) || beta_eff <= 0) return(Inf)
      gap <- beta_eff - alphai
      if (!is.finite(gap) || gap <= 1e-8) return(Inf)
      Ai * beta_eff / gap
    }
    M <- matrix(
      c(
        eta_comp(A_00, alpha_m_00), eta_comp(A_01, alpha_m_01),
        eta_comp(A_10, alpha_m_10), eta_comp(A_11, alpha_m_11)
      ),
      nrow = 2, byrow = TRUE
    )
    rho <- if (all(is.finite(M))) {
      max(Re(eigen(M, only.values = TRUE)$values))
    } else {
      Inf
    }
    if (!is.finite(rho)) return(-1e15)
    barrier_start <- suppressWarnings(as.numeric(stability_barrier_start))
    if (!is.finite(barrier_start) || is.na(barrier_start)) barrier_start <- 0.95
    barrier_power <- suppressWarnings(as.numeric(stability_barrier_power))
    if (!is.finite(barrier_power) || is.na(barrier_power) || barrier_power <= 0) barrier_power <- 2
    excess <- max(0, rho - barrier_start)
    if (excess > 0) {
      loglik <- loglik - stability_barrier_weight * (excess ^ barrier_power)
    }
  }
  loglik
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
#' @param maxit Maximum optim iterations.
#' @param fixed_params Named list of parameters to hold fixed.
#' @param symmetric If TRUE, constrain A_01=A_10 and alpha_m_01=alpha_m_10.
#' @param trace Trace level for optim.
#' @param t_trunc Temporal truncation.
#' @param ... Passed to loglik_etas_bivariate.
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
                               maxit = 5000,
                               fixed_params = NULL,
                               symmetric = FALSE,
                               trace = 0,
                               t_trunc = NULL,
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
  free_init <- full_init[free_idx]

  realiz <- realiz[order(realiz$t), ]
  t_idx <- realiz$t >= windowT[1] & realiz$t <= windowT[2]
  if (!all(t_idx)) realiz <- realiz[t_idx, ]
  n <- nrow(realiz)

  W_0 <- rep(1.0, n); W_1 <- rep(1.0, n)
  total_area <- spatstat.geom::area(as.owin(windowS))
  if (!is.null(treated_state_space)) {
    areaS_1 <- spatstat.geom::area(as.owin(treated_state_space))
    W_0[inside.owin(realiz$x, realiz$y, treated_state_space)] <- 0
  } else { areaS_1 <- total_area / 2 }
  if (!is.null(control_state_space)) {
    areaS_0 <- spatstat.geom::area(as.owin(control_state_space))
    W_1[inside.owin(realiz$x, realiz$y, control_state_space)] <- 0
  } else { areaS_0 <- total_area - areaS_1 }
  if (areaS_0 <= 0) areaS_0 <- 1
  if (areaS_1 <= 0) areaS_1 <- 1

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
  precomp <- list(W_0 = W_0, W_1 = W_1, areaS_0 = areaS_0, areaS_1 = areaS_1)

  profile_fn <- function(free_par) {
    par15 <- full_init
    par15[free_idx] <- free_par

    if (symmetric) {
      par15["A_10"] <- par15["A_01"]
      par15["alpha_m_10"] <- par15["alpha_m_01"]
    }

    ll_args <- list(
      params = par15, realiz = realiz,
      windowT = windowT, windowS = windowS,
      m0 = m0,
      control_state_space = control_state_space,
      treated_state_space = treated_state_space,
      background_rate_var = background_rate_var,
      treated_background_zero_before = treated_background_zero_before,
      t_trunc = t_trunc,
      precomp = precomp
    )
    # Only forward optional barrier controls when explicitly set so
    # defaults in loglik_etas_bivariate remain intact.
    if (!is.null(dots$beta_gr)) ll_args$beta_gr <- dots$beta_gr
    if (!is.null(dots$stability_barrier_start)) ll_args$stability_barrier_start <- dots$stability_barrier_start
    if (!is.null(dots$stability_barrier_weight)) ll_args$stability_barrier_weight <- dots$stability_barrier_weight
    if (!is.null(dots$stability_barrier_power)) ll_args$stability_barrier_power <- dots$stability_barrier_power
    do.call(loglik_etas_bivariate, ll_args)
  }

  fit <- stats::optim(
    par     = free_init,
    fn      = profile_fn,
    method  = "Nelder-Mead",
    control = list(fnscale = -1, trace = trace, maxit = maxit)
  )

  par15 <- full_init
  par15[free_idx] <- fit$par
  if (symmetric) {
    par15["A_10"] <- par15["A_01"]
    par15["alpha_m_10"] <- par15["alpha_m_01"]
  }
  fit$par <- par15
  fit$m0 <- m0
  return(fit)
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

  # Background for each process in its own state space
  gen_bg <- function(mu, ss, proc_id) {
    if (mu < 1e-10 || is.null(ss)) {
      return(list(x = numeric(0), y = numeric(0), t = numeric(0),
                  mag = numeric(0), proc = integer(0)))
    }
    n_bg <- rpois(1, mu * (windowT[2] - windowT[1]))
    if (n_bg == 0) {
      return(list(x = numeric(0), y = numeric(0), t = numeric(0),
                  mag = numeric(0), proc = integer(0)))
    }
    bb <- ss$xrange
    bby <- ss$yrange
    bx <- runif(n_bg * 3, bb[1], bb[2])
    by <- runif(n_bg * 3, bby[1], bby[2])
    ok <- inside.owin(bx, by, ss)
    bx <- bx[ok]; by <- by[ok]
    if (length(bx) > n_bg) { bx <- bx[1:n_bg]; by <- by[1:n_bg] }
    n_ok <- length(bx)
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

  bg0 <- gen_bg(mu_0, ctrl_ss, 0L)
  bg1 <- gen_bg(mu_1, treat_ss, 1L)

  p_x   <- c(bg0$x, bg1$x)
  p_y   <- c(bg0$y, bg1$y)
  p_t   <- c(bg0$t, bg1$t)
  p_mag <- c(bg0$mag, bg1$mag)
  p_proc <- c(bg0$proc, bg1$proc)
  p_bg  <- rep(TRUE, length(p_x))

  # Add filtration
  if (!is.null(filtration)) {
    if (is.data.frame(filtration)) {
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

  any_excitation <- max(pv[["A_00"]], pv[["A_11"]], pv[["A_01"]], pv[["A_10"]]) > 1e-10

  if (any_excitation && length(p_t) > 0) {
    children <- sim_etas_bivariate_children_cpp(
      parent_x = p_x, parent_y = p_y,
      parent_t = p_t, parent_mag = p_mag,
      parent_process = as.integer(p_proc),
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
      mag_pool = pool
    )

    n_ch <- length(children$x)
    if (n_ch > 0) {
      p_x <- c(p_x, children$x);     p_y <- c(p_y, children$y)
      p_t <- c(p_t, children$t);     p_mag <- c(p_mag, children$mag)
      p_proc <- c(p_proc, children$process_id)
      p_bg <- c(p_bg, rep(FALSE, n_ch))
    }
  }

  # Filter to window
  valid <- p_t >= windowT[1] & p_t <= windowT[2]
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

  result <- sim_etas_bivariate(
    params = etas_bivariate_params,
    windowT = time_window,
    windowS = Omega,
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

  tile_idx <- as.integer(tileindex(result$x, result$y, partition))
  result$tile_index <- tile_idx
  result$location_process <- partition_processes[pmin(pmax(tile_idx, 1), length(partition_processes))]
  result$location_process[is.na(tile_idx)] <- "control"
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
