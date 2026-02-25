#' @import spatstat
#' @import spatstat.geom
#' @import spatstat.random
#' @importFrom spatstat.geom owin as.owin inside.owin area.owin tileindex
#'   tilenames tiles as.polygonal quadrats
#' @importFrom data.table data.table rbindlist setDT
#' @importFrom Rcpp sourceCpp
#' @importFrom stats optim pnorm runif rpois rexp sd ks.test rbinom rmultinom
#' @importFrom dplyr `%>%` filter arrange pull group_by summarize n
#' @importFrom rlang .data
#' @useDynLib PPDisentangle, .registration = TRUE
NULL

#' Log-likelihood for a spatio-temporal Hawkes process (R version)
#'
#' @param params Numeric vector c(mu, alpha, beta, K)
#' @param realiz Data frame with columns x, y, t
#' @param windowT Numeric vector c(start, end) for the time window
#' @param windowS An owin object or convertible to one
#' @param background_rate_var Optional column name in realiz for inhomogeneous background
#' @param dists Optional precomputed distance list (space_dist, time_dist)
#' @param zero_background_region Optional owin where background intensity is zero
#' @param density_approx Logical; use density approximation for integral
#' @param numeric_integral Logical; use numeric integration
#' @param cl Optional parallel cluster
#' @param optimized Logical; use vectorized computation
#' @param ... Additional arguments
#' @return Scalar log-likelihood value
#' @export
loglik_hawk <- function(params,
                        realiz,
                        windowT,
                        windowS,
                        background_rate_var = NULL,
                        dists = NULL,
                        zero_background_region = NULL,
                        density_approx = TRUE,
                        numeric_integral = FALSE,
                        cl = NULL,
                        optimized = TRUE,
                        ...) {
  if (!inherits(windowS, "owin")) {
    windowS <- as.owin(windowS)
  }
  windowS <- as.owin(windowS)
  tval <- windowT[2] - windowT[1]
  max_t <- max(realiz$t)
  if (windowT[1] > min(realiz$t) | windowT[2] < max(realiz$t)) {
    warning("realization has points outside time window - dropping them")
    realiz <- realiz[realiz$t >= windowT[1] & realiz$t <= windowT[2], ]
  }
  mu <- params[1]
  alpha <- params[2]
  beta <- params[3]
  K <- params[4]

  if (min(mu, K, alpha, beta) < 0) return(-999999)
  if (K > .99999) return(-999999)

  if (is.null(dists)) {
    realiz <- realiz[order(realiz$t), ]
    x_diff <- outer(realiz$x, realiz$x, "-")
    y_diff <- outer(realiz$y, realiz$y, "-")
    space_dist <- x_diff^2 + y_diff^2
    time_dist <- outer(realiz$t, realiz$t, "-")
  } else {
    space_dist <- dists$space_dist
    time_dist <- dists$time_dist
  }
  exp_decay_mat <- exp(-beta * time_dist - alpha * space_dist)
  exp_decay_mat[upper.tri(exp_decay_mat)] <- 0

  if (!is.null(zero_background_region)) {
    area_zero_background <- spatstat.geom::area(owin(zero_background_region))
    adjust_factor <- (spatstat.geom::area.owin(windowS) - area_zero_background) / spatstat.geom::area.owin(windowS)
    in_zero_background <- 1 * inside.owin(realiz[, c("x", "y")], w = zero_background_region)
  } else {
    area_zero_background <- 0
    adjust_factor <- 1
    in_zero_background <- rep(0, dim(realiz)[1])
  }
  if (density_approx) {
    mu_star <- mu / (spatstat.geom::area(windowS))
    intlam <- adjust_factor * mu * tval + K * dim(realiz)[1]
    const <- K * (alpha / pi) * beta
  } else {
    mu_star <- mu / (spatstat.geom::area(windowS))
    intlam <- adjust_factor * mu * tval
    const <- K * (alpha / pi) * beta
    if (numeric_integral) {
      const <- K * (alpha / pi) * beta
      int_func <- function(x, realiz) {
        t <- x[1]
        y <- x[3]
        x <- x[2]
        x_diff <- x - realiz$x
        y_diff <- y - realiz$y
        t_diff <- t - realiz$t
        space_dist <- x_diff^2 + y_diff^2
        time_dist <- t_diff
        result <- sum(exp(-beta * t_diff - alpha * space_dist))
        if ((result == 0)) {
          return(1e-10)
        }
        return(K * result)
      }
      func <- function(i) {
        max_t <- realiz$t[i]
        min_t <- realiz$t[i - 1]
        cubature::hcubature(int_func,
          realiz = realiz[realiz$t == min_t, ],
          lowerLimit = c(min_t, windowS$xrange[1], windowS$yrange[1]),
          upperLimit = c(max_t, windowS$xrange[2], windowS$yrange[2])
        )$integral
      }
      if (!is.null(cl)) {
        parallel::clusterExport(cl, c(
          "func", "int_func", "hcubature",
          "realiz", "x_diff", "y_diff", "t_diff",
          "space_dist", "time_dist"
        ), envir = environment())
        pieces <- parallel::parSapply(cl = cl, X = 2:dim(realiz)[1], FUN = func)
      } else {
        pieces <- sapply(2:dim(realiz)[1], func)
      }
    } else {
      func_3 <- function(x, y, t) {
        t_comp <- (1 - exp(-beta * (max_t - t)))
        s_1 <- (pnorm(windowS$xrange[2], mean = x, sd = 1 / sqrt(2 * alpha)) -
          pnorm(0, mean = x, sd = 1 / sqrt(2 * alpha)))
        s_2 <- (pnorm(windowS$yrange[2], mean = y, sd = 1 / sqrt(2 * alpha)) -
          pnorm(0, mean = y, sd = 1 / sqrt(2 * alpha)))
        return(t_comp * s_1 * s_2)
      }
      if (optimized) {
        pieces <- func_3(realiz$x, realiz$y, realiz$t)
      } else {
        func_2 <- function(i) {
          t_comp <- (1 - exp(-beta * (max_t - realiz$t[i])))
          s_1 <- (pnorm(windowS$xrange[2], mean = realiz$x[i], sd = 1 / sqrt(2 * alpha)) -
            pnorm(0, mean = realiz$x[i], sd = 1 / sqrt(2 * alpha)))
          s_2 <- (pnorm(windowS$yrange[2], mean = realiz$y[i], sd = 1 / sqrt(2 * alpha)) -
            pnorm(0, mean = realiz$y[i], sd = 1 / sqrt(2 * alpha)))
          return(t_comp * s_1 * s_2)
        }
        pieces <- sapply(1:dim(realiz)[1], func_2)
      }
    }
    if (K != 0 & (dim(realiz)[1] != 1)) {
      intlam <- intlam + K * sum(pieces)
    }
  }

  if (!is.null(background_rate_var)) {
    mu_star <- mu_star * realiz[[background_rate_var]]
  } else {
    mu_star <- rep(mu_star, dim(realiz)[1])
  }
  sum_log <- log(mu_star[1])
  lamjs <- c()
  if (optimized) {
    gij_vec <- rowSums(
      lower.tri(exp_decay_mat, diag = FALSE) * exp_decay_mat
    )
    lamjs <- (1 - in_zero_background) * mu_star + const * gij_vec
    lamjs <- lamjs[2:length(lamjs)]
  } else {
    if (nrow(realiz) >= 2) {
      for (j in 2:nrow(realiz)) {
        gij <- sum(exp_decay_mat[j, 1:(j - 1)])
        lamjs[j - 1] <- (1 - in_zero_background[j]) * mu_star[j] + const * gij
      }
    } else {
      sum_log <- 0
    }
  }
  if (any(lamjs == 0)) {
    warning("Some points have 0 intensity, setting to min value")
    lamjs[lamjs == 0] <- min(lamjs[lamjs > 0])
  }
  if (any(is.na(lamjs)) || any(lamjs < 0)) {
    return(-999999)
  } else {
    sum_log <- sum_log + sum(log(lamjs))
  }
  loglik <- sum_log - intlam
  if (loglik == -Inf) {
    return(-999999)
  }
  if (K == 0) {
    n <- sum(1 - in_zero_background)
    stirling_approx <- n * log(n) - n + 0.5 * log(2 * pi * n)
    loglik <- loglik + n * log(tval * spatstat.geom::area(windowS) * adjust_factor) - stirling_approx
  }
  return(loglik)
}

#' Fast log-likelihood for a spatio-temporal Hawkes process (C++ accelerated)
#'
#' @param params Numeric vector c(mu, alpha, beta, K) or named list
#' @param realiz Data frame with columns x, y, t (and optionally W)
#' @param windowT Numeric vector c(start, end) for the time window
#' @param windowS An owin object or convertible to one
#' @param zero_background_region Optional owin where background intensity is zero
#' @param background_rate_var Column name for inhomogeneous background (default "W")
#' @param precomp Optional list from \code{precompute_loglik_args} to skip
#'   redundant area/inside.owin calculations when called in a tight loop.
#' @param ... Additional arguments
#' @return Scalar log-likelihood value
#' @export
loglik_hawk_fast <- function(params,
                             realiz,
                             windowT,
                             windowS,
                             zero_background_region = NULL,
                             background_rate_var = "W",
                             precomp = NULL,
                             ...) {
  if (is.list(params)) {
    mu <- params$mu; alpha <- params$alpha; beta <- params$beta; K <- params$K
  } else {
    mu <- params[1]; alpha <- params[2]; beta <- params[3]; K <- params[4]
  }

  if (is.unsorted(realiz$t)) realiz <- realiz[order(realiz$t), ]

  t_idx <- realiz$t >= windowT[1] & realiz$t <= windowT[2]
  if (!all(t_idx)) realiz <- realiz[t_idx, ]
  n <- nrow(realiz)
  if (n == 0) return(-Inf)
  if (min(mu, K, alpha, beta) < 0 || K >= 1) return(-Inf)

  W_vec <- if (!is.null(background_rate_var) && background_rate_var %in% names(realiz)) {
    realiz[[background_rate_var]]
  } else {
    rep(1.0, n)
  }

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
      if (active_area <= 0) warning("Zero background region covers entire window!")
      in_zero_bg <- inside.owin(realiz[, c("x", "y")], w = zero_background_region)
      W_vec[in_zero_bg] <- 0
    } else {
      active_area <- total_area
    }
  }

  tval <- windowT[2] - windowT[1]
  intlam <- (mu * tval) + (K * n)

  sum_log <- hawkes_loglik_inhom_cpp(
    t = realiz$t,
    x = realiz$x,
    y = realiz$y,
    W_val = W_vec,
    mu = mu,
    alpha = alpha,
    beta = beta,
    K = K,
    areaS = active_area,
    t_max = windowT[2]
  )

  loglik <- sum_log - intlam

  if (K == 0) {
    stirling_approx <- n * log(n) - n + 0.5 * log(2 * pi * n)
    loglik <- loglik + n * log(tval * active_area) - stirling_approx
  }

  return(loglik)
}

#' Pre-compute reusable values for repeated loglik_hawk_fast calls
#'
#' Computes owin areas and inside.owin membership once so they can be
#' reused across many likelihood evaluations on the same point set.
#'
#' @param all_points Data frame with x, y columns for all candidate points
#' @param windowS Full observation window (owin)
#' @param zero_background_region Optional owin where background is zero
#' @return A list with \code{active_area} and optionally
#'   \code{in_zero_bg_all} (logical vector for all_points)
#' @export
precompute_loglik_args <- function(all_points, windowS, zero_background_region = NULL) {
  if (!inherits(windowS, "owin")) windowS <- as.owin(windowS)
  total_area <- spatstat.geom::area(windowS)

  if (!is.null(zero_background_region)) {
    if (!inherits(zero_background_region, "owin")) {
      zero_background_region <- as.owin(zero_background_region)
    }
    zero_area <- spatstat.geom::area(zero_background_region)
    active_area <- total_area - zero_area
    in_zero_bg_all <- inside.owin(all_points[, c("x", "y")],
                                  w = zero_background_region)
  } else {
    active_area <- total_area
    in_zero_bg_all <- NULL
  }

  list(active_area = active_area, in_zero_bg_all = in_zero_bg_all)
}

#' Compensator of spatio-temporal Hawkes process (for RCT theorem)
#'
#' @param params Numeric vector c(mu, alpha, beta, K)
#' @param realiz Data frame with columns x, y, t
#' @param windowT Numeric vector c(start, end)
#' @param windowS An owin object
#' @param zero_background_region Optional owin where background is zero
#' @param optimized Logical; use vectorized computation
#' @return Numeric vector of compensator increments
#' @export
compensator_hawkes <- function(params,
                               realiz,
                               windowT,
                               windowS,
                               zero_background_region = NULL,
                               optimized = TRUE) {
  realiz <- realiz[order(realiz$t), ]
  realiz$t <- realiz$t - windowT[1]
  realiz <- realiz[realiz$t >= 0, ]

  windowS <- as.owin(windowS)
  tval <- windowT[2] - windowT[1]
  max_t <- max(realiz$t)

  mu <- params[1]
  alpha <- params[2]
  beta <- params[3]
  K <- params[4]

  if (is.null(zero_background_region)) {
    adjust_factor <- 1
    in_zero_background <- rep(0, dim(realiz)[1])
  } else {
    area_zero_background <- spatstat.geom::area(owin(zero_background_region))
    adjust_factor <- (spatstat.geom::area.owin(windowS) - area_zero_background) / spatstat.geom::area.owin(windowS)
    in_zero_background <- 1 * inside.owin(realiz[, c("x", "y")], w = zero_background_region)
  }
  mu_star <- mu / (spatstat.geom::area(windowS))
  intlam <- adjust_factor * mu

  func_2 <- function(i, max_t) {
    t_comp <- (1 - exp(-beta * (max_t - realiz$t[i])))
    s_1 <- (pnorm(windowS$xrange[2], mean = realiz$x[i], sd = 1 / sqrt(2 * alpha)) -
      pnorm(0, mean = realiz$x[i], sd = 1 / sqrt(2 * alpha)))
    s_2 <- (pnorm(windowS$yrange[2], mean = realiz$y[i], sd = 1 / sqrt(2 * alpha)) -
      pnorm(0, mean = realiz$y[i], sd = 1 / sqrt(2 * alpha)))
    return(t_comp * s_1 * s_2)
  }

  if (optimized) {
    t_comps <- outer(
      X = realiz$t, Y = realiz$t,
      FUN = function(x, y) {
        (x >= y) * (1 - exp(-beta * (x - y)))
      }
    )
    s_1 <- (pnorm(windowS$xrange[2], mean = realiz$x, sd = 1 / sqrt(2 * alpha)) -
      pnorm(0, mean = realiz$x, sd = 1 / sqrt(2 * alpha)))
    s_2 <- (pnorm(windowS$yrange[2], mean = realiz$y, sd = 1 / sqrt(2 * alpha)) -
      pnorm(0, mean = realiz$y, sd = 1 / sqrt(2 * alpha)))
    t_comps <- t_comps * matrix(s_1, byrow = TRUE, nrow = dim(t_comps)[1], ncol = dim(t_comps)[2])
    t_comps <- t_comps * matrix(s_2, byrow = TRUE, nrow = dim(t_comps)[1], ncol = dim(t_comps)[2])
    pieces <- rowSums(t_comps)
    incremental <- realiz$t * intlam + K * pieces
  } else {
    pieces <- lapply(1:dim(realiz)[1], function(i) {
      sapply(1:i, function(j) { func_2(j, realiz$t[i]) })
    })
    incremental <- sapply(1:length(pieces), function(i) {
      realiz$t[i] * (intlam) + K * sum(pieces[[i]])
    })
  }
  return(incremental)
}

#' KS test p-value for Hawkes process goodness-of-fit
#'
#' @param realiz Data frame with columns x, y, t
#' @param windowT Numeric vector c(start, end)
#' @param windowS An owin object
#' @param hawkes_par List with elements mu, alpha, beta, K
#' @param zero_background_region An owin object
#' @return Numeric p-value from the KS test
#' @export
ks_test_pval <- function(realiz, windowT, windowS, hawkes_par, zero_background_region) {
  compensators <- compensator_hawkes(
    params = unlist(hawkes_par),
    realiz = realiz,
    windowT = windowT,
    windowS = windowS,
    zero_background_region = zero_background_region,
    optimized = TRUE
  )
  compensator_incs <- diff(compensators)
  test_dist <- 1 - exp(-compensator_incs)
  test <- ks.test(test_dist, "punif")
  return(test$p.value)
}

#' Fit a spatio-temporal Hawkes process via MLE
#'
#' @param params_init Initial parameter values (list or numeric vector)
#' @param realiz Data frame with columns x, y, t
#' @param windowT Numeric vector c(start, end)
#' @param windowS An owin object
#' @param trace Trace level for optim (0 = silent)
#' @param maxit Maximum iterations for optim
#' @param poisson_flag If TRUE, return Poisson MLE directly
#' @param zero_background_region Optional owin where background is zero
#' @param ... Additional arguments passed to loglik_hawk
#' @return An optim result list with element par
#' @export
fit_hawkes <- function(params_init,
                       realiz,
                       windowT,
                       windowS,
                       trace = 0,
                       maxit,
                       poisson_flag = FALSE,
                       zero_background_region = NULL,
                       ...) {
  if (inherits(params_init, "list")) { params_init <- unlist(params_init) }
  if (poisson_flag) {
    if (is.null(zero_background_region)) {
      return(list(par = list(mu = dim(realiz)[1] / windowT[2], alpha = 1, beta = 1, K = 0), converged = TRUE))
    } else {
      win_area <- spatstat.geom::area(windowS)
      non_zero_area <- spatstat.geom::area.owin(zero_background_region)
      adjust_factor <- (win_area - non_zero_area) / win_area
      return(list(par = list(mu = dim(realiz)[1] / (adjust_factor * windowT[2]), alpha = 1, beta = 1, K = 0), converged = TRUE))
    }
  } else {
    realiz <- realiz[order(realiz$t), ]
    x_diff <- outer(realiz$x, realiz$x, "-")
    y_diff <- outer(realiz$y, realiz$y, "-")
    space_dist <- x_diff^2 + y_diff^2
    time_dist <- outer(realiz$t, realiz$t, "-")
    dists <- list(space_dist = space_dist, time_dist = time_dist)

    fit <- optim(
      par = params_init,
      fn = loglik_hawk,
      method = "Nelder-Mead",
      control = list(fnscale = -1, trace = trace, maxit = maxit),
      realiz = realiz,
      windowT = windowT,
      windowS = windowS,
      dists = dists,
      zero_background_region = zero_background_region,
      ...
    )
  }
  return(fit)
}

#' Simulate a spatio-temporal Hawkes process
#'
#' @param params List with elements mu, alpha, beta, K
#' @param windowT Numeric vector c(start, end)
#' @param windowS An owin object or numeric vector c(xmin, xmax, ymin, ymax)
#' @param background_realization Optional list to use as background events
#' @param filtration Optional data frame of history points that can trigger
#' @param optimized Logical; use vectorized aftershock generation
#' @param covariate_lookup Function(x,y) or Raster for inhomogeneous background
#' @param background_rate_var Column name for covariate values
#' @return List with elements x, y, t, n, background, W
#' @export
sim_hawkes <- function(params,
                       windowT,
                       windowS,
                       background_realization = NULL,
                       filtration = NULL,
                       optimized = TRUE,
                       covariate_lookup = NULL,
                       background_rate_var = NULL) {
  if (!inherits(windowS, "owin")) {
    windowS <- owin(
      xrange = c(windowS[1], windowS[2]),
      yrange = c(windowS[3], windowS[4])
    )
  }
  mu <- params$mu
  alpha <- params$alpha
  beta <- params$beta
  K <- params$K
  events <- list()
  events$n <- 0
  events$t <- numeric()

  areaS <- spatstat.geom::area(windowS)
  mu_star <- mu / areaS

  get_covariate <- function(x, y) {
    if (is.null(covariate_lookup)) return(rep(0, length(x)))
    if (is.function(covariate_lookup)) {
      res <- covariate_lookup(x, y)
    } else if (inherits(covariate_lookup, "Raster")) {
      res <- raster::extract(covariate_lookup, cbind(x, y))
    } else {
      stop("covariate_lookup must be a function or a Raster* object")
    }
    res[is.na(res)] <- 0
    return(res)
  }

  if (is.null(background_realization)) {
    if (is.null(covariate_lookup)) {
      n_bg <- rpois(1, mu * (windowT[2] - windowT[1]))
      events$n <- n_bg
      bg_pp <- runifpoint(n = n_bg, win = windowS)
    } else {
      pixel_fun <- function(x, y) { (windowT[2] - windowT[1]) * mu_star * get_covariate(x, y) }
      bg_pp <- spatstat.random::rpoispp(pixel_fun, win = windowS)
      n_bg <- bg_pp$n
    }
    events$x <- bg_pp$x
    events$y <- bg_pp$y
    events$t <- sort(runif(n_bg, min = windowT[1], max = windowT[2]))
    events$background <- rep(TRUE, n_bg)
  } else {
    events <- background_realization
    events$n <- length(events$t)
    events$t <- sort(events$t)
    events$background <- rep(TRUE, length(events$t))
    inside_space <- inside.owin(x = events$x, y = events$y, w = windowS)
    inside_time <- events$t <= windowT[2]
    valid <- inside_space & inside_time
    events <- lapply(events, function(x) { x[valid] })
  }
  if (!is.null(filtration)) {
    events <- dplyr::bind_rows(as.data.frame(events), as.data.frame(filtration))
    events <- as.list(events)
  }
  events$W <- get_covariate(events$x, events$y)

  event_queue <- data.table(time = events$t, x = events$x, y = events$y, background = events$background)
  new_events_list <- vector("list", nrow(event_queue) * (1 / (1 - K)) * 2)
  list_index <- 1
  tot <- 0
  tot_attempt <- 0
  if (!optimized) {
    while (nrow(event_queue) > 0) {
      current_event <- event_queue[1, ]
      event_queue <- event_queue[-1, ]

      num_aftershocks <- rpois(1, K)
      if (num_aftershocks > 0) {
        delay <- rexp(num_aftershocks, rate = beta)
        dist <- sqrt(rexp(num_aftershocks, rate = alpha))
        angle <- runif(num_aftershocks, min = 0, max = 2 * pi)
        new_lon <- current_event$x + dist * cos(angle)
        new_lat <- current_event$y + dist * sin(angle)
        new_time <- current_event$time + delay
        new_w <- get_covariate(new_lon, new_lat)

        inside_space <- inside.owin(x = new_lon, y = new_lat, w = windowS)
        inside_time <- (new_time >= windowT[1]) & (new_time <= windowT[2])
        valid <- inside_space & inside_time
        valid <- which(valid)
        tot_attempt <- tot_attempt + num_aftershocks
        tot <- tot + length(valid)

        if (length(valid) != 0) {
          events$x <- c(events$x, new_lon[valid])
          events$y <- c(events$y, new_lat[valid])
          events$t <- c(events$t, new_time[valid])
          events$W <- c(events$W, new_w[valid])
          events$background <- c(events$background, rep(FALSE, length(valid)))
          events$n <- events$n + length(valid)

          new_events_list[[list_index]] <- data.table(
            time = new_time[valid], x = new_lon[valid],
            y = new_lat[valid], background = rep(FALSE, length(valid)),
            W = new_w[valid]
          )
          list_index <- list_index + 1
        }
      }
      if (nrow(event_queue)[1] <= 1 & length(new_events_list) > 1) {
        event_queue <- rbindlist(list(event_queue, rbindlist(new_events_list, use.names = TRUE)))
        new_events_list <- vector("list", length(new_events_list))
        list_index <- 1
      }
    }
  } else {
    while (nrow(event_queue) > 0) {
      num_events <- nrow(event_queue)
      num_aftershocks <- rpois(num_events, K)
      idx_nonzero <- which(num_aftershocks > 0)
      if (length(idx_nonzero) == 0) break

      parent_x <- rep(event_queue$x[idx_nonzero], num_aftershocks[idx_nonzero])
      parent_y <- rep(event_queue$y[idx_nonzero], num_aftershocks[idx_nonzero])
      parent_time <- rep(event_queue$time[idx_nonzero], num_aftershocks[idx_nonzero])
      total_aftershocks <- length(parent_x)

      delay <- rexp(total_aftershocks, rate = beta)
      dist <- sqrt(rexp(total_aftershocks, rate = alpha))
      angle <- runif(total_aftershocks, min = 0, max = 2 * pi)
      new_lon <- parent_x + dist * cos(angle)
      new_lat <- parent_y + dist * sin(angle)
      new_time <- parent_time + delay
      new_w <- get_covariate(new_lon, new_lat)

      inside_space <- inside.owin(x = new_lon, y = new_lat, w = windowS)
      inside_time <- (new_time >= windowT[1]) & (new_time <= windowT[2])
      valid <- inside_space & inside_time

      if (any(valid)) {
        valid_new_lon <- new_lon[valid]
        valid_new_lat <- new_lat[valid]
        valid_new_time <- new_time[valid]
        valid_new_w <- new_w[valid]
        num_valid_new_events <- length(valid_new_lon)

        events$x <- c(events$x, valid_new_lon)
        events$y <- c(events$y, valid_new_lat)
        events$t <- c(events$t, valid_new_time)
        events$W <- c(events$W, valid_new_w)
        events$background <- c(events$background, rep(FALSE, num_valid_new_events))
        events$n <- events$n + num_valid_new_events

        event_queue <- data.table(
          time = valid_new_time, x = valid_new_lon,
          y = valid_new_lat, background = rep(FALSE, num_valid_new_events)
        )
      } else {
        event_queue <- data.table()
      }
    }
  }
  events$n <- rep(length(events$t), length(events$t))
  events
}

#' Fast Hawkes process simulation using C++ branching
#'
#' @param params List with elements mu, alpha, beta, K
#' @param windowT Numeric vector c(start, end)
#' @param windowS An owin object
#' @param background_realization Optional list to use as background events
#' @param filtration Optional data frame of history points
#' @param covariate_lookup Function(x,y) or Raster for inhomogeneous background
#' @param ... Additional arguments (ignored)
#' @return List with elements x, y, t, n, background, W
#' @export
sim_hawkes_fast <- function(params,
                            windowT,
                            windowS,
                            background_realization = NULL,
                            filtration = NULL,
                            covariate_lookup = NULL,
                            ...) {
  mu <- params$mu
  K  <- params$K

  empty_result <- list(x = numeric(0), y = numeric(0), t = numeric(0),
                       n = integer(0), background = logical(0), W = numeric(0))

  if (mu < 1e-10 && K < 1e-6 && is.null(background_realization)) {
    return(empty_result)
  }

  if (!inherits(windowS, "owin")) {
    windowS <- as.owin(windowS)
  }

  is_rect <- (windowS$type == "rectangle")
  x_min <- windowS$xrange[1]; x_max <- windowS$xrange[2]
  y_min <- windowS$yrange[1]; y_max <- windowS$yrange[2]

  alpha <- params$alpha
  beta  <- params$beta

  bg_generated <- FALSE

  if (is.null(background_realization)) {
    if (is.null(covariate_lookup)) {
      n_bg <- rpois(1, mu * (windowT[2] - windowT[1]))
      if (n_bg > 0) {
        bg_x <- runif(n_bg, x_min, x_max)
        bg_y <- runif(n_bg, y_min, y_max)
        if (!is_rect) {
          ok <- inside.owin(bg_x, bg_y, windowS)
          bg_x <- bg_x[ok]
          bg_y <- bg_y[ok]
        }
        n_ok <- length(bg_x)
        if (n_ok > 0) {
          bg_t <- sort(runif(n_ok, windowT[1], windowT[2]))
          p_x <- bg_x; p_y <- bg_y; p_t <- bg_t
          p_bg <- rep(TRUE, n_ok)
        } else {
          return(empty_result)
        }
      } else {
        if (K < 1e-6 && is.null(filtration)) return(empty_result)
        p_x <- numeric(0); p_y <- numeric(0)
        p_t <- numeric(0); p_bg <- logical(0)
      }
      bg_generated <- TRUE
    } else {
      areaS <- if (is_rect) (x_max - x_min) * (y_max - y_min) else spatstat.geom::area(windowS)
      mu_star <- mu / areaS
      pixel_fun <- function(x, y) { (windowT[2] - windowT[1]) * mu_star * covariate_lookup(x, y) }
      bg_pp <- spatstat.random::rpoispp(pixel_fun, win = windowS)
      n_bg <- bg_pp$n
      if (n_bg > 0) {
        bg_t <- sort(runif(n_bg, windowT[1], windowT[2]))
        p_x <- bg_pp$x; p_y <- bg_pp$y; p_t <- bg_t
        p_bg <- rep(TRUE, n_bg)
      } else {
        p_x <- numeric(0); p_y <- numeric(0)
        p_t <- numeric(0); p_bg <- logical(0)
      }
      bg_generated <- TRUE
    }
  } else {
    if (is.data.frame(background_realization)) {
      p_x <- background_realization$x
      p_y <- background_realization$y
      p_t <- background_realization$t
    } else {
      br <- if (is.list(background_realization)) background_realization else as.list(background_realization)
      p_x <- br$x; p_y <- br$y; p_t <- br$t
    }
    p_bg <- rep(TRUE, length(p_x))
    if (length(p_t) > 1L && is.unsorted(p_t)) {
      ord <- order(p_t)
      p_x <- p_x[ord]; p_y <- p_y[ord]; p_t <- p_t[ord]
    }
  }

  if (!is.null(filtration)) {
    if (is.data.frame(filtration)) {
      f_x <- filtration$x; f_y <- filtration$y; f_t <- filtration$t
      n_filt <- nrow(filtration)
    } else if (is.list(filtration)) {
      f_x <- filtration$x; f_y <- filtration$y; f_t <- filtration$t
      n_filt <- length(f_x)
    } else {
      f <- as.data.frame(filtration)
      f_x <- f$x; f_y <- f$y; f_t <- f$t
      n_filt <- nrow(f)
    }
    p_x <- c(p_x, f_x); p_y <- c(p_y, f_y); p_t <- c(p_t, f_t)
    p_bg <- c(p_bg, rep(TRUE, n_filt))
  }

  if (K > 1e-6 && length(p_t) > 0) {
    children <- sim_hawkes_children_cpp(
      parent_x = p_x, parent_y = p_y, parent_t = p_t,
      alpha = alpha, beta = beta, K = K,
      t_min = windowT[1], t_max = windowT[2],
      x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max
    )
    n_ch <- length(children$x)
    if (n_ch > 0) {
      combined_x <- c(p_x, children$x)
      combined_y <- c(p_y, children$y)
      combined_t <- c(p_t, children$t)
      combined_bg <- c(p_bg, rep(FALSE, n_ch))
    } else {
      combined_x <- p_x; combined_y <- p_y
      combined_t <- p_t; combined_bg <- p_bg
    }
  } else {
    combined_x <- p_x; combined_y <- p_y
    combined_t <- p_t; combined_bg <- p_bg
  }

  n_total <- length(combined_t)
  if (n_total == 0) return(empty_result)

  valid_t <- combined_t >= windowT[1] & combined_t <= windowT[2]
  if (!all(valid_t)) {
    combined_x <- combined_x[valid_t]; combined_y <- combined_y[valid_t]
    combined_t <- combined_t[valid_t]; combined_bg <- combined_bg[valid_t]
  }

  need_spatial <- !is_rect && length(combined_x) > 0
  if (need_spatial && !(bg_generated && K < 1e-6 && is.null(filtration))) {
    inside <- inside.owin(combined_x, combined_y, windowS)
    if (!all(inside)) {
      combined_x <- combined_x[inside]; combined_y <- combined_y[inside]
      combined_t <- combined_t[inside]; combined_bg <- combined_bg[inside]
    }
  }

  n_final <- length(combined_x)
  if (!is.null(covariate_lookup)) {
    if (is.function(covariate_lookup)) {
      W_vals <- covariate_lookup(combined_x, combined_y)
    } else {
      W_vals <- raster::extract(covariate_lookup, cbind(combined_x, combined_y))
    }
    W_vals[is.na(W_vals)] <- 0
  } else {
    W_vals <- numeric(n_final)
  }

  list(
    x = combined_x, y = combined_y, t = combined_t,
    n = rep(n_final, n_final),
    background = combined_bg, W = W_vals
  )
}

#' Generate inhomogeneous Hawkes process with different parameters per region
#'
#' @param Omega The full state space (owin or convertible)
#' @param partition A spatstat tess object
#' @param time_window Numeric vector c(start, end)
#' @param partition_processes Character vector of process names per tile
#' @param hawkes_params Named list of parameter lists (e.g. list(control=..., treated=...))
#' @param state_spaces Optional list of precomputed owin state spaces per process
#' @param filtration Optional data frame of history points
#' @param space_triggering Logical; if TRUE triggering depends on spatial location
#' @param ... Additional arguments passed to sim_hawkes_fast
#' @return A data.table with columns x, y, t, n, background, process, location_process
#' @export
generate_inhomogeneous_hawkes <- function(Omega,
                                          partition,
                                          time_window,
                                          partition_processes,
                                          hawkes_params,
                                          state_spaces = NULL,
                                          filtration = NULL,
                                          space_triggering = FALSE,
                                          ...) {
  processes <- unique(partition_processes)

  if (!inherits(Omega, "owin")) {
    Omega <- as.owin(Omega)
  }

  if (is.null(state_spaces)) {
    state_spaces <- lapply(processes, function(p) {
      idx <- which(partition_processes == p)
      as.owin(partition[idx])
    })
  }

  dots <- list(...)
  has_covariate <- !is.null(dots$covariate_lookup) && is.function(dots$covariate_lookup)
  filt_by_proc_precomputed <- dots$filt_by_proc

  all_bg_x <- list(); all_bg_y <- list(); all_bg_t <- list()
  all_bg_w <- list(); all_bg_proc <- list()
  for (k in seq_along(processes)) {
    p <- processes[k]
    hp <- hawkes_params[[p]]
    if (hp$mu < 1e-10) next
    tmp_hp <- hp; tmp_hp$K <- 0
    ev <- sim_hawkes_fast(
      params = tmp_hp, windowT = time_window, windowS = state_spaces[[k]],
      background_realization = NULL, filtration = NULL, ...
    )
    n_ev <- length(ev$x)
    if (n_ev == 0) next
    all_bg_x[[p]] <- ev$x
    all_bg_y[[p]] <- ev$y
    all_bg_t[[p]] <- ev$t
    all_bg_w[[p]] <- if (length(ev$W) == n_ev) ev$W else numeric(n_ev)
    all_bg_proc[[p]] <- p
  }

  filt_by_proc <- filt_by_proc_precomputed
  if (is.null(filt_by_proc) && !is.null(filtration)) {
    if (is.null(filtration$location_process)) {
      filtration$location_process <- partition_processes[tileindex(filtration$x, filtration$y, partition)]
    }
    if (is.data.frame(filtration)) {
      filt_by_proc <- split(filtration, filtration$location_process)
    } else {
      filt_df <- as.data.frame(filtration)
      filt_by_proc <- split(filt_df, filt_df$location_process)
    }
  }

  out_x <- list(); out_y <- list(); out_t <- list()
  out_bg <- list(); out_w <- list(); out_proc <- list()
  out_loc <- list(); out_tile <- list()
  idx <- 0L

  for (k in seq_along(processes)) {
    p <- processes[k]
    f <- if (!is.null(filt_by_proc) && p %in% names(filt_by_proc)) filt_by_proc[[p]] else NULL

    if (is.null(all_bg_x[[p]])) {
      hp <- hawkes_params[[p]]
      if (hp$K < 1e-6 && is.null(f)) next
      bg_realization <- list(x = numeric(0), y = numeric(0), t = numeric(0))
    } else {
      bg_realization <- list(x = all_bg_x[[p]], y = all_bg_y[[p]], t = all_bg_t[[p]])
    }

    new_events <- sim_hawkes_fast(
      params = hawkes_params[[p]], windowT = time_window, windowS = Omega,
      background_realization = bg_realization,
      filtration = f, ...
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
    out_x[[idx]] <- new_events$x
    out_y[[idx]] <- new_events$y
    out_t[[idx]] <- new_events$t
    out_bg[[idx]] <- new_events$background
    out_w[[idx]] <- w_vals
    out_proc[[idx]] <- rep(p, n_new)
    out_loc[[idx]] <- loc_proc
    out_tile[[idx]] <- tile_idx
  }

  if (idx == 0L) {
    events <- list(x = numeric(0), y = numeric(0), t = numeric(0),
                   n = integer(0), background = logical(0), W = numeric(0),
                   process = character(0), location_process = character(0),
                   tile_index = integer(0))
    setDT(events)
    return(events)
  }

  events <- list(
    x = unlist(out_x, use.names = FALSE),
    y = unlist(out_y, use.names = FALSE),
    t = unlist(out_t, use.names = FALSE),
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

#' Generate superposition of independent Hawkes processes
#'
#' @param Omega The full state space
#' @param partition A spatstat tess object
#' @param time_window Numeric vector c(start, end)
#' @param partition_processes Character vector of process names per tile
#' @param hawkes_params Named list of parameter lists
#' @param ... Additional arguments passed to sim_hawkes
#' @return List with superposed data and complete data
#' @export
generate_hawkes_superposition <- function(Omega,
                                          partition,
                                          time_window,
                                          partition_processes,
                                          hawkes_params,
                                          ...) {
  sims <- lapply(unique(partition_processes), function(x) {
    events <- sim_hawkes(hawkes_params[[x]], time_window, Omega, ...)
    events$tile <- tileindex(events$x, events$y, partition)
    events <- as.data.frame(events)
    return(events)
  })
  names(sims) <- unique(partition_processes)
  superposed <- lapply(1:length(partition_processes), function(i) {
    process <- partition_processes[i]
    sim <- sims[[process]]
    keep <- as.numeric(sim$tile) == i
    return(sim[keep, ])
  })
  superposed <- rbindlist(superposed)
  return(list(superposed = superposed, complete_data = sims))
}
