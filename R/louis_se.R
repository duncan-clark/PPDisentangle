#' @importFrom numDeriv grad hessian
NULL

#' Complete-data log-likelihood for a single labelling
#'
#' Evaluates the joint log-likelihood of control and treated components
#' given a labelling that splits observed data into the two processes.
#'
#' @param theta Numeric vector of length 8: c(mu_c, alpha_c, beta_c, K_c,
#'   mu_t, alpha_t, beta_t, K_t)
#' @param labelling Data frame with columns x, y, t, inferred_process and
#'   optionally W (background covariate)
#' @param treatment_time Scalar; only events with t >= treatment_time are used
#' @param statespace An owin observation window
#' @param partition A spatstat tess object
#' @param partition_processes Character vector of process names per tile
#' @param background_rate_var Column name for inhomogeneous background (or NULL)
#' @return Scalar log-likelihood value
#' @export
complete_data_loglik <- function(theta,
                                 labelling,
                                 treatment_time,
                                 statespace,
                                 partition,
                                 partition_processes,
                                 background_rate_var = NULL) {
  theta_c <- theta[1:4]
  theta_t <- theta[5:8]

  treated_idx <- (partition_processes == "treated")
  treated_state_space <- as.owin(partition[treated_idx])
  control_state_space <- as.owin(partition[!treated_idx])

  post <- labelling[labelling$t >= treatment_time, ]
  windowT <- c(treatment_time, max(post$t))

  control_data <- post[post$inferred_process == "control", ]
  treated_data <- post[post$inferred_process == "treated", ]

  ll_c <- -Inf
  ll_t <- -Inf

  if (nrow(control_data) > 1) {
    ll_c <- tryCatch(
      loglik_hawk_fast(
        params = theta_c,
        realiz = control_data,
        windowT = windowT,
        windowS = statespace,
        zero_background_region = treated_state_space,
        background_rate_var = background_rate_var
      ),
      error = function(e) -Inf
    )
  }

  if (nrow(treated_data) > 1) {
    ll_t <- tryCatch(
      loglik_hawk_fast(
        params = theta_t,
        realiz = treated_data,
        windowT = windowT,
        windowS = statespace,
        zero_background_region = control_state_space,
        background_rate_var = background_rate_var
      ),
      error = function(e) -Inf
    )
  }

  ll <- ll_c + ll_t
  if (!is.finite(ll)) ll <- -1e12
  return(ll)
}

#' Standard errors via Louis's method from SEM output
#'
#' Computes the observed information matrix using Louis's identity applied
#' to the importance-weighted labellings from the final SEM step. Returns
#' standard errors and Wald confidence intervals for all 8 Hawkes parameters
#' (4 control + 4 treated).
#'
#' Louis's identity:
#'   \code{I_obs(theta) = -E[H] - (E[g g^T] - E[g] E[g]^T)}
#' where \code{H} is the complete-data Hessian, \code{g} the complete-data
#' score, and expectations are over the conditional distribution of labels
#' given data, approximated by importance-weighted SEM samples.
#'
#' @param sem_result Output from \code{\link{adaptive_SEM}}
#' @param treatment_time Scalar treatment time
#' @param statespace An owin observation window
#' @param partition A spatstat tess object
#' @param partition_processes Character vector of process names per tile
#' @param background_rate_var Column name for inhomogeneous background (or NULL)
#' @param alpha Significance level for confidence intervals (default 0.05)
#' @param verbose Print progress
#' @return List with components:
#'   \item{se}{Named numeric vector of standard errors}
#'   \item{ci}{Data frame with Estimate, SE, Lower, Upper columns}
#'   \item{I_obs}{8x8 observed information matrix}
#'   \item{vcov}{8x8 variance-covariance matrix}
#'   \item{weights}{Importance weights used}
#' @export
louis_standard_errors <- function(sem_result,
                                  treatment_time,
                                  statespace,
                                  partition,
                                  partition_processes,
                                  background_rate_var = NULL,
                                  alpha = 0.05,
                                  verbose = TRUE) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Package 'numDeriv' is required for Louis's method. Install it with install.packages('numDeriv').")
  }

  theta_c <- unlist(sem_result$hawkes_params_control)
  theta_t <- unlist(sem_result$hawkes_params_treated)
  theta_hat <- c(theta_c, theta_t)
  names(theta_hat) <- c("mu_c", "alpha_c", "beta_c", "K_c",
                         "mu_t", "alpha_t", "beta_t", "K_t")
  d <- length(theta_hat)

  labellings <- sem_result$labellings
  M <- length(labellings)

  if (verbose) message("Computing importance weights for ", M, " labellings...")

  treated_idx <- (partition_processes == "treated")
  treated_state_space <- as.owin(partition[treated_idx])
  control_state_space <- as.owin(partition[!treated_idx])

  log_weights <- sapply(labellings, function(lab) {
    complete_data_loglik(
      theta_hat, lab, treatment_time, statespace,
      partition, partition_processes, background_rate_var
    )
  })

  log_weights <- log_weights - max(log_weights)
  weights <- exp(log_weights)
  weights <- weights / sum(weights)

  if (verbose) {
    ess <- 1 / sum(weights^2)
    message("Effective sample size: ", round(ess, 1), " / ", M)
  }

  if (verbose) message("Computing scores and Hessians for each labelling...")

  grads <- vector("list", M)
  hessians <- vector("list", M)

  for (i in seq_len(M)) {
    if (verbose && (i %% 10 == 0 || i == 1)) {
      message("  labelling ", i, "/", M)
    }

    lab_i <- labellings[[i]]

    ll_func <- function(th) {
      complete_data_loglik(
        th, lab_i, treatment_time, statespace,
        partition, partition_processes, background_rate_var
      )
    }

    grads[[i]] <- tryCatch(
      numDeriv::grad(ll_func, theta_hat),
      error = function(e) { warning("grad failed for labelling ", i); rep(0, d) }
    )

    hessians[[i]] <- tryCatch(
      numDeriv::hessian(ll_func, theta_hat),
      error = function(e) { warning("hessian failed for labelling ", i); matrix(0, d, d) }
    )
  }

  if (verbose) message("Assembling Louis's formula...")

  bar_H <- Reduce("+", mapply(function(H, w) w * H, hessians, weights, SIMPLIFY = FALSE))
  bar_G <- Reduce("+", mapply(function(g, w) w * g, grads, weights, SIMPLIFY = FALSE))

  GG_list <- mapply(function(g, w) w * (g %*% t(g)), grads, weights, SIMPLIFY = FALSE)
  bar_GG <- Reduce("+", GG_list)

  I_obs <- (-bar_H) - (bar_GG - bar_G %*% t(bar_G))

  rownames(I_obs) <- names(theta_hat)
  colnames(I_obs) <- names(theta_hat)

  vcov <- tryCatch(
    solve(I_obs),
    error = function(e) {
      warning("Information matrix is singular; using pseudoinverse.")
      eig <- eigen(I_obs, symmetric = TRUE)
      tol <- max(abs(eig$values)) * .Machine$double.eps * d
      pos <- eig$values > tol
      if (!any(pos)) return(matrix(NA, d, d))
      eig$vectors[, pos, drop = FALSE] %*%
        diag(1 / eig$values[pos], nrow = sum(pos)) %*%
        t(eig$vectors[, pos, drop = FALSE])
    }
  )
  rownames(vcov) <- names(theta_hat)
  colnames(vcov) <- names(theta_hat)

  se_raw <- diag(vcov)
  se <- ifelse(se_raw > 0, sqrt(se_raw), NA_real_)
  names(se) <- names(theta_hat)

  z <- qnorm(1 - alpha / 2)
  ci <- data.frame(
    Parameter = names(theta_hat),
    Estimate = theta_hat,
    SE = se,
    Lower = theta_hat - z * se,
    Upper = theta_hat + z * se,
    row.names = NULL
  )

  if (verbose) {
    message("\n--- Louis's Method Standard Errors ---")
    print(ci)
  }

  list(
    se = se,
    ci = ci,
    I_obs = I_obs,
    vcov = vcov,
    weights = weights
  )
}
