#' Normalize a vector of log-weights to proper weights summing to 1
#'
#' Uses the log-sum-exp trick for numerical stability.
#' Non-finite inputs (e.g. -Inf) receive zero weight.
#' Returns uniform weights if all finite values are identical.
#'
#' @param w Numeric vector of log-weights (e.g. log-likelihoods)
#' @return Numeric vector of same length as \code{w}, summing to 1
#' @export
normalize_weights <- function(w) {
  n <- length(w)
  if (n == 0) return(numeric(0))
  if (n == 1) return(1)
  ok <- is.finite(w)
  if (!any(ok)) return(rep(1 / n, n))
  w_max <- max(w[ok])
  if (w_max - min(w[ok]) < .Machine$double.eps * 100) {
    return(rep(1 / n, n))
  }
  log_w <- ifelse(ok, w - w_max, -Inf)
  weights <- exp(log_w)
  weights / sum(weights)
}

#' Generate Poisson point process data under a treatment partition
#'
#' @param Omega The state space (owin)
#' @param partition A spatstat tess object
#' @param treat_prop Proportion of tiles to treat
#' @param lambda_1 Intensity for control process
#' @param lambda_2 Intensity for treated process
#' @return List with obs_data, partition, and treated indicator
#' @export
generate_poisson_data <- function(Omega, partition, treat_prop = 0.5, lambda_1, lambda_2) {
  treated <- rbinom(partition$n, 1, treat_prop)
  control_data <- rpoispp(lambda_1, win = as.owin(partition[which(treated == 0)]))
  treated_data <- rpoispp(lambda_2, win = as.owin(partition[which(treated == 1)]))
  marks(control_data) <- 0
  marks(treated_data) <- 1
  obs_data <- superimpose(control_data, treated_data)
  obs_data <- data.frame(
    x = obs_data$x, y = obs_data$y,
    process = c("control", "treated")[marks(obs_data) + 1],
    location_process = c("control", "treated")[marks(obs_data) + 1],
    background = TRUE
  )
  return(list(obs_data = obs_data, partition = partition, treated = treated))
}
