#' Multivariate normal probability of a rectangle
#'
#' Uses \code{mvtnorm::pmvnorm()} with the Miwa algorithm, which is
#' deterministic, so that root finding on the bounds is stable. The number
#' of grid points can be changed with
#' \code{options(DesignComplexSurvival.miwa_steps = ...)} (default 128).
#'
#' When some coordinates have two finite limits and others only one, the
#' Miwa algorithm of mvtnorm replaces the infinite limits by plus or minus
#' 1000 on the standardized scale and issues a warning (\code{probval.Miwa()}
#' in mvtnorm). The same replacement is made here beforehand, which gives the
#' same probability without the warning.
#'
#' @param lower,upper Limits (may be infinite).
#' @param mean Mean vector.
#' @param sigma Covariance matrix.
#' @return The probability.
#' @noRd
.pmvn <- function(lower, upper, mean, sigma) {
  k <- length(lower)
  if (k == 1L) {
    sd <- sqrt(sigma[1L, 1L])
    if (is.infinite(upper)) {
      return(stats::pnorm(lower, mean, sd, lower.tail = FALSE))
    }
    if (is.infinite(lower)) {
      return(stats::pnorm(upper, mean, sd))
    }
    return(max(stats::pnorm(upper, mean, sd) - stats::pnorm(lower, mean, sd), 0))
  }
  lo_inf <- is.infinite(lower)
  up_inf <- is.infinite(upper)
  if (any(!lo_inf & !up_inf) && any(xor(lo_inf, up_inf))) {
    sd <- sqrt(diag(sigma))
    i <- lo_inf & !up_inf
    lower[i] <- mean[i] - 1000 * sd[i]
    i <- up_inf & !lo_inf
    upper[i] <- mean[i] + 1000 * sd[i]
  }
  steps <- getOption("DesignComplexSurvival.miwa_steps", 128L)
  p <- mvtnorm::pmvnorm(lower = lower, upper = upper, mean = mean,
                        sigma = sigma, algorithm = mvtnorm::Miwa(steps = steps))
  max(as.numeric(p), 0)
}
