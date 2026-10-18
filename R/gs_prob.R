#' Stopping probabilities of a group sequential design
#'
#' For statistics with a multivariate normal distribution, computes at each
#' look the probability of continuing to the look and then crossing the
#' efficacy bound (\code{prob_upper}) or falling below the futility bound
#' (\code{prob_lower}).
#'
#' @param lower Futility bounds (\code{-Inf} for none). At the last look the
#'   bound is usually equal to the efficacy bound.
#' @param upper Efficacy bounds (\code{Inf} for none).
#' @param mean Mean vector (recycled if of length 1).
#' @param sigma Covariance matrix.
#'
#' @return A data frame with \code{look}, \code{prob_upper},
#'   \code{prob_lower} and their cumulative sums \code{cum_upper} and
#'   \code{cum_lower}.
#'
#' @seealso \code{\link{gs_bounds}}
#'
#' @examples
#' b <- gs_bounds(c(0.5, 1))
#' corr <- matrix(c(1, sqrt(0.5), sqrt(0.5), 1), 2)
#' gs_prob(b$lower, b$upper, mean = 0, sigma = corr)
#' gs_prob(b$lower, b$upper, mean = c(2, 2.9), sigma = corr)
#'
#' @export
gs_prob <- function(lower, upper, mean = 0, sigma) {
  k_max <- length(upper)
  if (!is.numeric(upper) || !is.numeric(lower) || length(lower) != k_max ||
      anyNA(lower) || anyNA(upper)) {
    stop("'lower' and 'upper' must be numeric vectors of the same length.",
         call. = FALSE)
  }
  if (length(mean) == 1L) {
    mean <- rep(mean, k_max)
  }
  if (!is.numeric(mean) || length(mean) != k_max || anyNA(mean)) {
    stop("'mean' must have length 1 or one value per look.", call. = FALSE)
  }
  .check_sigma(sigma, k_max, "sigma")
  pu <- numeric(k_max)
  pl <- numeric(k_max)
  for (k in seq_len(k_max)) {
    prev <- seq_len(k - 1L)
    sig <- sigma[seq_len(k), seq_len(k), drop = FALSE]
    mu <- mean[seq_len(k)]
    if (is.finite(upper[k])) {
      pu[k] <- .pmvn(c(lower[prev], upper[k]), c(upper[prev], Inf), mu, sig)
    }
    if (is.finite(lower[k])) {
      pl[k] <- .pmvn(c(lower[prev], -Inf), c(upper[prev], lower[k]), mu, sig)
    }
  }
  data.frame(look = seq_len(k_max), prob_upper = pu, prob_lower = pl,
             cum_upper = cumsum(pu), cum_lower = cumsum(pl))
}
