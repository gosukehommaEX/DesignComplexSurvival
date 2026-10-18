#' Futility bound at one look by beta spending
#'
#' Solves for a so that the probability, under the alternative, of
#' continuing to look k and then falling below a equals \code{spend}. When
#' even a = c_k spends less than \code{spend}, returns c_k.
#'
#' @param k Look index.
#' @param lower_prev Futility bounds at looks 1 to k - 1.
#' @param upper Efficacy bounds at looks 1 to k.
#' @param spend Beta to spend at look k.
#' @param mean Mean vector under the alternative.
#' @param sigma Covariance matrix under the alternative.
#' @return The futility bound (-Inf when nothing is spent).
#' @noRd
.solve_lower <- function(k, lower_prev, upper, spend, mean, sigma) {
  if (!(spend > 0)) {
    return(-Inf)
  }
  sd_k <- sqrt(sigma[k, k])
  if (k == 1L) {
    return(min(mean[1L] + sd_k * stats::qnorm(spend), upper[1L]))
  }
  sig <- sigma[seq_len(k), seq_len(k), drop = FALSE]
  mu <- mean[seq_len(k)]
  up_prev <- upper[seq_len(k - 1L)]
  f <- function(x) {
    .pmvn(c(lower_prev, -Inf), c(up_prev, x), mu, sig) - spend
  }
  hi <- if (is.finite(upper[k])) upper[k] else mean[k] + 12 * sd_k
  f_hi <- f(hi)
  if (f_hi <= 0) {
    return(upper[k])
  }
  lo <- min(mean[k] + sd_k * stats::qnorm(spend), hi - 1)
  f_lo <- f(lo)
  it <- 0L
  while (f_lo > 0) {
    lo <- lo - 1
    f_lo <- f(lo)
    it <- it + 1L
    if (it > 50L) {
      stop("Cannot find a futility bound at look ", k, ".", call. = FALSE)
    }
  }
  stats::uniroot(f, c(lo, hi), f.lower = f_lo, f.upper = f_hi, tol = 1e-10)$root
}
