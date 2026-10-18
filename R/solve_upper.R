#' Efficacy bound at one look under the null hypothesis
#'
#' Solves for c so that the probability, under the null hypothesis, of
#' continuing to look k and then crossing c equals \code{spend}.
#'
#' @param k Look index.
#' @param lower_prev,upper_prev Bounds at looks 1 to k - 1 (lower bounds are
#'   -Inf for non-binding futility).
#' @param spend Alpha to spend at look k.
#' @param corr Correlation matrix of the statistics under the null.
#' @return The efficacy bound (Inf when nothing is spent).
#' @noRd
.solve_upper <- function(k, lower_prev, upper_prev, spend, corr) {
  if (!(spend > 0)) {
    return(Inf)
  }
  if (k == 1L) {
    return(stats::qnorm(spend, lower.tail = FALSE))
  }
  sig <- corr[seq_len(k), seq_len(k), drop = FALSE]
  f <- function(x) {
    .pmvn(c(lower_prev, x), c(upper_prev, Inf), rep(0, k), sig) - spend
  }
  hi <- stats::qnorm(spend, lower.tail = FALSE)
  f_hi <- f(hi)
  it <- 0L
  while (f_hi > 0 && it < 50L) {
    hi <- hi + 0.05
    f_hi <- f(hi)
    it <- it + 1L
  }
  lo <- hi - 1
  f_lo <- f(lo)
  it <- 0L
  while (f_lo < 0) {
    lo <- lo - 1
    f_lo <- f(lo)
    it <- it + 1L
    if (it > 30L) {
      stop("Cannot find an efficacy bound at look ", k,
           ": the alpha to spend exceeds the probability of continuing.",
           call. = FALSE)
    }
  }
  if (f_hi == 0) {
    return(hi)
  }
  stats::uniroot(f, c(lo, hi), f.lower = f_lo, f.upper = f_hi, tol = 1e-10)$root
}
