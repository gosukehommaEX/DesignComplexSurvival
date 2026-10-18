#' Fraction of the sample size enrolled by calendar time
#'
#' @param acc Output of \code{.acc_setup()}.
#' @param x Calendar times.
#' @return Numeric vector of enrolled fractions.
#' @noRd
.acc_frac <- function(acc, x) {
  a <- acc$alpha
  lo <- a[-length(a)]
  hi <- a[-1L]
  vapply(x, function(z) sum(acc$rho * pmax(pmin(z, hi) - lo, 0)), numeric(1))
}
