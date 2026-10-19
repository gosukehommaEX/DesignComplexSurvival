#' Enrolled fraction at calendar times (vectorized)
#'
#' Vectorized version of \code{.acc_frac()}: the fraction of the sample size
#' enrolled by each calendar time in \code{x} (0 for negative times).
#'
#' @param acc Output of \code{.acc_setup()}.
#' @param x Calendar times.
#' @return A numeric vector.
#' @noRd
.acc_cum <- function(acc, x) {
  nx <- length(x)
  if (nx == 0L) {
    return(numeric(0))
  }
  a <- acc$alpha
  lo <- a[-length(a)]
  len <- diff(a)
  d <- outer(x, lo, "-")
  d <- pmin(pmax(d, 0), matrix(len, nrow = nx, ncol = length(len), byrow = TRUE))
  as.vector(d %*% acc$rho)
}
