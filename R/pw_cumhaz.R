#' Integral of a piecewise constant function from 0
#'
#' @param x Evaluation points; negative values are treated as 0.
#' @param time Interval start times (first value 0, increasing).
#' @param rate Values in each interval (same length as \code{time}).
#' @return The integral from 0 to \code{x}.
#' @noRd
.pw_cumhaz <- function(x, time, rate) {
  x <- pmax(x, 0)
  idx <- pmax(findInterval(x, time), 1L)
  k <- length(time)
  cum <- c(0, cumsum(rate[-k] * diff(time)))
  cum[idx] + rate[idx] * (x - time[idx])
}
