#' Value of a piecewise constant function
#'
#' @param x Evaluation points (nonnegative).
#' @param time Interval start times (first value 0, increasing).
#' @param value Values in each interval (same length as \code{time}).
#' @return Values at \code{x}; intervals are closed on the left.
#' @noRd
.pw_value <- function(x, time, value) {
  value[pmax(findInterval(x, time), 1L)]
}
