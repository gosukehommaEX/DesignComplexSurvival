#' Break points of a Gauss-Legendre rule, without near duplicates
#'
#' The break points of the rules of the calendar-time computations come from
#' several formulas, and the same point can appear with rounding
#' differences; a piece of length 1e-15 would get its own nodes. Returns the
#' sorted break points strictly inside (a, b) more than \code{1e-9 max(|b|, 1)}
#' away from the ends and from the previous point.
#'
#' @param x Break points.
#' @param a,b Ends of the range.
#' @return Numeric vector.
#' @noRd
.cal_breaks <- function(x, a, b) {
  tol <- 1e-9 * max(abs(b), 1)
  x <- sort(x[x > a + tol & x < b - tol])
  if (length(x) < 2L) {
    return(x)
  }
  x[c(TRUE, diff(x) > tol)]
}
