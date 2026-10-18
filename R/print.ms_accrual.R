#' Print an ms_accrual object
#'
#' @param x An \code{ms_accrual} object.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' print(ms_accrual(time = c(0, 3), rate = c(10, 30)))
#'
#' @export
print.ms_accrual <- function(x, ...) {
  if (x$type == "rate") {
    cat("Piecewise uniform accrual by rates (patients per unit time):\n")
    print(data.frame(time = x$time, rate = x$rate), row.names = FALSE)
  } else {
    cat(sprintf("Piecewise uniform accrual by proportions, duration %s:\n",
                format(x$duration)))
    print(data.frame(time = x$time, prop = x$prop), row.names = FALSE)
  }
  invisible(x)
}
