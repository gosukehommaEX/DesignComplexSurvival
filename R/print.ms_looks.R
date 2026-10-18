#' Print an ms_looks object
#'
#' @param x An \code{ms_looks} object.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' print(ms_looks("os", fraction = c(0.5, 1)))
#'
#' @export
print.ms_looks <- function(x, ...) {
  if (x$trigger == "calendar") {
    cat(sprintf("%d look(s) at calendar times: %s\n", x$k,
                paste(format(x$time), collapse = ", ")))
  } else {
    cat(sprintf("%d look(s) driven by %s events at fractions: %s\n", x$k,
                toupper(x$trigger), paste(format(x$fraction), collapse = ", ")))
  }
  invisible(x)
}
