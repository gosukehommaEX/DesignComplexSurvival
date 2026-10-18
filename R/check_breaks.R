#' Check a vector of interval start times
#'
#' @param x Value to check.
#' @param name Argument name for messages.
#' @return \code{TRUE} invisibly; an error otherwise.
#' @noRd
.check_breaks <- function(x, name) {
  if (!is.numeric(x) || length(x) < 1L || anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("'%s' must be a finite numeric vector.", name), call. = FALSE)
  }
  if (x[1L] != 0) {
    stop(sprintf("'%s' must start at 0.", name), call. = FALSE)
  }
  if (length(x) > 1L && any(diff(x) <= 0)) {
    stop(sprintf("'%s' must be strictly increasing.", name), call. = FALSE)
  }
  invisible(TRUE)
}
