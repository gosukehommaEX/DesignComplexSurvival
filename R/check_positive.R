#' Check a positive finite scalar
#'
#' @param x Value to check.
#' @param name Argument name for messages.
#' @return \code{TRUE} invisibly; an error otherwise.
#' @noRd
.check_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    stop(sprintf("'%s' must be a positive number.", name), call. = FALSE)
  }
  invisible(TRUE)
}
