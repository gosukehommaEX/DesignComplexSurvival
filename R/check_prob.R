#' Check a probability strictly between 0 and 1
#'
#' @param x Value to check.
#' @param name Argument name for messages.
#' @return \code{TRUE} invisibly; an error otherwise.
#' @noRd
.check_prob <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x >= 1) {
    stop(sprintf("'%s' must be a number between 0 and 1.", name), call. = FALSE)
  }
  invisible(TRUE)
}
