#' Check a vector of nonnegative times
#'
#' @param x Value to check.
#' @param name Argument name for messages.
#' @param positive Whether the times must be positive.
#' @return \code{TRUE} invisibly; an error otherwise.
#' @noRd
.check_times <- function(x, name, positive = FALSE) {
  if (!is.numeric(x) || length(x) < 1L || anyNA(x) || any(!is.finite(x)) ||
      any(x < 0) || (positive && any(x <= 0))) {
    stop(sprintf("'%s' must be a vector of finite %s numbers.", name,
                 if (positive) "positive" else "nonnegative"), call. = FALSE)
  }
  invisible(TRUE)
}
