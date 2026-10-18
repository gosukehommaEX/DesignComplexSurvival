#' Check information fractions
#'
#' @param x Value to check.
#' @return \code{TRUE} invisibly; an error otherwise.
#' @noRd
.check_info_frac <- function(x) {
  if (!is.numeric(x) || length(x) < 1L || anyNA(x) || any(x <= 0) ||
      any(x > 1 + 1e-12) || (length(x) > 1L && any(diff(x) <= 0)) ||
      abs(x[length(x)] - 1) > 1e-10) {
    stop("'info_frac' must be strictly increasing, positive and end at 1.",
         call. = FALSE)
  }
  invisible(TRUE)
}
