#' Check a covariance matrix
#'
#' @param x Value to check.
#' @param k Required dimension.
#' @param name Argument name for messages.
#' @return \code{TRUE} invisibly; an error otherwise.
#' @noRd
.check_sigma <- function(x, k, name) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != k || ncol(x) != k ||
      anyNA(x) || any(diag(x) <= 0) || max(abs(x - t(x))) > 1e-10) {
    stop(sprintf("'%s' must be a symmetric %d by %d matrix with a positive diagonal.",
                 name, k, k), call. = FALSE)
  }
  invisible(TRUE)
}
