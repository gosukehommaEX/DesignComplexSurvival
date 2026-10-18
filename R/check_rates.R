#' Check and expand a vector of rates
#'
#' @param x Value to check (length 1 or \code{n}).
#' @param n Required length.
#' @param name Argument name for messages.
#' @param upper Upper limit of the values.
#' @return The vector expanded to length \code{n}.
#' @noRd
.check_rates <- function(x, n, name, upper = Inf) {
  if (!is.numeric(x) || length(x) < 1L || anyNA(x) || any(!is.finite(x)) ||
      any(x < 0) || any(x > upper)) {
    stop(sprintf("'%s' must contain finite values between 0 and %s.", name,
                 format(upper)), call. = FALSE)
  }
  if (length(x) == 1L) {
    x <- rep(x, n)
  }
  if (length(x) != n) {
    stop(sprintf("'%s' must have length 1 or %d.", name, n), call. = FALSE)
  }
  as.numeric(x)
}
