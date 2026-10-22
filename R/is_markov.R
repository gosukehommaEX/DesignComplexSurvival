#' Whether a covariance matrix is that of a Markov sequence
#'
#' A normal vector is a Markov sequence (each coordinate depends on the
#' earlier ones only through the previous one) when the inverse of its
#' correlation matrix is tridiagonal, as for the correlation
#' \code{sqrt(t_j / t_k)} of statistics with independent increments. Used by
#' \code{.pmvn()} to choose between the Miwa and the Genz-Bretz algorithms.
#'
#' @param sigma Covariance matrix.
#' @param tol Tolerance for the elements of the inverse correlation matrix
#'   outside the tridiagonal band, relative to its largest element.
#' @return \code{TRUE} or \code{FALSE} (\code{TRUE} for at most two
#'   coordinates, \code{FALSE} when the correlation matrix cannot be
#'   inverted).
#' @noRd
.is_markov <- function(sigma, tol = 1e-8) {
  k <- nrow(sigma)
  if (k <= 2L) {
    return(TRUE)
  }
  sd <- sqrt(diag(sigma))
  prec <- tryCatch(solve(sigma / outer(sd, sd)), error = function(e) NULL)
  if (is.null(prec)) {
    return(FALSE)
  }
  off <- abs(row(prec) - col(prec)) >= 2L
  max(abs(prec[off])) <= tol * max(abs(prec))
}
