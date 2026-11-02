#' Conditional distribution of some statistics given others
#'
#' For a normal vector with mean \code{mean} and covariance \code{sigma},
#' the distribution of the coordinates \code{future} given that the
#' coordinates \code{given} equal \code{z}.
#'
#' @param mean Mean vector.
#' @param sigma Covariance matrix.
#' @param given,future Indices of the conditioning and of the conditioned
#'   coordinates.
#' @param z Values of the conditioning coordinates.
#' @return A list with the conditional mean vector \code{mean} and covariance
#'   matrix \code{sigma} (made exactly symmetric) of the coordinates
#'   \code{future}.
#' @noRd
.cond_moments <- function(mean, sigma, given, future, z) {
  s_fg <- sigma[future, given, drop = FALSE]
  b <- s_fg %*% solve(sigma[given, given, drop = FALSE])
  m <- mean[future] + as.vector(b %*% (z - mean[given]))
  s <- sigma[future, future, drop = FALSE] - b %*% t(s_fg)
  list(mean = m, sigma = (s + t(s)) / 2)
}
