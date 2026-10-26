#' Values or integrals of piecewise polynomials
#'
#' Evaluates the polynomials of \code{.pw_poly()} at the points \code{x}, or
#' their integrals from the start of the first piece to \code{x}. Points
#' outside the pieces are moved to the nearest end.
#'
#' @param pp Output of \code{.pw_poly()}.
#' @param x Points.
#' @param integral Whether to return the integrals.
#' @return Matrix with one row per point and one column per function.
#' @noRd
.pw_poly_eval <- function(pp, x, integral = FALSE) {
  np <- length(pp$lo)
  n <- pp$n
  nf <- dim(pp$coef)[3L]
  x <- pmin(pmax(x, pp$lo[1L]), pp$hi[np])
  k <- pmin(pmax(findInterval(x, pp$lo), 1L), np)
  half <- (pp$hi[k] - pp$lo[k]) / 2
  xi <- pmin(pmax((x - (pp$lo[k] + pp$hi[k]) / 2) / half, -1), 1)
  if (integral) {
    # integral of P_j from -1 to xi: xi + 1 for j = 0, (P_{j+1} - P_{j-1}) / (2j + 1)
    pm <- .leg_basis(xi, n)
    basis <- matrix(0, length(x), n)
    basis[, 1L] <- xi + 1
    for (j in seq_len(n - 1L)) {
      basis[, j + 1L] <- (pm[, j + 2L] - pm[, j]) / (2 * j + 1)
    }
  } else {
    basis <- .leg_basis(xi, n - 1L)
  }
  out <- matrix(0, length(x), nf)
  for (f in seq_len(nf)) {
    cf <- base::t(matrix(pp$coef[, k, f], nrow = n))
    out[, f] <- rowSums(basis * cf)
    if (integral) {
      out[, f] <- pp$cum_lo[k, f] + half * out[, f]
    }
  }
  out
}
