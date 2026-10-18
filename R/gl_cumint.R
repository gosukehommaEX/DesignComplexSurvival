#' Cumulative integration matrix on Gauss-Legendre nodes
#'
#' Returns the matrix Q such that the integral of f from -1 to the i-th
#' Gauss-Legendre node equals the sum over j of Q[i, j] f(x_j), exactly for
#' polynomials of degree below n. The matrix is built from Legendre
#' polynomials and cached.
#'
#' @param n Number of nodes (at least 2).
#' @return An n by n matrix.
#' @noRd
.gl_cumint <- function(n) {
  key <- paste0("cq", n)
  if (!is.null(.dcs_env[[key]])) {
    return(.dcs_env[[key]])
  }
  x <- .gl_rule(n)$x
  # Legendre polynomials P_0, ..., P_n evaluated at the nodes
  pm <- matrix(0, n, n + 1L)
  pm[, 1L] <- 1
  pm[, 2L] <- x
  for (k in seq_len(n - 1L)) {
    pm[, k + 2L] <- ((2 * k + 1) * x * pm[, k + 1L] - k * pm[, k]) / (k + 1)
  }
  vand <- pm[, seq_len(n), drop = FALSE]
  # integral of P_k from -1 to x: x + 1 for k = 0, (P_{k+1} - P_{k-1}) / (2k + 1)
  wint <- matrix(0, n, n)
  wint[, 1L] <- x + 1
  for (k in seq_len(n - 1L)) {
    wint[, k + 1L] <- (pm[, k + 2L] - pm[, k]) / (2 * k + 1)
  }
  q <- wint %*% solve(vand)
  assign(key, q, envir = .dcs_env)
  q
}
