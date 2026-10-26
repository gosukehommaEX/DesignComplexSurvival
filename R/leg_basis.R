#' Legendre polynomials
#'
#' @param x Points in [-1, 1].
#' @param deg Largest degree.
#' @return Matrix with one row per point and the columns P_0, ..., P_deg.
#' @noRd
.leg_basis <- function(x, deg) {
  out <- matrix(0, length(x), deg + 1L)
  out[, 1L] <- 1
  if (deg >= 1L) {
    out[, 2L] <- x
  }
  for (k in seq_len(max(deg - 1L, 0L))) {
    out[, k + 2L] <- ((2 * k + 1) * x * out[, k + 1L] - k * out[, k]) / (k + 1)
  }
  out
}
