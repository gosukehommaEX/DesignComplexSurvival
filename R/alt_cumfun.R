#' Cumulative functions of the per-patient influence function
#'
#' For a patient in group z, the influence function of the log-rank score
#' at follow-up X is delta g_z(X) - L_z(X), where g_1 = 1 - q1, g_0 = -q1
#' and L_z(x) is the integral from 0 to x of g_z(s) lbar(s). This function
#' returns, at the points \code{x}, the values of L_z, of H_z (integral of
#' (g_z - L_z) e_z - L_z c_z) and of J_z (integral of (g_z - L_z)^2 e_z +
#' L_z^2 c_z), and the probabilities R_z of being at risk. L_z is integrated
#' inside each piece with the Gauss-Legendre cumulative integration matrix.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param x Sorted positive evaluation points.
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @param breaks Break points of the intensities.
#' @param hmax Maximum piece length.
#' @param n Number of nodes per piece.
#' @return A matrix with columns \code{L0}, \code{L1}, \code{H0}, \code{H1},
#'   \code{J0}, \code{J1}, \code{R0} and \code{R1}.
#' @noRd
.alt_cumfun <- function(scenario, x, endpoint, breaks, hmax, n = 8L) {
  g <- .gl_pieces(0, max(x), c(x, breaks), hmax = hmax, n = n)
  np <- length(g$lo)
  f <- .alt_funs(scenario, g$x, endpoint)
  qm <- .gl_cumint(n)
  wref <- .gl_rule(n)$w
  half <- (g$hi - g$lo) / 2
  idx <- match(x, g$hi)
  if (anyNA(idx)) {
    stop("Internal error: integration grid mismatch.", call. = FALSE)
  }
  out <- matrix(0, length(x), 8L,
                dimnames = list(NULL, c("L0", "L1", "H0", "H1", "J0", "J1", "R0", "R1")))
  for (z in 0:1) {
    gz <- if (z == 1L) 1 - f[, "q1"] else -f[, "q1"]
    fm <- matrix(gz * f[, "lbar"], nrow = n)
    tot <- half * colSums(wref * fm)
    l_lo <- c(0, cumsum(tot))[seq_len(np)]
    l_node <- rep(l_lo, each = n) + as.vector(qm %*% fm) * rep(half, each = n)
    ez <- f[, if (z == 1L) "e1" else "e0"]
    cz <- f[, if (z == 1L) "c1" else "c0"]
    hh <- (gz - l_node) * ez - l_node * cz
    hj <- (gz - l_node)^2 * ez + l_node^2 * cz
    h_cum <- cumsum(colSums(matrix(g$w * hh, nrow = n)))
    j_cum <- cumsum(colSums(matrix(g$w * hj, nrow = n)))
    out[, paste0("L", z)] <- cumsum(tot)[idx]
    out[, paste0("H", z)] <- h_cum[idx]
    out[, paste0("J", z)] <- j_cum[idx]
  }
  fx <- .alt_funs(scenario, x, endpoint)
  out[, "R0"] <- fx[, "R0"]
  out[, "R1"] <- fx[, "R1"]
  out
}
