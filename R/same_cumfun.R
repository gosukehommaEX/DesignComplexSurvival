#' Cumulative functions for the moments of the per-patient functionals
#'
#' For each group and each functional r in (\code{psi}, \code{phi},
#' \code{dlt}) with r1 = rho + B (event observed at x) and r0 = B (censored at
#' x), returns at the points \code{x}: B(x); A_ab(x), the integral from 0 to
#' x of r1_a r1_b e + r0_a r0_b c; M_a(x), the integral of r1_a e + r0_a c;
#' and R(x), the probability of being at risk, where e and c are the
#' densities of an observed event and of an observed dropout. This extends
#' \code{.alt_cumfun()} (the score only: J = A_psi.psi, H = M_psi, L = -B_psi)
#' to the observed variance and the event indicator. Integrals inside each
#' piece use the Gauss-Legendre cumulative integration matrix.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param x Sorted positive evaluation points.
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @param breaks Break points of the intensities.
#' @param hmax Maximum piece length.
#' @param n Number of nodes per piece.
#' @return A list of two lists (control and treatment group) with \code{B}
#'   (matrix, columns \code{psi}, \code{phi}, \code{dlt}), \code{A} (matrix,
#'   columns such as \code{"psi.phi"}), \code{M} (matrix) and \code{R}
#'   (vector).
#' @noRd
.same_cumfun <- function(scenario, x, endpoint, breaks, hmax, n = 8L) {
  keys <- c("psi", "phi", "dlt")
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
  # integral from 0 of the node values v: at the nodes and at the piece ends
  cumint <- function(v) {
    fm <- matrix(v, nrow = n)
    tot <- half * colSums(wref * fm)
    start <- c(0, cumsum(tot))[seq_len(np)]
    list(node = rep(start, each = n) + as.vector(qm %*% fm) * rep(half, each = n),
         end = cumsum(tot))
  }
  # integral from 0 to each piece end of the node values v
  cumend <- function(v) {
    cumsum(colSums(matrix(g$w * v, nrow = n)))
  }
  q1 <- f[, "q1"]
  lbar <- f[, "lbar"]
  fx <- .alt_funs(scenario, x, endpoint)
  out <- vector("list", 2L)
  for (z in 0:1) {
    gz <- if (z == 1L) 1 - q1 else -q1
    lc <- cumint(gz * lbar)
    kc <- cumint((1 - 2 * q1) * gz * lbar)
    b_node <- cbind(psi = -lc$node, phi = kc$node, dlt = 0)
    r1 <- cbind(psi = gz, phi = q1 * (1 - q1), dlt = 1) + b_node
    r0 <- b_node
    ez <- f[, if (z == 1L) "e1" else "e0"]
    cz <- f[, if (z == 1L) "c1" else "c0"]
    a_mat <- matrix(0, length(x), 9L,
                    dimnames = list(NULL, as.vector(outer(keys, keys, paste, sep = "."))))
    for (k1 in keys) {
      for (k2 in keys) {
        a_mat[, paste(k1, k2, sep = ".")] <-
          cumend(r1[, k1] * r1[, k2] * ez + r0[, k1] * r0[, k2] * cz)[idx]
      }
    }
    m_mat <- matrix(0, length(x), 3L, dimnames = list(NULL, keys))
    for (k in keys) {
      m_mat[, k] <- cumend(r1[, k] * ez + r0[, k] * cz)[idx]
    }
    b_end <- cbind(psi = -lc$end[idx], phi = kc$end[idx], dlt = 0)
    out[[z + 1L]] <- list(B = b_end, A = a_mat, M = m_mat,
                          R = fx[, if (z == 1L) "R1" else "R0"])
  }
  out
}
