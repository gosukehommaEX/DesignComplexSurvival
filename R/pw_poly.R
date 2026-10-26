#' Piecewise polynomials through values at Gauss-Legendre nodes
#'
#' On each piece of a rule from \code{.gl_pieces()} with \code{n} nodes per
#' piece, the polynomial of degree \code{n - 1} through the values at the
#' nodes, in the Legendre basis of the piece. Used to evaluate functions known
#' at the nodes of one rule, and their integrals, at other points
#' (\code{.pw_poly_eval()}).
#'
#' @param g Output of \code{.gl_pieces()}.
#' @param vals Values at the nodes \code{g$x}: a vector or a matrix with one
#'   column per function.
#' @param n Number of nodes per piece of \code{g}.
#' @return A list with the piece ends \code{lo} and \code{hi}, the
#'   coefficients \code{coef} (array: degree, piece, function), the integrals
#'   \code{cum_lo} from \code{lo[1]} to the start of each piece (matrix:
#'   piece, function) and \code{n}.
#' @noRd
.pw_poly <- function(g, vals, n) {
  vals <- as.matrix(vals)
  np <- length(g$lo)
  nf <- ncol(vals)
  vinv <- solve(.leg_basis(.gl_rule(n)$x, n - 1L))
  coef <- array(0, c(n, np, nf))
  tot <- matrix(0, np, nf)
  for (f in seq_len(nf)) {
    cm <- vinv %*% matrix(vals[, f], nrow = n)
    coef[, , f] <- cm
    tot[, f] <- (g$hi - g$lo) * cm[1L, ]
  }
  cum_lo <- rbind(0, .col_cumsum(tot))[seq_len(np), , drop = FALSE]
  list(lo = g$lo, hi = g$hi, coef = coef, cum_lo = cum_lo, n = n)
}
