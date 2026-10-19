#' Derivative of an accrual-weighted integral with respect to calendar time
#'
#' For each calendar time t, the derivative with respect to t of the integral
#' from 0 to t of h(s) Abar(t - s) ds computed by \code{.acc_integral()},
#' where Abar is the enrolled fraction. It equals the integral of
#' h(s) alpha(t - s) ds with the accrual density alpha, that is the sum over
#' the accrual pieces of rho_j (H1(t - a_j) - H1(t - a_(j+1))), where H1(w) is
#' the integral of h from 0 to w.
#'
#' @param acc Output of \code{.acc_setup()}.
#' @param t Calendar times.
#' @param hfun Function of the follow-up time returning a matrix (one column
#'   per integrand).
#' @param breaks Break points of the integrands.
#' @param hmax Maximum piece length.
#' @param n_gl Number of Gauss-Legendre nodes per piece.
#' @return A matrix with one row per time and one column per integrand.
#' @noRd
.acc_integral_deriv <- function(acc, t, hfun, breaks, hmax, n_gl = 12L) {
  a <- acc$alpha
  nj <- length(a) - 1L
  wm <- pmax(outer(t, a, "-"), 0)
  pts <- sort(unique(wm[wm > 0]))
  if (length(pts) == 0L) {
    h0 <- as.matrix(hfun(0))
    out <- matrix(0, length(t), ncol(h0))
    colnames(out) <- colnames(h0)
    return(out)
  }
  g <- .gl_pieces(0, max(pts), c(pts, breaks), hmax = hmax, n = n_gl)
  hv <- as.matrix(hfun(g$x))
  h1 <- .col_cumsum(rowsum(g$w * hv, g$piece, reorder = TRUE))
  h1_at <- function(w) {
    out <- matrix(0, length(w), ncol(hv))
    pos <- w > 0
    if (any(pos)) {
      idx <- match(w[pos], g$hi)
      if (anyNA(idx)) {
        stop("Internal error: integration grid mismatch.", call. = FALSE)
      }
      out[pos, ] <- h1[idx, , drop = FALSE]
    }
    out
  }
  res <- matrix(0, length(t), ncol(hv))
  for (j in seq_len(nj)) {
    res <- res + acc$rho[j] * (h1_at(wm[, j]) - h1_at(wm[, j + 1L]))
  }
  colnames(res) <- colnames(hv)
  res
}
