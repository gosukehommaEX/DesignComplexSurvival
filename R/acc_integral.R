#' Integrals weighted by the enrolled fraction
#'
#' For each calendar time t computes the integral from 0 to t of
#' Abar(t - s) h(s), where Abar is the enrolled fraction. Because Abar is
#' piecewise linear with density rho on the accrual intervals, the integral
#' equals the sum over intervals m of rho_m [H2(t - alpha_(m-1)) -
#' H2(t - alpha_m)], where H2(w) is the integral from 0 to w of (w - s) h(s)
#' (zero for w <= 0). H2 is accumulated over pieces with Gauss-Legendre
#' rules.
#'
#' @param acc Output of \code{.acc_setup()}.
#' @param t Calendar times.
#' @param hfun Function of the follow-up time returning a matrix with one
#'   column per integrand.
#' @param breaks Break points of the integrands.
#' @param hmax Maximum piece length.
#' @param n_gl Number of nodes per piece.
#' @return A matrix with one row per calendar time and one column per
#'   integrand.
#' @noRd
.acc_integral <- function(acc, t, hfun, breaks, hmax, n_gl = 12L) {
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
  np <- length(g$lo)
  i1 <- rowsum(g$w * hv, g$piece, reorder = TRUE)
  i2 <- rowsum((g$w * (g$hi[g$piece] - g$x)) * hv, g$piece, reorder = TRUE)
  h1 <- .col_cumsum(i1)
  h1_prev <- rbind(0, h1[-np, , drop = FALSE])
  h2 <- .col_cumsum((g$hi - g$lo) * h1_prev + i2)
  h2_at <- function(w) {
    out <- matrix(0, length(w), ncol(hv))
    pos <- w > 0
    if (any(pos)) {
      idx <- match(w[pos], g$hi)
      if (anyNA(idx)) {
        stop("Internal error: integration grid mismatch.", call. = FALSE)
      }
      out[pos, ] <- h2[idx, , drop = FALSE]
    }
    out
  }
  res <- matrix(0, length(t), ncol(hv))
  for (j in seq_len(nj)) {
    res <- res + acc$rho[j] * (h2_at(wm[, j]) - h2_at(wm[, j + 1L]))
  }
  colnames(res) <- colnames(hv)
  res
}
