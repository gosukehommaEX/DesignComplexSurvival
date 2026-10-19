#' Moments of the per-patient functionals of one endpoint at several times
#'
#' Per-patient means E[r_a(t_i)] and products E[r_a(t_i) r_b(t_j)] of the
#' score (\code{psi}), the observed variance (\code{phi}) and the event
#' indicator (\code{dlt}) of one endpoint, summed over the groups with the
#' allocation probabilities. For a patient enrolled at e and calendar times
#' t <= t', let W = t - e and W' = t' - e. Splitting by an event or dropout
#' before W, between W and W', or after W' gives
#' E[r_a(t) r_b(t') | e] = A_ab(W) + B_a(W) (M_b(W') - M_b(W) + B_b(W') R(W'))
#' and E[r_a(t) | e] = M_a(W) + B_a(W) R(W), with the functions of
#' \code{.same_cumfun()}; these are integrated over the accrual density.
#' For the score this is the computation of \code{.lr_alt_cov()}.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param acc Output of \code{.acc_setup()}.
#' @param t Calendar times (positive).
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @return A list with \code{mean} (matrix, one row per time, columns
#'   \code{psi}, \code{phi}, \code{dlt}) and \code{prod} (array with
#'   dimensions time, time, functional, functional).
#' @noRd
.same_moments <- function(scenario, acc, t, endpoint) {
  keys <- c("psi", "phi", "dlt")
  k <- length(t)
  brk <- .scn_breaks(scenario)
  a <- acc$alpha
  hmax <- .scn_hmax(scenario, max(c(t, 1)))
  pairs <- which(upper.tri(diag(k), diag = TRUE), arr.ind = TRUE)
  enodes <- lapply(seq_len(nrow(pairs)), function(q) {
    tl <- min(t[pairs[q, ]])
    th <- max(t[pairs[q, ]])
    hi <- min(tl, acc$duration)
    g <- .gl_pieces(0, hi, c(a, th - brk, tl - brk), hmax = hmax, n = 10L)
    rho_e <- .pw_value(g$x, a[-length(a)], acc$rho)
    list(w_hi = th - g$x, w_lo = tl - g$x, w = g$w * rho_e)
  })
  pts <- sort(unique(unlist(lapply(enodes, function(z) c(z$w_hi, z$w_lo)),
                            use.names = FALSE)))
  avg <- matrix(0, k, 3L, dimnames = list(NULL, keys))
  prod <- array(0, c(k, k, 3L, 3L), dimnames = list(NULL, NULL, keys, keys))
  if (length(pts) == 0L) {
    return(list(mean = avg, prod = prod))
  }
  cf <- .same_cumfun(scenario, pts, endpoint, brk, hmax)
  p <- scenario$prob
  for (q in seq_len(nrow(pairs))) {
    i <- pairs[q, 1L]
    j <- pairs[q, 2L]
    early <- if (t[i] <= t[j]) i else j
    late <- if (t[i] <= t[j]) j else i
    z <- enodes[[q]]
    if (length(z$w) == 0L) {
      next
    }
    i_hi <- match(z$w_hi, pts)
    i_lo <- match(z$w_lo, pts)
    e_mat <- matrix(0, 3L, 3L, dimnames = list(keys, keys))
    for (grp in 0:1) {
      cg <- cf[[grp + 1L]]
      pg <- p[[grp + 1L]]
      for (k1 in keys) {
        for (k2 in keys) {
          val <- cg$A[i_lo, paste(k1, k2, sep = ".")] +
            cg$B[i_lo, k1] * (cg$M[i_hi, k2] - cg$M[i_lo, k2] +
                                cg$B[i_hi, k2] * cg$R[i_hi])
          e_mat[k1, k2] <- e_mat[k1, k2] + pg * sum(z$w * val)
        }
      }
      if (i == j) {
        for (k1 in keys) {
          val <- cg$M[i_lo, k1] + cg$B[i_lo, k1] * cg$R[i_lo]
          avg[i, k1] <- avg[i, k1] + pg * sum(z$w * val)
        }
      }
    }
    prod[early, late, , ] <- e_mat
    if (early != late) {
      prod[late, early, , ] <- base::t(e_mat)
    }
  }
  list(mean = avg, prod = prod)
}
