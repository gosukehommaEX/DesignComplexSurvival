#' Covariance of the log-rank score under the alternative
#'
#' Per-patient covariance of the score at calendar times t (scaled so that
#' the covariance of the score is n times this value). For calendar times
#' t' <= t and a patient enrolled at e <= t', let W = t - e and W' = t' - e.
#' Conditioning on whether the event or dropout time falls before W',
#' between W' and W, or after W gives
#' E[psi(t) psi(t') | e, z] = J_z(W') - L_z(W') (H_z(W) - H_z(W')) +
#' L_z(W) L_z(W') R_z(W), with the functions of \code{.alt_cumfun()}. The
#' covariance is the integral of this over the accrual density, summed over
#' the groups with the allocation probabilities, minus U(t) U(t'). On the
#' diagonal this equals V(t) of Luo et al. (2019).
#'
#' @param scenario An \code{ms_scenario} object.
#' @param acc Output of \code{.acc_setup()}.
#' @param t Calendar times.
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @param u_pp Per-patient score means at \code{t}.
#' @return A symmetric matrix.
#' @noRd
.lr_alt_cov <- function(scenario, acc, t, endpoint, u_pp) {
  k <- length(t)
  brk <- .scn_breaks(scenario)
  a <- acc$alpha
  hmax <- .scn_hmax(scenario, max(t))
  pairs <- which(upper.tri(diag(k), diag = TRUE), arr.ind = TRUE)
  enodes <- lapply(seq_len(nrow(pairs)), function(q) {
    tl <- min(t[pairs[q, 1L]], t[pairs[q, 2L]])
    th <- max(t[pairs[q, 1L]], t[pairs[q, 2L]])
    hi <- min(tl, acc$duration)
    g <- .gl_pieces(0, hi, c(a, th - brk, tl - brk), hmax = hmax, n = 10L)
    rho_e <- .pw_value(g$x, a[-length(a)], acc$rho)
    list(w_hi = th - g$x, w_lo = tl - g$x, w = g$w * rho_e)
  })
  pts <- sort(unique(unlist(lapply(enodes, function(z) c(z$w_hi, z$w_lo)),
                            use.names = FALSE)))
  cf <- .alt_cumfun(scenario, pts, endpoint, brk, hmax)
  p <- scenario$prob
  e_mat <- matrix(0, k, k)
  for (q in seq_len(nrow(pairs))) {
    z <- enodes[[q]]
    i_hi <- match(z$w_hi, pts)
    i_lo <- match(z$w_lo, pts)
    tot <- 0
    for (grp in 0:1) {
      l_hi <- cf[i_hi, paste0("L", grp)]
      l_lo <- cf[i_lo, paste0("L", grp)]
      h_hi <- cf[i_hi, paste0("H", grp)]
      h_lo <- cf[i_lo, paste0("H", grp)]
      j_lo <- cf[i_lo, paste0("J", grp)]
      r_hi <- cf[i_hi, paste0("R", grp)]
      kz <- j_lo - l_lo * (h_hi - h_lo) + l_hi * l_lo * r_hi
      tot <- tot + p[[grp + 1L]] * sum(z$w * kz)
    }
    e_mat[pairs[q, 1L], pairs[q, 2L]] <- tot
    e_mat[pairs[q, 2L], pairs[q, 1L]] <- tot
  }
  e_mat - outer(u_pp, u_pp)
}
