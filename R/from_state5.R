#' Overall survival from the switched progression-free state
#'
#' For a patient in state 5 (switched before progression, progression-free)
#' at the time \code{w} since randomization, the probability of being alive
#' at the time \code{u >= w} and the death density at \code{u}. The patient
#' dies in state 5 with hazard \code{death_switch_pre}, or progresses at v
#' with hazard \code{prog_switch} and then dies with hazard
#' \code{death_switch} on the post-progression clock (\code{.post_switched()}).
#' The integral over v is split into the pieces on which v stays in one
#' interval j of \code{time} and the post-progression clock stays in one
#' interval k of \code{time_post} (v itself with \code{clock = "forward"},
#' u - v with \code{clock = "reset"}). On a piece the integrand is the
#' exponential of a linear function with slope d_k - h5_j, where h5 is the
#' hazard of leaving state 5, so its integral has a closed form. The
#' exponential is evaluated at the end of the piece where it is larger, so
#' that its exponent is nonpositive and the remaining factor lies between 0
#' and 1.
#'
#' @param model An \code{ms_model} object.
#' @param w,u Numeric vectors of equal length with \code{w <= u}.
#' @return A list with \code{surv} and \code{dens}.
#' @noRd
.from_state5 <- function(model, w, u) {
  tm <- model$time
  tp <- model$time_post
  a5 <- model$prog_switch
  h5 <- a5 + model$death_switch_pre
  d <- model$death_switch
  forward <- model$clock == "forward"
  nw <- length(w)
  ch5_w <- .pw_cumhaz(w, tm, h5)
  stay <- exp(-(.pw_cumhaz(u, tm, h5) - ch5_w))
  surv <- stay
  dens <- stay * .pw_value(u, tm, model$death_switch_pre)
  if (all(a5 == 0) || nw == 0L) {
    return(list(surv = surv, dens = dens))
  }
  nt <- length(tm)
  np <- length(tp)
  up_t <- c(tm[-1L], Inf)
  up_p <- c(tp[-1L], Inf)
  # pieces (j, k) for all pairs; j varies fastest within k
  ii <- rep.int(seq_len(nw), nt * np)
  jj <- rep.int(rep(seq_len(nt), each = nw), np)
  kk <- rep(seq_len(np), each = nw * nt)
  uu <- u[ii]
  if (forward) {
    lk <- tp[kk]
    uk <- up_p[kk]
  } else {
    lk <- uu - up_p[kk]
    uk <- uu - tp[kk]
  }
  lo <- pmax(w[ii], tm[jj], lk)
  hi <- pmin(uu, up_t[jj], uk)
  ok <- hi > lo & a5[jj] > 0
  if (!any(ok)) {
    return(list(surv = surv, dens = dens))
  }
  ii <- ii[ok]
  jj <- jj[ok]
  kk <- kk[ok]
  lo <- lo[ok]
  uu <- uu[ok]
  ln <- hi[ok] - lo
  # log of the survival in state 5 from w to lo and in state 3 from lo to u
  a_lo <- -(.pw_cumhaz(lo, tm, h5) - ch5_w[ii])
  if (forward) {
    a_lo <- a_lo - (.pw_cumhaz(uu, tp, d) - .pw_cumhaz(lo, tp, d))
    hd <- .pw_value(uu, tp, d)
  } else {
    a_lo <- a_lo - .pw_cumhaz(uu - lo, tp, d)
    hd <- d[kk]
  }
  slope <- d[kk] - h5[jj]
  piece <- a5[jj] * exp(a_lo + pmax(slope, 0) * ln) * ln * .expm1_ratio(-abs(slope) * ln)
  m <- .sum_by(cbind(piece, piece * hd), ii, nw)
  list(surv = surv + m[, 1L], dens = dens + m[, 2L])
}
