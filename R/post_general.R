#' Transitions after progression with a switch clock or an AFT effect
#'
#' Quantities of the post-progression process (states progressed without
#' switch 2, switched 3 and dead 4) between the clock values \code{from}
#' (progression) and \code{to}, for \code{switch_clock = "switch"} and for an
#' accelerated failure time effect of the switch (\code{switch_af} not 1).
#' After a switch at the clock value sigma, the death hazard at the clock
#' value tau is d(tau - sigma) with \code{"switch"} (d = \code{death_switch})
#' and phi c(sigma + phi (tau - sigma)) with the acceleration factor phi
#' (c = \code{death_prog}); S3(sigma, tau) is the survival after the switch.
#'
#' The probability of having switched after progression and being alive at
#' tau is the integral over sigma from \code{from} to tau of
#' r(sigma) exp(-(Y(sigma) - Y(from))) S3(sigma, tau), where r is the switch
#' hazard and Y the cumulative hazard of leaving state 2; the death density
#' of these patients has the death hazard after the switch as an extra
#' factor. The range of sigma is split into the pieces on which sigma stays
#' in one interval j of \code{time_post} and the argument of the death hazard
#' after the switch (tau - sigma, or sigma + phi (tau - sigma)) stays in one
#' interval k. On a piece the integrand is r_j times the exponential of a
#' linear function with slope d_k - y_j (\code{"switch"}) or
#' c_j - (1 - phi) c_k - y_j (AFT effect), so its integral has a closed form.
#' The exponential is evaluated at the end of the piece where it is larger,
#' so that its exponent is nonpositive and the remaining factor lies between
#' 0 and 1.
#'
#' @param model An \code{ms_model} object with \code{switch_clock = "switch"}
#'   or with \code{switch_af} different from 1.
#' @param from,to Numeric vectors of equal length with \code{from <= to}.
#' @return A list as \code{.post_trans()}: \code{r22}, \code{r23},
#'   \code{r33}, the death hazards \code{c_to} (state 2) and \code{d_to}
#'   (after a switch at progression) at \code{to}, and the death density
#'   \code{f23} of the patients who switched after progression.
#' @noRd
.post_general <- function(model, from, to) {
  tp <- model$time_post
  r <- model$switch_rate
  cc <- model$death_prog
  y <- r + cc
  aft <- !is.null(model$switch_af)
  phi <- if (aft) model$switch_af else 1
  d <- model$death_switch
  np <- length(tp)
  nf <- length(from)
  # log survival and death hazard at tau after a switch at sig
  log_s3 <- function(sig, tau) {
    if (aft) {
      -(.pw_cumhaz(sig + phi * (tau - sig), tp, cc) - .pw_cumhaz(sig, tp, cc))
    } else {
      -.pw_cumhaz(tau - sig, tp, d)
    }
  }
  h3 <- function(sig, tau) {
    if (aft) {
      phi * .pw_value(sig + phi * (tau - sig), tp, cc)
    } else {
      .pw_value(tau - sig, tp, d)
    }
  }
  cy_from <- .pw_cumhaz(from, tp, y)
  r22 <- exp(-(.pw_cumhaz(to, tp, y) - cy_from))
  r33 <- exp(log_s3(from, to))
  out <- list(r22 = r22, r23 = numeric(nf), r33 = r33,
              c_to = .pw_value(to, tp, cc), d_to = h3(from, to), f23 = numeric(nf))
  if (all(r == 0) || nf == 0L) {
    return(out)
  }
  # pieces (j, k) for all pairs; j varies fastest within k
  up <- c(tp[-1L], Inf)
  ii <- rep.int(seq_len(nf), np * np)
  jj <- rep.int(rep(seq_len(np), each = nf), np)
  kk <- rep(seq_len(np), each = nf * np)
  fr <- from[ii]
  ta <- to[ii]
  if (aft) {
    e1 <- (tp[kk] - phi * ta) / (1 - phi)
    e2 <- (up[kk] - phi * ta) / (1 - phi)
    lk <- pmin(e1, e2)
    uk <- pmax(e1, e2)
  } else {
    lk <- ta - up[kk]
    uk <- ta - tp[kk]
  }
  lo <- pmax(fr, tp[jj], lk)
  hi <- pmin(ta, up[jj], uk)
  ok <- hi > lo & r[jj] > 0
  if (!any(ok)) {
    return(out)
  }
  ii <- ii[ok]
  jj <- jj[ok]
  kk <- kk[ok]
  lo <- lo[ok]
  ta <- ta[ok]
  ln <- hi[ok] - lo
  slope <- if (aft) cc[jj] - (1 - phi) * cc[kk] - y[jj] else d[kk] - y[jj]
  a_lo <- -(.pw_cumhaz(lo, tp, y) - cy_from[ii]) + log_s3(lo, ta)
  piece <- r[jj] * exp(a_lo + pmax(slope, 0) * ln) * ln * .expm1_ratio(-abs(slope) * ln)
  hk <- if (aft) phi * cc[kk] else d[kk]
  m <- .sum_by(cbind(piece, piece * hk), ii, nf)
  out$r23 <- m[, 1L]
  out$f23 <- m[, 2L]
  out
}
