#' Group functions at a calendar time when switching starts at a calendar time
#'
#' For the follow-up times \code{s}, the integrals over the entry time
#' \code{e <= t - s} (with respect to the accrual distribution) of the OS
#' survival function and density of a patient who can switch from the
#' follow-up time \code{o(e) = max(tc - e, 0)} on. The entry range is split at
#' \code{tc - s} and \code{tc}: before \code{tc - s} the patient cannot switch
#' by \code{s} (the model without switching), from \code{tc} on the patient
#' can switch from randomization (the model itself), and only in between is
#' the integral over the entry time computed by Gauss-Legendre rules with
#' \code{.os_avail()}. The rules end at the accrual breaks and at the entry
#' times where (s, o) lies on a line of \code{.cal_lines()}. A group without
#' switching gives \code{S(s) A(t - s)} as in the other computations.
#'
#' @param model An \code{ms_model} object.
#' @param acc Output of \code{.acc_setup()}.
#' @param t Calendar time (one value).
#' @param s Follow-up times, not larger than \code{t}.
#' @param tc Calendar time from which switching is possible (positive).
#' @param hmax Largest piece of the rules over the entry time.
#' @return A list with \code{Y} (probability of being enrolled by \code{t - s}
#'   and alive at \code{s}) and \code{D} (the same with the OS density), both
#'   without the dropout and the allocation probability.
#' @noRd
.cal_group <- function(model, acc, t, s, tc, hmax) {
  top <- pmin(t - s, acc$duration)
  ep <- .ep_eval(model, s, "os")
  a_all <- .acc_cum(acc, pmax(top, 0))
  if (!.has_switch(model)) {
    return(list(Y = ep$surv * a_all, D = ep$dens * a_all))
  }
  m0 <- model
  m0$switch_prob[] <- 0
  m0$switch_rate[] <- 0
  e0 <- .ep_eval(m0, s, "os")
  a1 <- .acc_cum(acc, pmax(pmin(top, tc - s), 0))
  a3 <- pmax(a_all - .acc_cum(acc, tc), 0)
  y <- e0$surv * a1 + ep$surv * a3
  d <- e0$dens * a1 + ep$dens * a3
  lo <- pmax(tc - s, 0)
  hi <- pmin(top, tc)
  a <- acc$alpha
  ln <- .cal_lines(model)
  nodes <- lapply(seq_along(s), function(k) {
    if (!(hi[k] > lo[k])) {
      return(list(e = numeric(0), w = numeric(0)))
    }
    g <- .gl_pieces(lo[k], hi[k],
                    .cal_breaks(c(a, tc - (ln[, 1L] + ln[, 2L] * s[k])), lo[k], hi[k]),
                    hmax = hmax, n = 10L)
    list(e = g$x, w = g$w)
  })
  len <- vapply(nodes, function(z) length(z$e), integer(1))
  if (sum(len) > 0L) {
    idx <- rep.int(seq_along(s), len)
    e <- unlist(lapply(nodes, `[[`, "e"), use.names = FALSE)
    w <- unlist(lapply(nodes, `[[`, "w"), use.names = FALSE) *
      .pw_value(e, a[-length(a)], acc$rho)
    av <- .os_avail(model, s[idx], tc - e)
    mid <- .sum_by(cbind(w * av$surv, w * av$dens), idx, length(s))
    y <- y + mid[, 1L]
    d <- d + mid[, 2L]
  }
  list(Y = y, D = d)
}
