#' Contributions of progressed patients
#'
#' Integrates over the time of progression u from 0 to s. With
#' \code{what = "os"} returns the parts of the overall survival function and
#' density that come from patients who progressed by s; with
#' \code{what = "states"} returns the occupation probabilities of states 2
#' (progressed without switch) and 3 (switched). Patients who progress after
#' a switch before progression (state 5) enter state 3 at progression
#' (\code{.post_switched()}). The pieces of the Gauss-Legendre rule end at
#' the break points of the integrand in u, which include, with
#' \code{switch_clock = "switch"} or an accelerated failure time effect of
#' the switch, the points where the argument of the death hazard after the
#' switch crosses a break point of \code{time_post}.
#'
#' @param model An \code{ms_model} object.
#' @param s Follow-up times since randomization.
#' @param what \code{"os"} or \code{"states"}.
#' @return A list of two numeric vectors.
#' @noRd
.os_conv <- function(model, s, what = "os") {
  ns <- length(s)
  if (ns == 0L) {
    return(if (what == "os") list(surv = numeric(0), dens = numeric(0))
           else list(p12 = numeric(0), p13 = numeric(0)))
  }
  tm <- model$time
  tp <- model$time_post
  lam <- max(c(model$prog + model$death + model$switch_pre,
               model$switch_rate + model$death_prog, model$death_switch,
               model$prog_switch + model$death_switch_pre), 1e-12)
  forward <- model$clock == "forward"
  has5 <- any(model$switch_pre > 0) && any(model$prog_switch > 0)
  sw <- identical(model$switch_clock, "switch")
  phi <- if (is.null(model$switch_af) || model$switch_af == 1) NULL else model$switch_af
  # values of the post-progression clock at which the integrand in u has breaks
  extra <- if (sw) {
    as.vector(outer(tp, tp, "+"))
  } else if (!is.null(phi)) {
    as.vector(outer(tp, tp, function(a, b) (a - (1 - phi) * b) / phi))
  } else {
    numeric(0)
  }
  nodes <- lapply(s, function(si) {
    if (!(si > 0)) {
      return(list(u = numeric(0), w = numeric(0)))
    }
    brk <- if (forward) c(tm, tp) else c(tm, si - tp)
    if (sw) {
      brk <- c(brk, si - (if (forward) tp else extra))
    } else if (!is.null(phi)) {
      brk <- c(brk, if (forward) (tp - phi * si) / (1 - phi) else si - extra)
    }
    g <- .gl_pieces(0, si, brk, hmax = max(2 / lam, si / 400), n = 10L)
    list(u = g$x, w = g$w)
  })
  len <- vapply(nodes, function(z) length(z$u), integer(1))
  idx <- rep.int(seq_len(ns), len)
  u <- as.numeric(unlist(lapply(nodes, `[[`, "u"), use.names = FALSE))
  w <- as.numeric(unlist(lapply(nodes, `[[`, "w"), use.names = FALSE))
  sv <- s[idx]
  pre <- .pre_states(model, u)
  kappa <- pre$p1 * .pw_value(u, tm, model$prog)
  pp <- .pw_value(u, tm, model$switch_prob)
  o0 <- if (forward) u else numeric(length(u))
  o1 <- if (forward) sv else sv - u
  tr <- .post_trans(model, o0, o1)
  if (has5) {
    k5 <- pre$p5 * .pw_value(u, tm, model$prog_switch)
    t5 <- .post_switched(model, o0, o1)
  }
  if (what == "os") {
    gs <- (1 - pp) * (tr$r22 + tr$r23) + pp * tr$r33
    gd <- (1 - pp) * (tr$r22 * tr$c_to + tr$f23) + pp * tr$r33 * tr$d_to
    xs <- w * kappa * gs
    xd <- w * kappa * gd
    if (has5) {
      xs <- xs + w * k5 * t5$surv
      xd <- xd + w * k5 * t5$dens
    }
    m <- .sum_by(cbind(xs, xd), idx, ns)
    list(surv = m[, 1L], dens = m[, 2L])
  } else {
    x3 <- w * kappa * ((1 - pp) * tr$r23 + pp * tr$r33)
    if (has5) {
      x3 <- x3 + w * k5 * t5$surv
    }
    m <- .sum_by(cbind(w * kappa * (1 - pp) * tr$r22, x3), idx, ns)
    list(p12 = m[, 1L], p13 = m[, 2L])
  }
}
