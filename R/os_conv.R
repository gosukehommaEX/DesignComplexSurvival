#' Contributions of progressed patients
#'
#' Integrates over the time of progression u from 0 to s. With
#' \code{what = "os"} returns the parts of the overall survival function and
#' density that come from patients who progressed by s; with
#' \code{what = "states"} returns the occupation probabilities of states 2
#' (progressed without switch) and 3 (switched).
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
  lam <- max(c(model$prog + model$death,
               model$switch_rate + model$death_prog, model$death_switch), 1e-12)
  forward <- model$clock == "forward"
  nodes <- lapply(s, function(si) {
    if (!(si > 0)) {
      return(list(u = numeric(0), w = numeric(0)))
    }
    brk <- if (forward) c(tm, tp) else c(tm, si - tp)
    g <- .gl_pieces(0, si, brk, hmax = max(2 / lam, si / 400), n = 10L)
    list(u = g$x, w = g$w)
  })
  len <- vapply(nodes, function(z) length(z$u), integer(1))
  idx <- rep.int(seq_len(ns), len)
  u <- as.numeric(unlist(lapply(nodes, `[[`, "u"), use.names = FALSE))
  w <- as.numeric(unlist(lapply(nodes, `[[`, "w"), use.names = FALSE))
  sv <- s[idx]
  kappa <- exp(-.pw_cumhaz(u, tm, model$prog + model$death)) *
    .pw_value(u, tm, model$prog)
  pp <- .pw_value(u, tm, model$switch_prob)
  tr <- if (forward) {
    .post_trans(model, u, sv)
  } else {
    .post_trans(model, numeric(length(u)), sv - u)
  }
  if (what == "os") {
    gs <- (1 - pp) * (tr$r22 + tr$r23) + pp * tr$r33
    gd <- (1 - pp) * (tr$r22 * tr$c_to + tr$r23 * tr$d_to) + pp * tr$r33 * tr$d_to
    m <- .sum_by(cbind(w * kappa * gs, w * kappa * gd), idx, ns)
    list(surv = m[, 1L], dens = m[, 2L])
  } else {
    m <- .sum_by(cbind(w * kappa * (1 - pp) * tr$r22,
                       w * kappa * ((1 - pp) * tr$r23 + pp * tr$r33)), idx, ns)
    list(p12 = m[, 1L], p13 = m[, 2L])
  }
}
