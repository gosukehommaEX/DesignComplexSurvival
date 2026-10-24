#' Survival function and density of an endpoint for one group
#'
#' With a switch before progression (\code{switch_pre} positive), PFS is the
#' time to leave the progression-free states 1 and 5 by progression or death
#' when \code{model$pfs_switch} is \code{"event"}. When it is
#' \code{"censor"}, PFS is censored at the switch; \code{surv} is then the
#' probability of being in state 1 and \code{dens} the density of leaving
#' state 1 by progression or death, so that their ratio is the hazard
#' \code{prog + death} compared by the log-rank test.
#'
#' @param model An \code{ms_model} object.
#' @param s Follow-up times since randomization (nonnegative).
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @return A list with \code{surv} and \code{dens} (without dropout).
#' @noRd
.ep_eval <- function(model, s, endpoint) {
  s <- pmax(s, 0)
  tm <- model$time
  if (any(model$switch_pre > 0)) {
    pre <- .pre_states(model, s)
    if (endpoint == "pfs") {
      e1 <- pre$p1 * .pw_value(s, tm, model$prog + model$death)
      if (identical(model$pfs_switch, "censor")) {
        return(list(surv = pre$p1, dens = e1))
      }
      e5 <- pre$p5 * .pw_value(s, tm, model$prog_switch + model$death_switch_pre)
      return(list(surv = pre$p1 + pre$p5, dens = e1 + e5))
    }
    surv <- pre$p1 + pre$p5
    dens <- pre$p1 * .pw_value(s, tm, model$death) +
      pre$p5 * .pw_value(s, tm, model$death_switch_pre)
    if (any(model$prog > 0) || any(model$prog_switch > 0)) {
      conv <- .os_conv(model, s, what = "os")
      surv <- surv + conv$surv
      dens <- dens + conv$dens
    }
    return(list(surv = surv, dens = dens))
  }
  p11 <- exp(-.pw_cumhaz(s, tm, model$prog + model$death))
  if (endpoint == "pfs") {
    return(list(surv = p11, dens = .pw_value(s, tm, model$prog + model$death) * p11))
  }
  b_s <- .pw_value(s, tm, model$death)
  if (all(model$prog == 0)) {
    return(list(surv = p11, dens = b_s * p11))
  }
  conv <- .os_conv(model, s, what = "os")
  list(surv = p11 + conv$surv, dens = b_s * p11 + conv$dens)
}
