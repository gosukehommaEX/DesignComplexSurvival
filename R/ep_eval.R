#' Survival function and density of an endpoint for one group
#'
#' @param model An \code{ms_model} object.
#' @param s Follow-up times since randomization (nonnegative).
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @return A list with \code{surv} and \code{dens} (without dropout).
#' @noRd
.ep_eval <- function(model, s, endpoint) {
  s <- pmax(s, 0)
  tm <- model$time
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
