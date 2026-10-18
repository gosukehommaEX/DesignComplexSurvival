#' Quantities of one group used in the log-rank calculations
#'
#' @param model An \code{ms_model} object.
#' @param s Follow-up times since randomization.
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @return A list with the survival function \code{S}, density \code{f},
#'   probability of being at risk \code{R} (no event and no dropout), density
#'   of an observed event \code{e} and density of an observed dropout
#'   \code{c}.
#' @noRd
.grp_eval <- function(model, s, endpoint) {
  ep <- .ep_eval(model, s, endpoint)
  g <- exp(-.pw_cumhaz(s, model$time, model$dropout))
  r <- ep$surv * g
  list(S = ep$surv, f = ep$dens, R = r, e = ep$dens * g,
       c = r * .pw_value(pmax(s, 0), model$time, model$dropout))
}
