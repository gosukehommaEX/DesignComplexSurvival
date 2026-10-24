#' Survival after progression from the switched progression-free state
#'
#' Survival and death density between the post-progression clock values
#' \code{from} (progression) and \code{to} of a patient who progressed after a
#' switch before progression (from state 5 to state 3). The death hazard is
#' \code{death_switch} on the post-progression clock, so the survival is
#' exp(-(D(to) - D(from))) with D the cumulative hazard of
#' \code{death_switch}.
#'
#' @param model An \code{ms_model} object.
#' @param from,to Numeric vectors of equal length with \code{from <= to}.
#' @return A list with \code{surv} and \code{dens}.
#' @noRd
.post_switched <- function(model, from, to) {
  tp <- model$time_post
  d <- model$death_switch
  surv <- exp(-(.pw_cumhaz(to, tp, d) - .pw_cumhaz(from, tp, d)))
  list(surv = surv, dens = surv * .pw_value(to, tp, d))
}
