#' Whether a scenario switches from a calendar time
#'
#' The calendar-time computations (\code{.cal_core()}, \code{.cal_zeta()}) are
#' used for OS when \code{switch_start} of the scenario is positive and at
#' least one group can switch after progression. Scenarios created before
#' version 0.3.0 have no \code{switch_start} and use the other computations.
#'
#' @param scenario An \code{ms_scenario} object.
#' @return \code{TRUE} or \code{FALSE}.
#' @noRd
.is_calendar <- function(scenario) {
  tc <- scenario$switch_start
  if (is.null(tc) || !(tc > 0)) {
    return(FALSE)
  }
  .has_switch(scenario$control) || .has_switch(scenario$treatment)
}
