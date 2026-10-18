#' Scenario with the accrual duration changed
#'
#' Scales the interval limits of an accrual given by proportions so that
#' the accrual duration becomes \code{duration}.
#'
#' @param scenario An \code{ms_scenario} object with \code{prop} accrual.
#' @param duration New accrual duration.
#' @return An \code{ms_scenario} object.
#' @noRd
.scale_accrual <- function(scenario, duration) {
  acc <- scenario$accrual
  f <- duration / acc$duration
  acc$time <- acc$time * f
  acc$duration <- duration
  scenario$accrual <- acc
  scenario
}
