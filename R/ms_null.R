#' Null scenario
#'
#' Replaces the model of the treatment group by the model of the control
#' group, keeping the allocation and the accrual.
#'
#' @param scenario An \code{ms_scenario} object.
#'
#' @return An \code{ms_scenario} object under the null hypothesis.
#'
#' @seealso \code{\link{ms_scenario}}
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' lr_moments(ms_null(scn), n = 400, calendar_time = 30)
#'
#' @export
ms_null <- function(scenario) {
  .check_scenario(scenario)
  scenario$treatment <- scenario$control
  scenario
}
