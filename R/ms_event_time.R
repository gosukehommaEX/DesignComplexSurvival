#' Calendar time at which a number of events is expected
#'
#' Solves for the calendar time at which the expected total number of
#' events of an endpoint reaches each target.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param n Total sample size.
#' @param events Target numbers of events (positive).
#' @param endpoint \code{"os"} or \code{"pfs"}.
#'
#' @return Numeric vector of calendar times. An error is raised when a
#'   target cannot be reached.
#'
#' @seealso \code{\link{ms_events}}
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' ms_event_time(scn, n = 400, events = c(150, 300))
#'
#' @export
ms_event_time <- function(scenario, n, events, endpoint = c("os", "pfs")) {
  endpoint <- match.arg(endpoint)
  .check_scenario(scenario)
  .check_positive(n, "n")
  .check_times(events, "events", positive = TRUE)
  acc <- .acc_setup(scenario$accrual, n)
  .event_times(scenario, acc, events, endpoint)
}
