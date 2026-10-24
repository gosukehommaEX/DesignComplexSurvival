#' Expected numbers of patients, events and dropouts over calendar time
#'
#' @details
#' For group g with allocation probability p_g, the expected number of
#' events by calendar time t is \code{n p_g} times the integral from 0 to t
#' of A(t - s) f_g(s) G_g(s), where A is the enrolled fraction, f_g the
#' density of the endpoint and G_g the probability of no dropout. The
#' expected number of dropouts is obtained in the same way with the dropout
#' density. For PFS, the dropouts include the censoring of PFS only
#' (\code{dropout_pfs} of \code{\link{ms_model}}) and, with
#' \code{pfs_switch = "censor"} in \code{\link{ms_scenario}}, the censoring
#' at a switch before progression.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param n Total sample size.
#' @param calendar_time Calendar times from the start of the study.
#' @param endpoint \code{"os"} or \code{"pfs"}.
#'
#' @return A data frame with \code{time}, the numbers enrolled
#'   (\code{enrolled}, \code{enrolled_control}, \code{enrolled_treatment}),
#'   the expected events (\code{events_control}, \code{events_treatment},
#'   \code{events}) and the expected dropouts (\code{dropouts_control},
#'   \code{dropouts_treatment}).
#'
#' @seealso \code{\link{ms_event_time}}, \code{\link{lr_moments}}
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12, dropout = 0.005)
#' trt <- ms_model(death = 0.7 * log(2) / 12, dropout = 0.005)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' ms_events(scn, n = 400, calendar_time = c(12, 24, 36))
#'
#' @export
ms_events <- function(scenario, n, calendar_time, endpoint = c("os", "pfs")) {
  endpoint <- match.arg(endpoint)
  .check_scenario(scenario)
  .check_positive(n, "n")
  .check_times(calendar_time, "calendar_time")
  acc <- .acc_setup(scenario$accrual, n)
  ev <- .events_core(scenario, acc, calendar_time, endpoint)
  enr <- n * .acc_frac(acc, calendar_time)
  p <- scenario$prob
  data.frame(time = calendar_time, enrolled = enr,
             enrolled_control = enr * p[[1L]], enrolled_treatment = enr * p[[2L]],
             events_control = ev$ev0, events_treatment = ev$ev1,
             events = ev$ev0 + ev$ev1,
             dropouts_control = ev$do0, dropouts_treatment = ev$do1)
}
