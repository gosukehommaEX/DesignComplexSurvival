#' Joint distribution of the log-rank statistic at several calendar times
#'
#' Returns the mean vector and covariance matrix of the standardized
#' log-rank statistics Z(t_1), ..., Z(t_K) of one endpoint.
#'
#' @details
#' With \code{variance = "null"} the statistics have unit variances and
#' correlations \code{sqrt(V(t_j) / V(t_k))} for t_j <= t_k, from the
#' independent increments of the score (Tsiatis 1981). With
#' \code{variance = "alternative"} the covariance of the score under the
#' alternative is that of Luo et al. (2019), computed here from the
#' per-patient influence functions of the score, and is divided by
#' \code{sqrt(V(t_j) V(t_k))}.
#'
#' @inheritParams lr_moments
#' @param calendar_time Strictly increasing calendar times (positive).
#'
#' @return A list with \code{time}, \code{events}, \code{info},
#'   \code{mean} (mean vector of Z) and \code{sigma} (covariance matrix of
#'   Z).
#'
#' @references
#' Tsiatis AA (1981). The asymptotic joint distribution of the efficient
#' scores test for the proportional hazards model calculated over time.
#' Biometrika 68:311-315.
#'
#' Luo X, Mao X, Chen X, Qiu J, Bai S, Quan H (2019). Design and monitoring
#' of survival trials in complex scenarios. Statistics in Medicine
#' 38:192-209.
#'
#' @seealso \code{\link{lr_moments}}, \code{\link{gs_prob}}
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' lr_cov(scn, n = 400, calendar_time = c(24, 36))
#'
#' @export
lr_cov <- function(scenario, n, calendar_time, endpoint = c("os", "pfs"),
                   mean = c("direct", "schoenfeld"),
                   variance = c("null", "alternative")) {
  endpoint <- match.arg(endpoint)
  mean <- match.arg(mean)
  variance <- match.arg(variance)
  .check_scenario(scenario)
  .check_positive(n, "n")
  .check_times(calendar_time, "calendar_time", positive = TRUE)
  if (length(calendar_time) > 1L && any(diff(calendar_time) <= 0)) {
    stop("'calendar_time' must be strictly increasing.", call. = FALSE)
  }
  acc <- .acc_setup(scenario$accrual, n)
  core <- .lr_core(scenario, acc, calendar_time, endpoint, mean)
  sigma <- .z_cov(scenario, acc, calendar_time, endpoint, core, variance)
  list(time = calendar_time, events = core$events, info = core$info,
       mean = core$z_mean, sigma = sigma)
}
