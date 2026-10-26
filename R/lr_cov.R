#' Joint distribution of the log-rank statistic at several calendar times
#'
#' Returns the mean vector and covariance matrix of the standardized
#' log-rank statistics Z(t_1), ..., Z(t_K) of one endpoint, or of PFS and OS
#' together.
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
#' With \code{variance = "studentized"} the covariance is the asymptotic
#' covariance of Z = -U / sqrt(V-hat) itself, where V-hat is the usual
#' variance estimator. The per-patient influence function of Z adds to that
#' of the score a term for the variability of V-hat and, for times determined
#' by the observed number of events of an endpoint (\code{trigger}), a term
#' for the variability of the analysis time. Both terms vanish under the null
#' hypothesis.
#'
#' For PFS and OS together (\code{endpoint = c("pfs", "os")}) the statistics
#' are ordered by endpoint and then by time. With \code{"studentized"} all
#' covariances come from the influence functions. With \code{"null"} the
#' statistics of one endpoint have the covariances described above and the
#' correlations between the endpoints are those of the studentized
#' statistics. With \code{"alternative"} the covariances between the
#' endpoints are the covariances of the scores divided by
#' \code{sqrt(V_P V_O)}. Under the null hypothesis the three choices agree,
#' and the covariance between the endpoints is that of Danzer et al. (2025,
#' Theorem 1).
#'
#' With switching from a calendar time (\code{switch_start} of
#' \code{\link{ms_scenario}}), the hazard of OS depends on the entry time
#' and the score has no independent increments, even under the null
#' hypothesis. The means, information and covariances of OS then integrate
#' over the entry time, with the proportion in the risk set and the pooled
#' hazard of each calendar time. With \code{"null"} the statistics have unit
#' variances and the correlations of the studentized statistics (which,
#' under the null hypothesis, are the exact correlations), with
#' \code{"alternative"} the covariances of the scores divided by
#' \code{sqrt(V_j V_k)}, and with \code{"studentized"} the covariances of
#' the studentized statistics; \code{trigger} is used with all three. Only
#' OS alone with \code{mean = "direct"} is covered; PFS is not affected.
#'
#' @inheritParams lr_moments
#' @param calendar_time Strictly increasing calendar times (positive).
#' @param endpoint \code{"os"}, \code{"pfs"} or both (for example
#'   \code{c("pfs", "os")}).
#'
#' @return A list with \code{time}, \code{endpoint}, \code{events},
#'   \code{info}, \code{mean} (mean vector of Z) and \code{sigma}
#'   (covariance matrix of Z), with one element (row) per endpoint and time.
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
#' Danzer MF, Rufibach K, Beyersmann J, Schmidt R (2025). Exhausting the type
#' I error level in event-driven group-sequential designs with a closed
#' testing procedure for progression-free and overall survival.
#' arXiv:2512.08658.
#'
#' @seealso \code{\link{lr_moments}}, \code{\link{gs_prob}}
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' lr_cov(scn, n = 400, calendar_time = c(24, 36))
#' fl <- ms_scenario(ms_fleischer(3, 9), ms_fleischer(4.3, 12.9),
#'                   accrual = ms_accrual(rate = 20))
#' lr_cov(fl, n = 400, calendar_time = 30, endpoint = c("pfs", "os"),
#'        variance = "studentized")
#'
#' @export
lr_cov <- function(scenario, n, calendar_time, endpoint = "os",
                   mean = c("direct", "schoenfeld"),
                   variance = c("null", "alternative", "studentized"),
                   trigger = NULL) {
  mean <- match.arg(mean)
  variance <- match.arg(variance)
  if (!is.character(endpoint) || !(length(endpoint) %in% 1:2) || anyNA(endpoint) ||
      !all(endpoint %in% c("os", "pfs")) || anyDuplicated(endpoint) > 0L) {
    stop("'endpoint' must be \"os\", \"pfs\" or both.", call. = FALSE)
  }
  .check_scenario(scenario)
  .check_positive(n, "n")
  .check_times(calendar_time, "calendar_time", positive = TRUE)
  if (length(calendar_time) > 1L && any(diff(calendar_time) <= 0)) {
    stop("'calendar_time' must be strictly increasing.", call. = FALSE)
  }
  k <- length(calendar_time)
  trigger <- .check_trigger(trigger, k)
  .check_calendar(scenario, endpoint, mean, variance, trigger)
  acc <- .acc_setup(scenario$accrual, n)
  if (length(endpoint) == 1L) {
    core <- .lr_core(scenario, acc, calendar_time, endpoint, mean)
    sigma <- .z_cov(scenario, acc, calendar_time, endpoint, core, variance, trigger)
    return(list(time = calendar_time, endpoint = rep(endpoint, k),
                events = core$events, info = core$info, mean = core$z_mean,
                sigma = sigma))
  }
  cores <- lapply(endpoint, function(e) .lr_core(scenario, acc, calendar_time, e, mean))
  zc <- .zeta_cov(scenario, acc, rep(endpoint, each = k), rep(calendar_time, 2L),
                  if (is.null(trigger)) NULL else rep(trigger, 2L), variance)
  list(time = rep(calendar_time, 2L), endpoint = rep(endpoint, each = k),
       events = unlist(lapply(cores, `[[`, "events"), use.names = FALSE),
       info = unlist(lapply(cores, `[[`, "info"), use.names = FALSE),
       mean = unlist(lapply(cores, `[[`, "z_mean"), use.names = FALSE),
       sigma = zc$sigma)
}
