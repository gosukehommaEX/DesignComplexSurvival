#' Mean and variance of the log-rank statistic over calendar time
#'
#' Computes, at each calendar time, the expected events, the mean of the
#' log-rank score U (observed minus expected events in the treatment group),
#' its information V (the limit of the usual variance estimator) and the
#' mean and variance of the standardized statistic Z = -U / sqrt(V), which
#' is positive when the treatment is better.
#'
#' @details
#' Let y0 and y1 be the probabilities of being at risk at follow-up time s
#' in each group, multiplied by the allocation probabilities and by the
#' enrolled fraction at t - s. Then E(U) = n times the integral of
#' y0 y1 / (y0 + y1) (lambda1 - lambda0) and V = n times the integral of
#' y0 y1 / (y0 + y1)^2 (y0 lambda0 + y1 lambda1), from 0 to t.
#'
#' The mean of Z is
#' \itemize{
#'   \item \code{mean = "direct"}: -E(U) / sqrt(V), the mean of Yung and Liu
#'     (2020, Section 2.3.3);
#'   \item \code{mean = "schoenfeld"}: -sqrt(p0 p1 D) times the
#'     event-weighted average of the log hazard ratio, where D is the
#'     expected number of events (Schoenfeld 1981).
#' }
#' The variance of Z is 1 for \code{variance = "null"}. For
#' \code{variance = "alternative"} it is the asymptotic variance of the
#' score under the alternative divided by V, from Luo et al. (2019), which
#' accounts for the variability of the numbers at risk. This variance
#' assumes simple randomization; Yung and Liu (2020) show that block
#' randomization gives a variance that is not larger. For
#' \code{variance = "studentized"} it is the asymptotic variance of Z itself,
#' which also accounts for the variability of V and, for times determined by
#' the observed number of events of an endpoint (\code{trigger}), of the
#' analysis time (see \code{\link{lr_cov}}).
#'
#' @param scenario An \code{ms_scenario} object.
#' @param n Total sample size.
#' @param calendar_time Calendar times (positive).
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @param mean \code{"direct"} or \code{"schoenfeld"}.
#' @param variance \code{"null"}, \code{"alternative"} or
#'   \code{"studentized"}.
#' @param trigger Endpoint whose observed number of events determines each
#'   calendar time (\code{"os"} or \code{"pfs"}; one value or one per time),
#'   or \code{NULL} or \code{NA} for fixed calendar times. Used only with
#'   \code{variance = "studentized"}.
#'
#' @return A data frame with \code{time}, \code{events_control},
#'   \code{events_treatment}, \code{events}, \code{score_mean},
#'   \code{info}, \code{z_mean} and \code{z_var}.
#'
#' @references
#' Luo X, Mao X, Chen X, Qiu J, Bai S, Quan H (2019). Design and monitoring
#' of survival trials in complex scenarios. Statistics in Medicine
#' 38:192-209.
#'
#' Schoenfeld D (1981). The asymptotic properties of nonparametric tests for
#' comparing survival distributions. Biometrika 68:316-319.
#'
#' Yung G, Liu Y (2020). Sample size and power for the weighted log-rank
#' test and Kaplan-Meier based tests with allowance for nonproportional
#' hazards. Biometrics 76:939-950.
#'
#' @seealso \code{\link{lr_cov}}, \code{\link{ms_events}}
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' lr_moments(scn, n = 400, calendar_time = c(24, 36))
#'
#' @export
lr_moments <- function(scenario, n, calendar_time, endpoint = c("os", "pfs"),
                       mean = c("direct", "schoenfeld"),
                       variance = c("null", "alternative", "studentized"),
                       trigger = NULL) {
  endpoint <- match.arg(endpoint)
  mean <- match.arg(mean)
  variance <- match.arg(variance)
  .check_scenario(scenario)
  .check_positive(n, "n")
  .check_times(calendar_time, "calendar_time", positive = TRUE)
  trigger <- .check_trigger(trigger, length(calendar_time))
  acc <- .acc_setup(scenario$accrual, n)
  core <- .lr_core(scenario, acc, calendar_time, endpoint, mean)
  if (variance == "null") {
    z_var <- rep(1, length(calendar_time))
  } else if (variance == "alternative") {
    cv <- .lr_alt_cov(scenario, acc, calendar_time, endpoint, core$score_mean / n)
    z_var <- diag(cv) / (core$info / n)
  } else {
    z_var <- diag(.zeta_cov(scenario, acc, rep(endpoint, length(calendar_time)),
                            calendar_time, trigger, "studentized")$sigma)
  }
  data.frame(time = calendar_time, events_control = core$events_control,
             events_treatment = core$events_treatment, events = core$events,
             score_mean = core$score_mean, info = core$info,
             z_mean = core$z_mean, z_var = z_var)
}
