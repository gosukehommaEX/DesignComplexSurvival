#' Expected events and log-rank moments of OS when switching starts at a
#' calendar time
#'
#' The counterpart of \code{.lr_core()} (with \code{mean = "direct"}) for
#' scenarios with \code{.is_calendar()}: per calendar time, the integrals over
#' the follow-up time of the functions of \code{.cal_look()}.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param acc Output of \code{.acc_setup()}.
#' @param t Calendar times.
#' @return A data frame as \code{.lr_core()}.
#' @noRd
.cal_core <- function(scenario, acc, t) {
  n <- acc$n
  m <- matrix(0, length(t), 6L)
  for (i in seq_along(t)) {
    if (t[i] > 0) {
      lk <- .cal_look(scenario, acc, t[i])
      m[i, ] <- c(lk$ev, lk$do, lk$mu, lk$v)
    }
  }
  ev0 <- n * m[, 1L]
  ev1 <- n * m[, 2L]
  score <- n * m[, 5L]
  info <- n * m[, 6L]
  data.frame(time = t, events_control = ev0, events_treatment = ev1,
             events = ev0 + ev1, dropouts_control = n * m[, 3L],
             dropouts_treatment = n * m[, 4L], score_mean = score,
             info = info, z_mean = -score / sqrt(info))
}
