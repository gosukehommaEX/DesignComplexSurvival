#' Calendar times of given numbers of OS events when switching starts at a
#' calendar time
#'
#' Newton's method on the expected number of OS events of \code{.cal_core()},
#' with its derivative with respect to the calendar time from
#' \code{.cal_dt()}, started at the time of the same number of events when
#' switching is possible from randomization (\code{switch_start = 0}, the
#' computations of \code{.event_times()}). A step that leaves the interval
#' known to contain the time is replaced by bisection (or by doubling the
#' time while no upper end is known). The iteration stops when a step is at
#' most 1e-11 max(t, 1) and returns the time after that step.
#'
#' @param scenario An \code{ms_scenario} object with \code{.is_calendar()}.
#' @param acc Output of \code{.acc_setup()}.
#' @param targets Expected numbers of OS events.
#' @return Calendar times; \code{NA} where no starting time is available
#'   (the target cannot be reached with \code{switch_start = 0}) or
#'   the iteration does not converge in 30 steps (\code{.event_times()} then
#'   uses its bracketing search).
#' @noRd
.cal_event_times <- function(scenario, acc, targets) {
  s0 <- scenario
  s0$switch_start <- 0
  start <- vapply(targets, function(d) {
    tryCatch(.event_times(s0, acc, d, "os"), error = function(e) NA_real_)
  }, numeric(1))
  n <- acc$n
  out <- rep(NA_real_, length(targets))
  for (i in seq_along(targets)) {
    t <- start[i]
    if (!is.finite(t) || !(t > 0)) {
      next
    }
    lo <- 0
    hi <- Inf
    for (it in seq_len(30L)) {
      lk <- .cal_look(scenario, acc, t)
      dd <- .cal_dt(scenario, acc, lk)$dd
      f <- n * sum(lk$ev) - targets[i]
      if (f == 0) {
        out[i] <- t
        break
      }
      if (f < 0) {
        lo <- t
      } else {
        hi <- t
      }
      der <- n * sum(lk$g$w * (dd[[1L]] + dd[[2L]]))
      t_new <- if (der > 0) t - f / der else NA_real_
      if (is.finite(t_new) && abs(t_new - t) <= 1e-11 * max(t, 1)) {
        out[i] <- t_new
        break
      }
      if (!is.finite(t_new) || t_new <= lo || t_new >= hi) {
        t_new <- if (is.finite(hi)) (lo + hi) / 2 else 2 * t
      }
      t <- t_new
    }
  }
  out
}
