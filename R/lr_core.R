#' Expected events and log-rank moments at calendar times
#'
#' @param scenario An \code{ms_scenario} object.
#' @param acc Output of \code{.acc_setup()}.
#' @param t Calendar times.
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @param mean \code{"direct"} or \code{"schoenfeld"}.
#' @return A data frame with the expected events and dropouts by group, the
#'   score mean (observed minus expected events in the treatment group), the
#'   information and the mean of the standardized statistic
#'   Z = -U / sqrt(V), which is positive when the treatment is better.
#'   With switching from a calendar time (\code{.is_calendar()}), OS is
#'   computed by \code{.cal_core()}.
#' @noRd
.lr_core <- function(scenario, acc, t, endpoint, mean = "direct") {
  if (endpoint == "os" && .is_calendar(scenario)) {
    if (mean != "direct") {
      stop("Switching from a calendar time ('switch_start') needs mean = \"direct\" for OS.",
           call. = FALSE)
    }
    return(.cal_core(scenario, acc, t))
  }
  sch <- mean == "schoenfeld"
  hfun <- function(s) .lr_integrands(scenario, s, endpoint, schoenfeld = sch)
  res <- .acc_integral(acc, t, hfun, breaks = .scn_breaks(scenario),
                       hmax = .scn_hmax(scenario, max(c(t, 1))))
  n <- acc$n
  ev0 <- n * res[, "ev0"]
  ev1 <- n * res[, "ev1"]
  score <- n * res[, "u"]
  info <- n * res[, "v"]
  if (sch) {
    p <- scenario$prob
    z <- -sqrt(p[[1L]] * p[[2L]]) * (n * res[, "lw"]) / sqrt(ev0 + ev1)
  } else {
    z <- -score / sqrt(info)
  }
  data.frame(time = t, events_control = ev0, events_treatment = ev1,
             events = ev0 + ev1, dropouts_control = n * res[, "do0"],
             dropouts_treatment = n * res[, "do1"], score_mean = score,
             info = info, z_mean = z)
}
