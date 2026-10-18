#' Number of events, sample size, accrual duration or follow-up for a target power
#'
#' Solves for one design parameter so that the design of
#' \code{\link{ms_power}} has the target power, with the other parameters
#' fixed.
#'
#' @details
#' The unknown is one of
#' \itemize{
#'   \item \code{"events"}: the final number of events of the trigger
#'     endpoint, with \code{n} and the accrual fixed (event-driven looks);
#'   \item \code{"followup"}: the time from the end of accrual to the final
#'     analysis, with \code{n} and the accrual fixed (event-driven looks);
#'   \item \code{"n"}: the sample size, with the accrual rates (or, for an
#'     accrual given by proportions, the accrual duration) fixed and the
#'     final analysis at \code{followup} after the end of accrual
#'     (event-driven looks) or at the calendar times of the looks;
#'   \item \code{"accrual_duration"}: the accrual duration of an accrual
#'     given by proportions, with \code{n} fixed and the final analysis at
#'     \code{followup} after the end of accrual (event-driven looks) or at
#'     the calendar times of the looks.
#' }
#' For \code{"events"} and \code{"followup"} the search is over the calendar
#' time of the final analysis, so that the two give the same design (before
#' rounding). The search doubles or halves the unknown until the power
#' brackets the target, evaluates the power on a grid within the bracket
#' (or within \code{interval}) and refines the first grid interval in which
#' the target is reached, so that the smallest solution is returned when the
#' power is not monotone. Beta spending uses one minus \code{power} as the
#' total type II error.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param looks An \code{ms_looks} object.
#' @param test An \code{ms_test} object.
#' @param power Target power.
#' @param solve_for \code{"events"}, \code{"n"}, \code{"accrual_duration"}
#'   or \code{"followup"}.
#' @param n Total sample size (fixed for \code{"events"}, \code{"followup"}
#'   and \code{"accrual_duration"}).
#' @param followup Time from the end of accrual to the final analysis
#'   (fixed for \code{"n"} and \code{"accrual_duration"} with event-driven
#'   looks).
#' @param interval Optional search interval for the unknown (for
#'   \code{"events"} and \code{"followup"}, an interval of the calendar time
#'   of the final analysis).
#' @param round_up Whether to round the number of events or the sample size
#'   up to an integer.
#'
#' @return An object of class \code{c("ms_size", "ms_power")}: the result of
#'   \code{\link{ms_power}} for the solved design, with \code{solve_for},
#'   the unrounded \code{solution} and the \code{target_power}.
#'
#' @seealso \code{\link{ms_power}}
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' ms_size(scn, looks = ms_looks("os", fraction = c(0.5, 1)),
#'         test = ms_test(), power = 0.9, solve_for = "events", n = 500)
#'
#' @export
ms_size <- function(scenario, looks, test, power = 0.9,
                    solve_for = c("events", "n", "accrual_duration", "followup"),
                    n = NULL, followup = NULL, interval = NULL,
                    round_up = TRUE) {
  solve_for <- match.arg(solve_for)
  .check_scenario(scenario)
  if (!inherits(looks, "ms_looks")) {
    stop("'looks' must be an 'ms_looks' object.", call. = FALSE)
  }
  if (!inherits(test, "ms_test")) {
    stop("'test' must be an 'ms_test' object.", call. = FALSE)
  }
  .check_prob(power, "power")
  if (power <= test$alpha) {
    stop("'power' must be larger than the type I error.", call. = FALSE)
  }
  beta <- 1 - power
  event_driven <- looks$trigger != "calendar"
  if (solve_for %in% c("events", "followup")) {
    if (!event_driven) {
      stop("solve_for = \"", solve_for, "\" needs event-driven looks.", call. = FALSE)
    }
    if (is.null(n)) {
      stop("'n' must be given.", call. = FALSE)
    }
    if (!is.null(followup)) {
      stop("'followup' must not be given when solving for the events or the follow-up.",
           call. = FALSE)
    }
    .check_positive(n, "n")
    a_dur <- .acc_setup(scenario$accrual, n)$duration
    fun <- function(x) {
      .design_core(scenario, n, looks, test, followup = x - a_dur, beta = beta,
                   type1 = FALSE)$power
    }
    t_final <- .solve_power(fun, power, start = max(a_dur, 1), interval = interval,
                            what = "the time of the final analysis")
    if (solve_for == "events") {
      acc <- .acc_setup(scenario$accrual, n)
      ev <- .events_core(scenario, acc, t_final, looks$trigger)
      solution <- ev$ev0 + ev$ev1
      d_final <- if (round_up) ceiling(solution - 1e-8) else solution
      res <- .design_core(scenario, n, looks, test, events = d_final, beta = beta)
    } else {
      solution <- t_final - a_dur
      if (solution < 0) {
        stop("The target power is reached before the end of accrual; ",
             "use solve_for = \"events\".", call. = FALSE)
      }
      res <- .design_core(scenario, n, looks, test, followup = solution, beta = beta)
    }
  } else {
    if (event_driven && is.null(followup)) {
      stop("'followup' must be given for event-driven looks.", call. = FALSE)
    }
    if (!event_driven && !is.null(followup)) {
      stop("'followup' is not used with calendar-driven looks.", call. = FALSE)
    }
    if (!is.null(followup) && (!is.numeric(followup) || length(followup) != 1L ||
                               is.na(followup) || followup < 0)) {
      stop("'followup' must be a nonnegative number.", call. = FALSE)
    }
    if (solve_for == "n") {
      fun <- function(x) {
        .design_core(scenario, x, looks, test, followup = followup, beta = beta,
                     type1 = FALSE)$power
      }
      solution <- .solve_power(fun, power, start = 100, interval = interval,
                               what = "n")
      n_final <- if (round_up) ceiling(solution - 1e-8) else solution
      res <- .design_core(scenario, n_final, looks, test, followup = followup,
                          beta = beta)
    } else {
      if (scenario$accrual$type != "prop") {
        stop("solve_for = \"accrual_duration\" needs an accrual given by 'prop' and 'duration'.",
             call. = FALSE)
      }
      if (is.null(n)) {
        stop("'n' must be given.", call. = FALSE)
      }
      .check_positive(n, "n")
      fun <- function(x) {
        .design_core(.scale_accrual(scenario, x), n, looks, test,
                     followup = followup, beta = beta, type1 = FALSE)$power
      }
      solution <- .solve_power(fun, power, start = scenario$accrual$duration,
                               interval = interval, what = "the accrual duration")
      res <- .design_core(.scale_accrual(scenario, solution), n, looks, test,
                          followup = followup, beta = beta)
    }
  }
  res$solve_for <- solve_for
  res$solution <- solution
  res$target_power <- power
  structure(res, class = c("ms_size", "ms_power"))
}
