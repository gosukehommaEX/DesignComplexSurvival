#' Timing of the analyses
#'
#' Defines when the interim and final analyses take place: when the
#' expected number of events of an endpoint reaches given fractions of its
#' final number of events, or at given calendar times, and which endpoints
#' are tested at each analysis.
#'
#' @details
#' With one value of \code{trigger}, all looks are driven by the events of
#' the same endpoint (\code{"os"} or \code{"pfs"}, with \code{fraction}) or
#' take place at calendar times (\code{"calendar"}, with \code{time}). The
#' final number of events is given in \code{\link{ms_power}} through
#' \code{events} or \code{followup}, and the look k takes place when the
#' expected number of events of the trigger endpoint reaches
#' \code{fraction[k]} times the final number.
#'
#' With one value of \code{trigger} per look, each look has its own trigger.
#' \code{fraction} then has one value per look, \code{NA} at the
#' calendar-driven looks, and \code{time} one value per look, \code{NA} at the
#' event-driven looks. The fractions of the looks driven by the events of an
#' endpoint refer to the final number of events of that endpoint (see
#' \code{\link{ms_power}}) and must be strictly increasing; an event-driven
#' final look must have the fraction 1. The calendar times of the looks, which
#' depend on the scenario, must be strictly increasing.
#'
#' By default every look tests all endpoints of the test
#' (\code{\link{ms_test}}). \code{test} gives the endpoints tested at each
#' look, for example \code{list("pfs", c("pfs", "os"), "os")} for an interim
#' analysis of PFS, the final analysis of PFS with an interim analysis of OS,
#' and the final analysis of OS.
#'
#' @param trigger \code{"os"} or \code{"pfs"} for event-driven looks, or
#'   \code{"calendar"}; one value for all looks or one value per look.
#' @param fraction Event fractions of the event-driven looks: with one
#'   trigger, strictly increasing and ending at 1 (defaults to 1, a single
#'   analysis); with one trigger per look, one value per look in (0, 1]
#'   (\code{NA} at calendar-driven looks).
#' @param time Strictly increasing calendar times of the calendar-driven
#'   looks (with one trigger per look, one value per look, \code{NA} at
#'   event-driven looks).
#' @param test \code{NULL} (all endpoints at all looks), or a list with one
#'   element per look giving the endpoints tested at the look (\code{"os"},
#'   \code{"pfs"} or both). A character vector with one endpoint per look is
#'   also accepted.
#'
#' @return An object of class \code{ms_looks}.
#'
#' @seealso \code{\link{ms_power}}, \code{\link{ms_size}}
#'
#' @examples
#' ms_looks("os", fraction = c(0.5, 0.75, 1))
#' ms_looks("calendar", time = c(24, 36))
#' # PFS interim and final analyses driven by PFS events, an OS interim
#' # analysis at the final PFS analysis and the final OS analysis driven by
#' # OS events
#' ms_looks(c("pfs", "pfs", "os"), fraction = c(0.6, 1, 1),
#'          test = list("pfs", c("pfs", "os"), "os"))
#'
#' @export
ms_looks <- function(trigger = c("os", "pfs", "calendar"), fraction = NULL,
                     time = NULL, test = NULL) {
  if (missing(trigger)) {
    trigger <- "os"
  }
  if (length(trigger) != 1L) {
    return(.looks_by_look(trigger, fraction, time, test))
  }
  trigger <- match.arg(trigger, c("os", "pfs", "calendar"))
  if (trigger == "calendar") {
    if (is.null(time) || !is.null(fraction)) {
      stop("Calendar-driven looks need 'time' and no 'fraction'.", call. = FALSE)
    }
    .check_times(time, "time", positive = TRUE)
    if (length(time) > 1L && any(diff(time) <= 0)) {
      stop("'time' must be strictly increasing.", call. = FALSE)
    }
    k <- length(time)
    time <- as.numeric(time)
  } else {
    if (!is.null(time)) {
      stop("Event-driven looks need 'fraction' and no 'time'.", call. = FALSE)
    }
    if (is.null(fraction)) {
      fraction <- 1
    }
    if (!is.numeric(fraction) || anyNA(fraction) || any(fraction <= 0) ||
        (length(fraction) > 1L && any(diff(fraction) <= 0)) ||
        abs(fraction[length(fraction)] - 1) > 1e-12) {
      stop("'fraction' must be strictly increasing, positive and end at 1.",
           call. = FALSE)
    }
    k <- length(fraction)
    fraction <- as.numeric(fraction)
  }
  if (k > 20L) {
    stop("At most 20 looks are supported.", call. = FALSE)
  }
  structure(list(trigger = trigger, fraction = fraction, time = time, k = k,
                 test = .check_look_test(test, k)),
            class = "ms_looks")
}
