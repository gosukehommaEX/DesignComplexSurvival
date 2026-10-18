#' Timing of the analyses
#'
#' Defines when the interim and final analyses take place: when the
#' expected number of events of an endpoint reaches given fractions of the
#' final number of events, or at given calendar times.
#'
#' @details
#' For event-driven looks (\code{trigger = "os"} or \code{"pfs"}), the final
#' number of events is given in \code{\link{ms_power}} through
#' \code{events} or \code{followup}, and the look k takes place when the
#' expected number of events of the trigger endpoint reaches
#' \code{fraction[k]} times the final number.
#'
#' @param trigger \code{"os"} or \code{"pfs"} for event-driven looks, or
#'   \code{"calendar"}.
#' @param fraction Strictly increasing event fractions ending at 1 (for
#'   event-driven looks; defaults to 1, a single analysis).
#' @param time Strictly increasing calendar times (for calendar-driven
#'   looks).
#'
#' @return An object of class \code{ms_looks}.
#'
#' @seealso \code{\link{ms_power}}, \code{\link{ms_size}}
#'
#' @examples
#' ms_looks("os", fraction = c(0.5, 0.75, 1))
#' ms_looks("calendar", time = c(24, 36))
#'
#' @export
ms_looks <- function(trigger = c("os", "pfs", "calendar"), fraction = NULL,
                     time = NULL) {
  trigger <- match.arg(trigger)
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
  structure(list(trigger = trigger, fraction = fraction, time = time, k = k),
            class = "ms_looks")
}
