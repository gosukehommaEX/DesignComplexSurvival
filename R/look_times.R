#' Calendar times of the looks
#'
#' An event-driven look takes place when the expected number of events of its
#' trigger endpoint reaches its fraction times the final number of events of
#' that endpoint. The final numbers of events are given by \code{events}; with
#' \code{followup} (event-driven final look) or without \code{events}
#' (calendar-driven final look), they are the expected numbers at the final
#' analysis, which takes place at \code{followup} after the end of accrual or
#' at its calendar time. With \code{targets}, the event-driven looks take
#' place when the expected numbers of events reach \code{targets} (used under
#' the null scenarios) and \code{events} and \code{followup} are not used.
#' The time of a final analysis fixed by \code{followup} is not searched for
#' from its expected number of events.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param acc Output of \code{.acc_setup()}.
#' @param lk Output of \code{.looks_table()}.
#' @param events Final numbers of events of the trigger endpoints
#'   (\code{.check_events()}), or \code{NULL}.
#' @param followup Time from the end of accrual to the final analysis, or
#'   \code{NULL}.
#' @param targets Event targets of the looks (\code{NA} for calendar-driven
#'   looks), or \code{NULL}.
#' @param what Name of the scenario for messages.
#' @return A list with \code{time}, \code{target} (event target of each look,
#'   \code{NA} for calendar-driven looks), \code{trigger} (endpoint whose
#'   observed events determine each time, \code{NA} for calendar times and for
#'   a final analysis fixed by \code{followup}) and \code{final} (final numbers
#'   of events of the trigger endpoints, named by endpoint).
#' @noRd
.look_times <- function(scenario, acc, lk, events = NULL, followup = NULL,
                        targets = NULL, what = "the scenario") {
  k <- nrow(lk)
  ev <- lk$trigger != "calendar"
  eps <- unique(lk$trigger[ev])
  times <- lk$time
  trig <- ifelse(ev, lk$trigger, NA_character_)
  final <- stats::setNames(numeric(0), character(0))
  t_fixed <- NULL
  if (is.null(targets)) {
    targets <- rep(NA_real_, k)
    if (!any(ev)) {
      if (!is.null(events) || !is.null(followup)) {
        stop("'events' and 'followup' are not used with calendar-driven looks.",
             call. = FALSE)
      }
    } else {
      t_ref <- NULL
      if (ev[k]) {
        if (is.null(events) == is.null(followup)) {
          stop("Give exactly one of 'events' and 'followup' for event-driven looks.",
               call. = FALSE)
        }
        if (!is.null(followup)) {
          t_ref <- acc$duration + followup
          if (!(t_ref > 0)) {
            stop("The final analysis must take place after the start of the study.",
                 call. = FALSE)
          }
          t_fixed <- t_ref
        }
      } else {
        if (!is.null(followup)) {
          stop("'followup' is not used when the final look is calendar-driven.",
               call. = FALSE)
        }
        if (is.null(events)) {
          t_ref <- lk$time[k]
        }
      }
      if (is.null(t_ref)) {
        final <- .match_events(events, eps)
      } else {
        final <- vapply(eps, function(e) {
          x <- .events_core(scenario, acc, t_ref, e)
          x$ev0 + x$ev1
        }, numeric(1))
      }
      targets[ev] <- lk$fraction[ev] * unname(final[lk$trigger[ev]])
    }
  }
  searched <- ev
  if (!is.null(t_fixed)) {
    searched[k] <- FALSE
    times[k] <- t_fixed
    trig[k] <- NA_character_
  }
  for (e in eps) {
    i <- which(lk$trigger == e & searched)
    if (length(i) > 0L) {
      times[i] <- .event_times(scenario, acc, targets[i], e)
    }
  }
  if (k > 1L && any(diff(times) <= 1e-6 * max(1, times))) {
    stop(sprintf("The analysis times under %s are not increasing (%s).", what,
                 paste(format(signif(times, 6)), collapse = ", ")), call. = FALSE)
  }
  list(time = times, target = targets, trigger = trig, final = final)
}
