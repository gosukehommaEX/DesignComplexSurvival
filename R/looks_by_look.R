#' Looks defined one by one
#'
#' The part of \code{ms_looks()} for one trigger per look.
#'
#' @param trigger,fraction,time,test See \code{ms_looks()}.
#' @return An object of class \code{ms_looks}.
#' @noRd
.looks_by_look <- function(trigger, fraction, time, test) {
  if (!is.character(trigger) || length(trigger) < 1L || anyNA(trigger) ||
      !all(trigger %in% c("os", "pfs", "calendar"))) {
    stop("'trigger' must contain \"os\", \"pfs\" or \"calendar\" (one value or one per look).",
         call. = FALSE)
  }
  k <- length(trigger)
  if (k > 20L) {
    stop("At most 20 looks are supported.", call. = FALSE)
  }
  ev <- trigger != "calendar"
  if (is.null(fraction)) {
    fraction <- rep(NA_real_, k)
  }
  if (is.null(time)) {
    time <- rep(NA_real_, k)
  }
  if (!(is.numeric(fraction) || all(is.na(fraction))) || length(fraction) != k ||
      any(is.na(fraction) != !ev)) {
    stop("'fraction' must have one value per look, NA exactly at the calendar-driven looks.",
         call. = FALSE)
  }
  if (!(is.numeric(time) || all(is.na(time))) || length(time) != k ||
      any(is.na(time) != ev)) {
    stop("'time' must have one value per look, NA exactly at the event-driven looks.",
         call. = FALSE)
  }
  fraction <- as.numeric(fraction)
  time <- as.numeric(time)
  if (any(ev) && (any(!is.finite(fraction[ev])) || any(fraction[ev] <= 0) ||
                  any(fraction[ev] > 1))) {
    stop("'fraction' must be in (0, 1] at the event-driven looks.", call. = FALSE)
  }
  for (e in c("os", "pfs")) {
    f <- fraction[trigger == e]
    if (length(f) > 1L && any(diff(f) <= 0)) {
      stop(sprintf("'fraction' must be strictly increasing over the looks driven by %s events.",
                   toupper(e)), call. = FALSE)
    }
  }
  if (ev[k] && abs(fraction[k] - 1) > 1e-12) {
    stop("'fraction' must be 1 at an event-driven final look.", call. = FALSE)
  }
  if (any(!ev)) {
    .check_times(time[!ev], "time", positive = TRUE)
    if (sum(!ev) > 1L && any(diff(time[!ev]) <= 0)) {
      stop("'time' must be strictly increasing.", call. = FALSE)
    }
  }
  structure(list(trigger = trigger, fraction = fraction, time = time, k = k,
                 test = .check_look_test(test, k)),
            class = "ms_looks")
}
