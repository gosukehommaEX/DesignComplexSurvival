#' Final numbers of events of the trigger endpoints
#'
#' @param events Value of the argument \code{events} (checked by
#'   \code{.check_events()}).
#' @param eps Trigger endpoints of the event-driven looks.
#' @return A numeric vector named by \code{eps}.
#' @noRd
.match_events <- function(events, eps) {
  if (is.null(names(events))) {
    if (length(events) != 1L || length(eps) != 1L) {
      stop(paste0("'events' must be named by endpoint (\"os\", \"pfs\") when the ",
                  "looks are driven by the events of two endpoints."), call. = FALSE)
    }
    return(stats::setNames(as.numeric(events), eps))
  }
  if (!setequal(names(events), eps)) {
    stop(sprintf("'events' must give the final number of events of %s.",
                 paste(toupper(eps), collapse = " and ")), call. = FALSE)
  }
  stats::setNames(as.numeric(events[eps]), eps)
}
