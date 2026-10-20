#' Check the final numbers of events
#'
#' @param events A positive number, or positive numbers named by endpoint
#'   (\code{"os"}, \code{"pfs"}).
#' @return \code{TRUE} invisibly; an error otherwise.
#' @noRd
.check_events <- function(events) {
  if (is.null(names(events))) {
    return(.check_positive(events, "events"))
  }
  if (!is.numeric(events) || !(length(events) %in% 1:2) || anyNA(events) ||
      any(!is.finite(events)) || any(events <= 0) ||
      !all(names(events) %in% c("os", "pfs")) || anyDuplicated(names(events)) > 0L) {
    stop("'events' must be a positive number or positive numbers named \"os\" and \"pfs\".",
         call. = FALSE)
  }
  invisible(TRUE)
}
