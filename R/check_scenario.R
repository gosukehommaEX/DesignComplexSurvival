#' Check an ms_scenario object
#'
#' @param x Value to check.
#' @return \code{TRUE} invisibly; an error otherwise.
#' @noRd
.check_scenario <- function(x) {
  if (!inherits(x, "ms_scenario")) {
    stop("'scenario' must be an 'ms_scenario' object (see ms_scenario()).",
         call. = FALSE)
  }
  invisible(TRUE)
}
