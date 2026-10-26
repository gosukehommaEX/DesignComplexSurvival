#' Check the scope of the calendar-time switching
#'
#' With switching from a calendar time (\code{.is_calendar()}), the
#' distribution of OS depends on the entry time. The package then covers the
#' log-rank statistic of OS alone with \code{mean = "direct"}; the joint
#' distribution with PFS (two endpoints, an OS test at looks driven by PFS
#' events, or a PFS test at looks driven by OS events with a variance other
#' than \code{"null"}) is not covered. PFS alone is not affected by the
#' switching and is computed as before.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param endpoints Endpoints of the test or of the statistic.
#' @param mean \code{"direct"} or \code{"schoenfeld"}.
#' @param variance \code{"null"}, \code{"alternative"} or
#'   \code{"studentized"}.
#' @param trigger Triggers of the looks (\code{"os"}, \code{"pfs"},
#'   \code{"calendar"} or \code{NA}), or \code{NULL}.
#' @return \code{TRUE}, invisibly; an error outside the scope.
#' @noRd
.check_calendar <- function(scenario, endpoints, mean = "direct", variance = "null",
                            trigger = NULL) {
  if (!.is_calendar(scenario)) {
    return(invisible(TRUE))
  }
  trig <- unique(trigger[!is.na(trigger) & trigger != "calendar"])
  if (length(endpoints) > 1L) {
    stop("Switching from a calendar time ('switch_start') supports one endpoint only.",
         call. = FALSE)
  }
  if (endpoints == "os") {
    if (mean != "direct") {
      stop("Switching from a calendar time ('switch_start') needs mean = \"direct\" for OS.",
           call. = FALSE)
    }
    if ("pfs" %in% trig) {
      stop(paste0("With switching from a calendar time ('switch_start'), the looks of an ",
                  "OS test cannot be driven by PFS events."), call. = FALSE)
    }
  } else if ("os" %in% trig && variance != "null") {
    stop(paste0("With switching from a calendar time ('switch_start'), a PFS test at looks ",
                "driven by OS events needs variance = \"null\"."), call. = FALSE)
  }
  invisible(TRUE)
}
