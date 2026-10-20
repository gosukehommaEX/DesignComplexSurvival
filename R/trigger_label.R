#' Description of the triggers of the looks
#'
#' @param trigger Triggers (\code{"os"}, \code{"pfs"} or \code{"calendar"}).
#' @return A character string such as \code{"OS events"} or
#'   \code{"PFS events and calendar time"}.
#' @noRd
.trigger_label <- function(trigger) {
  lab <- ifelse(trigger == "calendar", "calendar time", paste(toupper(trigger), "events"))
  paste(unique(lab), collapse = " and ")
}
