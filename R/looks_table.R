#' Per-look definition of the analyses
#'
#' Expands an \code{ms_looks} object into one row per look.
#'
#' @param looks An \code{ms_looks} object.
#' @return A data frame with \code{trigger} (\code{"os"}, \code{"pfs"} or
#'   \code{"calendar"}), \code{fraction} (\code{NA} for calendar-driven looks)
#'   and \code{time} (\code{NA} for event-driven looks).
#' @noRd
.looks_table <- function(looks) {
  k <- looks$k
  frac <- if (is.null(looks$fraction)) rep(NA_real_, k) else as.numeric(looks$fraction)
  tm <- if (is.null(looks$time)) rep(NA_real_, k) else as.numeric(looks$time)
  data.frame(trigger = rep_len(looks$trigger, k), fraction = frac, time = tm,
             stringsAsFactors = FALSE)
}
