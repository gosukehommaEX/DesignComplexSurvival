#' Check the argument trigger
#'
#' @param trigger \code{NULL}, or a vector of length 1 or \code{k} with
#'   \code{"os"}, \code{"pfs"} or \code{NA} (a fixed calendar time).
#' @param k Number of calendar times.
#' @return \code{NULL} or a character vector of length \code{k}.
#' @noRd
.check_trigger <- function(trigger, k) {
  if (is.null(trigger)) {
    return(NULL)
  }
  if (!(is.character(trigger) || all(is.na(trigger))) ||
      !(length(trigger) %in% c(1L, k)) ||
      !all(is.na(trigger) | trigger %in% c("os", "pfs"))) {
    stop("'trigger' must be NULL or \"os\", \"pfs\" or NA (one value or one per time).",
         call. = FALSE)
  }
  rep_len(as.character(trigger), k)
}
