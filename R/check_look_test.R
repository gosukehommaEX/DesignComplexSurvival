#' Check the endpoints tested at each look
#'
#' @param test \code{NULL}, a list with one element per look (\code{"os"},
#'   \code{"pfs"} or both), or a character vector with one endpoint per look.
#' @param k Number of looks.
#' @return \code{NULL} or a list of \code{k} character vectors.
#' @noRd
.check_look_test <- function(test, k) {
  if (is.null(test)) {
    return(NULL)
  }
  if (is.character(test) && length(test) == k && !anyNA(test)) {
    test <- as.list(test)
  }
  ok <- is.list(test) && length(test) == k &&
    all(vapply(test, function(x) {
      is.character(x) && length(x) %in% 1:2 && !anyNA(x) &&
        all(x %in% c("os", "pfs")) && anyDuplicated(x) == 0L
    }, logical(1)))
  if (!ok) {
    stop(paste0("'test' must be NULL or a list with one element per look, each ",
                "\"os\", \"pfs\" or both."), call. = FALSE)
  }
  lapply(test, as.character)
}
