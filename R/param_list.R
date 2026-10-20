#' Spending function parameters of two endpoints
#'
#' @param param \code{NULL}, a list with one element per endpoint
#'   (\code{NULL} or numeric), or a numeric vector with one value for both
#'   endpoints or one value per endpoint (not for \code{sf = "user"}).
#' @param sf Spending functions of the two endpoints.
#' @param name Argument name for messages.
#' @return A list of length 2.
#' @noRd
.param_list <- function(param, sf, name) {
  if (is.null(param)) {
    return(list(NULL, NULL))
  }
  if (is.list(param) && length(param) == 2L &&
      all(vapply(param, function(x) is.null(x) || (is.numeric(x) && !anyNA(x)),
                 logical(1)))) {
    return(unname(param))
  }
  if (!is.list(param) && is.numeric(param) && length(param) %in% 1:2 && !anyNA(param) &&
      !any(sf == "user")) {
    return(as.list(rep_len(as.numeric(param), 2L)))
  }
  stop(sprintf(paste0("'%s' must be a list with one element per endpoint, or one number ",
                      "or one number per endpoint."), name), call. = FALSE)
}
