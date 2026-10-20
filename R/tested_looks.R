#' Looks at which each endpoint is tested
#'
#' @param looks An \code{ms_looks} object.
#' @param endpoints Tested endpoints (\code{test$endpoints}).
#' @return A list named by endpoint with the indices of the looks at which
#'   the endpoint is tested.
#' @noRd
.tested_looks <- function(looks, endpoints) {
  k <- looks$k
  test <- looks$test
  if (is.null(test)) {
    out <- rep(list(seq_len(k)), length(endpoints))
  } else {
    extra <- setdiff(unlist(test), endpoints)
    if (length(extra) > 0L) {
      stop(sprintf("The looks test %s, which is not an endpoint of the test.",
                   toupper(extra[1L])), call. = FALSE)
    }
    out <- lapply(endpoints, function(e) {
      which(vapply(test, function(x) e %in% x, logical(1)))
    })
    empty <- lengths(out) == 0L
    if (any(empty)) {
      stop(sprintf("%s is not tested at any look.", toupper(endpoints[empty][1L])),
           call. = FALSE)
    }
  }
  names(out) <- endpoints
  if (sum(lengths(out)) > 20L) {
    stop("The two endpoints may be tested at most 20 times in total.", call. = FALSE)
  }
  out
}
