#' Power targeted by ms_size()
#'
#' @param res Output of \code{.design_core()} or \code{.design_core2()}.
#' @param target \code{"first"}, \code{"each"}, \code{"both"} or \code{"any"}.
#' @return The power of the first endpoint (\code{"first"}; the power of the
#'   design for one endpoint), the smaller of the powers of the two endpoints
#'   (\code{"each"}), the probability of rejecting both hypotheses
#'   (\code{"both"}) or of rejecting at least one (\code{"any"}).
#' @noRd
.target_power <- function(res, target) {
  pw <- res$power
  if (length(pw) == 1L) {
    return(pw)
  }
  switch(target,
    first = pw[[1L]],
    each = min(pw[[1L]], pw[[2L]]),
    both = pw[["both"]],
    any = pw[["any"]])
}
