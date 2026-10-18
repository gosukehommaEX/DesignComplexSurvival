#' Break points of the integrands of a scenario
#'
#' Returns the break points of the intensities of both groups. For a model
#' with progression and \code{clock = "reset"}, the OS density is a
#' convolution over the time of progression and its derivative jumps at the
#' sums of the pre-progression and post-progression break points, which are
#' added as well.
#'
#' @param scenario An \code{ms_scenario} object.
#' @return Sorted unique break points.
#' @noRd
.scn_breaks <- function(scenario) {
  m <- list(scenario$control, scenario$treatment)
  sort(unique(unlist(lapply(m, function(x) {
    b <- c(x$time, x$time_post)
    if (x$clock == "reset" && any(x$prog > 0)) {
      b <- c(b, as.vector(outer(x$time, x$time_post, "+")))
    }
    b
  }))))
}
