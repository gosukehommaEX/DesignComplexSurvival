#' Break points of the integrands of a scenario
#'
#' Returns the break points of the intensities of both groups. For a model
#' with progression and \code{clock = "reset"}, the OS density is a
#' convolution over the time of progression and its derivative jumps at the
#' sums of the pre-progression and post-progression break points, which are
#' added as well. With \code{switch_clock = "switch"} the post-progression
#' survival also has breaks at the sums a + b of two post-progression break
#' points, and with an accelerated failure time effect of the switch (factor
#' phi) at (a - (1 - phi) b) / phi. With \code{clock = "reset"} these are
#' added after each pre-progression break point; with
#' \code{clock = "forward"} they are added as they are, together with the
#' sums of the pre-progression and post-progression break points
#' (\code{"switch"}) or the values (a - (1 - phi) b) / phi for
#' pre-progression break points b (accelerated failure time effect).
#'
#' @param scenario An \code{ms_scenario} object.
#' @return Sorted unique break points.
#' @noRd
.scn_breaks <- function(scenario) {
  m <- list(scenario$control, scenario$treatment)
  sort(unique(unlist(lapply(m, function(x) {
    b <- c(x$time, x$time_post)
    has_post <- any(x$prog > 0) || (any(x$switch_pre > 0) && any(x$prog_switch > 0))
    if (x$clock == "reset" && has_post) {
      b <- c(b, as.vector(outer(x$time, x$time_post, "+")))
    }
    tp <- x$time_post
    phi <- x$switch_af
    if (identical(x$switch_clock, "switch") && any(x$prog > 0)) {
      e <- as.vector(outer(tp, tp, "+"))
      b <- c(b, if (x$clock == "reset") as.vector(outer(x$time, e, "+"))
                else c(e, as.vector(outer(x$time, tp, "+"))))
    } else if (!is.null(phi) && phi != 1 && any(x$prog > 0)) {
      aft <- function(a, v) (a - (1 - phi) * v) / phi
      e <- as.vector(outer(tp, tp, aft))
      b <- c(b, if (x$clock == "reset") as.vector(outer(x$time, e, "+"))
                else c(e, as.vector(outer(tp, x$time, aft))))
    }
    b[b >= 0]
  }))))
}
