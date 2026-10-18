#' Calendar times at which the expected number of events reaches targets
#'
#' @param scenario An \code{ms_scenario} object.
#' @param acc Output of \code{.acc_setup()}.
#' @param targets Target numbers of events (positive).
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @return Numeric vector of calendar times.
#' @noRd
.event_times <- function(scenario, acc, targets, endpoint) {
  ftot <- function(x) {
    e <- .events_core(scenario, acc, x, endpoint)
    e$ev0 + e$ev1
  }
  out <- numeric(length(targets))
  hi <- max(acc$duration, 1)
  dhi <- ftot(hi)
  for (i in seq_along(targets)) {
    d <- targets[i]
    iter <- 0L
    while (dhi < d) {
      hi_new <- 2 * hi
      d_new <- ftot(hi_new)
      if (d_new - dhi <= 1e-10 * d_new || iter >= 60L) {
        stop(sprintf(paste0("The expected number of %s events (at most about %s) ",
                            "cannot reach %s."),
                     toupper(endpoint), format(signif(d_new, 6)),
                     format(signif(d, 6))), call. = FALSE)
      }
      hi <- hi_new
      dhi <- d_new
      iter <- iter + 1L
    }
    out[i] <- stats::uniroot(function(x) ftot(x) - d, lower = 0, upper = hi,
                             f.lower = -d, f.upper = dhi - d,
                             tol = 1e-10)$root
  }
  out
}
