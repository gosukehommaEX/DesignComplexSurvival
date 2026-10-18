#' Expected events and dropouts by group at calendar times
#'
#' @param scenario An \code{ms_scenario} object.
#' @param acc Output of \code{.acc_setup()}.
#' @param t Calendar times.
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @return A data frame with \code{time}, \code{ev0}, \code{ev1}, \code{do0}
#'   and \code{do1} (control and treatment).
#' @noRd
.events_core <- function(scenario, acc, t, endpoint) {
  hfun <- function(s) {
    .lr_integrands(scenario, s, endpoint)[, c("ev0", "ev1", "do0", "do1"), drop = FALSE]
  }
  res <- .acc_integral(acc, t, hfun, breaks = .scn_breaks(scenario),
                       hmax = .scn_hmax(scenario, max(c(t, 1))))
  n <- acc$n
  data.frame(time = t, ev0 = n * res[, "ev0"], ev1 = n * res[, "ev1"],
             do0 = n * res[, "do0"], do1 = n * res[, "do1"])
}
