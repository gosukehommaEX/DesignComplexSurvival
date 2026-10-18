#' Print an ms_size object
#'
#' @param x An \code{ms_size} object.
#' @param digits Number of significant digits.
#' @param ... Passed to \code{\link{print.ms_power}}.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' res <- ms_size(scn, looks = ms_looks("os"), test = ms_test(),
#'                power = 0.9, solve_for = "events", n = 500)
#' print(res)
#'
#' @export
print.ms_size <- function(x, digits = 4, ...) {
  lab <- switch(x$solve_for,
    events = "final number of events",
    n = "sample size",
    accrual_duration = "accrual duration",
    followup = "follow-up after the end of accrual")
  cat(sprintf("Solved for the %s: %s (target power %s)\n\n", lab,
              format(signif(x$solution, digits + 2)), format(x$target_power)))
  print.ms_power(x, digits = digits, ...)
  invisible(x)
}
