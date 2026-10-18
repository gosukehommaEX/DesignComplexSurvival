#' Print an ms_scenario object
#'
#' @param x An \code{ms_scenario} object.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' print(ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20)))
#'
#' @export
print.ms_scenario <- function(x, ...) {
  cat(sprintf("Two-group scenario, allocation %s:1 (treatment:control)\n\n",
              format(x$allocation)))
  cat("Control group. ")
  print(x$control)
  cat("\nTreatment group. ")
  print(x$treatment)
  cat("\n")
  print(x$accrual)
  invisible(x)
}
