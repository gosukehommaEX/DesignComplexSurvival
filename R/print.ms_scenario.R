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
  if (any(c(x$control$switch_pre, x$treatment$switch_pre) > 0)) {
    cat(sprintf("\nPFS at a switch before progression: %s\n",
                if (identical(x$pfs_switch, "censor")) "censored" else "followed up"))
  }
  if (!is.null(x$switch_start) && x$switch_start > 0) {
    cat(sprintf("\nSwitching possible from calendar time %s (time since the start of accrual)\n",
                format(x$switch_start)))
  }
  cat("\n")
  print(x$accrual)
  invisible(x)
}
