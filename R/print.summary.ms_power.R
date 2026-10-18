#' Print a summary.ms_power object
#'
#' @param x A \code{summary.ms_power} object.
#' @param digits Number of significant digits.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' res <- ms_power(scn, n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
#'                 test = ms_test(), events = 300)
#' print(summary(res))
#'
#' @export
print.summary.ms_power <- function(x, digits = 4, ...) {
  s <- x$settings
  cat("Group sequential log-rank design: summary\n")
  cat(sprintf("  Tested endpoint: %s; trigger: %s; mean: \"%s\"; variance: \"%s\"\n",
              toupper(s$endpoint), s$trigger, s$mean, s$variance))
  cat(sprintf("  Futility: %s%s; spending time: \"%s\"; alpha: %s\n", s$futility,
              if (s$futility == "none") "" else if (s$binding) " (binding)" else " (non-binding)",
              s$spending_time, format(s$alpha)))
  cat("\nOverview:\n")
  print(signif(x$overview, digits))
  cat("\nBy look:\n")
  print(format(x$looks, digits = digits), row.names = FALSE)
  invisible(x)
}
