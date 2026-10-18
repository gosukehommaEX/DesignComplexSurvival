#' Print an ms_power object
#'
#' @param x An \code{ms_power} object.
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
#' print(res)
#'
#' @export
print.ms_power <- function(x, digits = 4, ...) {
  s <- x$settings
  trig <- if (s$trigger == "calendar") "calendar time" else paste(toupper(s$trigger), "events")
  cat("Group sequential log-rank design\n")
  cat(sprintf("  Tested endpoint: %s; looks driven by %s\n", toupper(s$endpoint), trig))
  cat(sprintf("  Sample size: %s; accrual duration: %s; final analysis at: %s\n",
              format(signif(x$n, digits)), format(signif(x$accrual_duration, digits)),
              format(signif(x$study_duration, digits))))
  cat(sprintf("  Power: %s\n", format(signif(x$power, digits))))
  if (!anyNA(x$type1)) {
    cat(sprintf("  Type I error: %s (non-binding), %s (binding)\n",
                format(signif(x$type1[["nonbinding"]], digits)),
                format(signif(x$type1[["binding"]], digits))))
  }
  cols <- c("look", "time", "events", "info_frac", "z_upper", "p_upper",
            "hr_upper", "z_lower", "hr_lower", "prob_upper", "prob_lower",
            "cum_power")
  print(format(x$looks[, cols], digits = digits), row.names = FALSE)
  invisible(x)
}
