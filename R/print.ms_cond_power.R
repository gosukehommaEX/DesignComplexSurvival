#' Print an ms_cond_power object
#'
#' @param x An \code{ms_cond_power} object.
#' @param digits Number of significant digits.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' des <- ms_power(scn, n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
#'                 test = ms_test(), events = 300)
#' print(ms_cond_power(des, look = 1, z = 1.2, under = "trend"))
#'
#' @export
print.ms_cond_power <- function(x, digits = 4, ...) {
  dist <- switch(x$under,
                 design = "design scenario",
                 null = "null hypothesis",
                 trend = sprintf("current trend (hazard ratio %s)", format(signif(x$hr, digits))),
                 hr = sprintf("hazard ratio %s", format(signif(x$hr, digits))),
                 scenario = "given scenario")
  zs <- paste(vapply(x$z, function(v) format(signif(v, digits)), ""), collapse = ", ")
  cat("Conditional power of a group sequential log-rank design\n")
  if (length(x$z) == 1L) {
    cat(sprintf("  Observed at look %d: z = %s\n", x$look, zs))
  } else {
    cat(sprintf("  Observed at looks 1 to %d: z = %s\n", x$look, zs))
  }
  cat(sprintf("  Distribution: %s\n", dist))
  cat(sprintf("  Conditional power: %s\n", format(signif(x$cond_power, digits))))
  cat(sprintf("  Without the futility bounds of the later interim looks: %s\n",
              format(signif(x$cond_power_nonbinding, digits))))
  print(format(x$looks, digits = digits), row.names = FALSE)
  invisible(x)
}
