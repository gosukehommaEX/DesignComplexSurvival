#' Print a design with two endpoints
#'
#' The part of \code{print.ms_power()} for two endpoints.
#'
#' @param x An \code{ms_power} object with two endpoints.
#' @param digits Number of significant digits.
#' @return \code{x}, invisibly.
#' @noRd
.print_power2 <- function(x, digits) {
  s <- x$settings
  ep <- toupper(s$endpoint)
  f <- function(v) format(signif(v, digits))
  proc <- if (s$mtp == "bonferroni") {
    sprintf("Bonferroni with weights %s and %s", format(s$weights[1L]), format(s$weights[2L]))
  } else {
    sprintf("fixed sequence, %s", s$hierarchical)
  }
  cat("Group sequential log-rank design with two endpoints\n")
  cat(sprintf("  Endpoints: %s, then %s (%s); looks driven by %s\n", ep[1L], ep[2L], proc,
              .trigger_label(s$trigger)))
  cat(sprintf("  Sample size: %s; accrual duration: %s; final analysis at: %s\n",
              f(x$n), f(x$accrual_duration), f(x$study_duration)))
  pw <- x$power
  cat(sprintf("  Power: %s %s, %s %s; both %s; at least one %s\n", ep[1L], f(pw[[1L]]),
              ep[2L], f(pw[[2L]]), f(pw[["both"]]), f(pw[["any"]])))
  if (!is.null(x$fwer)) {
    for (i in seq_len(nrow(x$fwer))) {
      r <- x$fwer[i, ]
      lab <- if (r$scenario == "global") {
        "global null"
      } else {
        sprintf("%s (%s null)", r$scenario, r$null_endpoints)
      }
      cat(sprintf("  FWER under %s: %s (non-binding), %s (binding)\n", lab,
                  f(r$nonbinding), f(r$binding)))
    }
  }
  lk <- x$looks
  lk$endpoint <- toupper(lk$endpoint)
  cols <- c("endpoint", "look", "time", "events", "info_frac", "z_upper", "p_upper",
            "hr_upper", "z_lower", "hr_lower", "prob_upper", "prob_lower", "cum_power")
  print(format(lk[, cols], digits = digits), row.names = FALSE)
  invisible(x)
}
