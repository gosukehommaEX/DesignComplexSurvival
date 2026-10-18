#' Print an ms_test object
#'
#' @param x An \code{ms_test} object.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' print(ms_test(futility = "beta"))
#'
#' @export
print.ms_test <- function(x, ...) {
  cat(sprintf("One-sided log-rank test of %s at level %s\n", toupper(x$endpoints),
              format(x$alpha)))
  cat(sprintf("  Efficacy: alpha spending \"%s\"\n", x$sf_eff))
  fut <- switch(x$futility,
    none = "none",
    fixed = sprintf("fixed bounds on the %s scale (%s)",
                    if (x$fut_scale == "z") "Z" else "hazard ratio",
                    paste(format(x$fut_bound), collapse = ", ")),
    beta = sprintf("beta spending \"%s\" with beta = %s", x$sf_fut, format(x$beta)))
  cat(sprintf("  Futility: %s%s\n", fut,
              if (x$futility == "none") "" else if (x$binding) " (binding)" else " (non-binding)"))
  cat(sprintf("  Mean: \"%s\"; variance: \"%s\"; spending time: \"%s\"\n",
              x$mean, x$variance, x$spending_time))
  invisible(x)
}
