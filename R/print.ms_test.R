#' Print an ms_test object
#'
#' @param x An \code{ms_test} object.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' print(ms_test(futility = "beta"))
#' print(ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence"))
#' print(ms_test(endpoints = c("pfs", "os"), mtp = "bonferroni",
#'               weights = c(0.8, 0.2), transition = c(1, 1)))
#'
#' @export
print.ms_test <- function(x, ...) {
  if (length(x$endpoints) == 2L) {
    ep <- toupper(x$endpoints)
    cat(sprintf("One-sided log-rank tests of %s and %s at familywise level %s\n", ep[1L],
                ep[2L], format(x$alpha)))
    if (x$mtp == "bonferroni") {
      cat(sprintf("  Bonferroni: %s at level %s, %s at level %s\n", ep[1L],
                  format(x$weights[1L] * x$alpha), ep[2L], format(x$weights[2L] * x$alpha)))
      g <- x$transition
      if (!is.null(g) && any(g > 0)) {
        lev <- x$alpha * (x$weights + rev(g * x$weights))
        cat(sprintf(paste0("  Recycling: fraction %s from %s to %s and %s from %s to %s; ",
                           "%s at level %s after %s is rejected, %s at level %s after %s ",
                           "is rejected\n"),
                    format(g[1L]), ep[1L], ep[2L], format(g[2L]), ep[2L], ep[1L],
                    ep[2L], format(lev[2L]), ep[1L], ep[1L], format(lev[1L]), ep[2L]))
      }
    } else {
      cat(sprintf("  Fixed sequence: %s, then %s (%s)\n", ep[1L], ep[2L], x$hierarchical))
    }
    cat(sprintf("  Efficacy: alpha spending \"%s\" (%s), \"%s\" (%s)\n", x$sf_eff[1L],
                ep[1L], x$sf_eff[2L], ep[2L]))
  } else {
    cat(sprintf("One-sided log-rank test of %s at level %s\n", toupper(x$endpoints),
                format(x$alpha)))
    cat(sprintf("  Efficacy: alpha spending \"%s\"\n", x$sf_eff))
  }
  fut <- switch(x$futility,
    none = "none",
    fixed = sprintf("fixed bounds on the %s scale (%s)",
                    if (x$fut_scale == "z") "Z" else "hazard ratio",
                    paste(format(x$fut_bound), collapse = ", ")),
    beta = sprintf("beta spending \"%s\" with beta = %s", x$sf_fut, format(x$beta)))
  cat(sprintf("  Futility%s: %s%s\n",
              if (length(x$endpoints) == 2L) sprintf(" (%s)", toupper(x$endpoints[1L])) else "",
              fut,
              if (x$futility == "none") "" else if (x$binding) " (binding)" else " (non-binding)"))
  cat(sprintf("  Mean: \"%s\"; variance: \"%s\"; spending time: \"%s\"\n",
              x$mean, x$variance, x$spending_time))
  invisible(x)
}
