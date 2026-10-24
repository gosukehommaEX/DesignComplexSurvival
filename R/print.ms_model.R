#' Print an ms_model object
#'
#' @param x An \code{ms_model} object.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' print(ms_fleischer(median_pfs = 6, median_os = 15))
#'
#' @export
print.ms_model <- function(x, ...) {
  five <- any(x$switch_pre > 0)
  cat(if (five) "Five" else "Four", "-state model for one group (clock after progression: \"",
      x$clock, "\")\n", sep = "")
  cat("Before progression (time since randomization):\n")
  pre <- data.frame(time = x$time, prog = x$prog, death = x$death,
                    switch_prob = x$switch_prob, dropout = x$dropout)
  if (any(x$dropout_pfs > 0)) {
    pre$dropout_pfs <- x$dropout_pfs
  }
  print(pre, row.names = FALSE)
  if (five) {
    cat("Switch before progression (time since randomization):\n")
    print(data.frame(time = x$time, switch_pre = x$switch_pre,
                     prog_switch = x$prog_switch,
                     death_switch_pre = x$death_switch_pre),
          row.names = FALSE)
  }
  if (any(x$prog > 0) || (five && any(x$prog_switch > 0))) {
    cat(sprintf("After progression (time since %s):\n",
                if (x$clock == "reset") "progression" else "randomization"))
    print(data.frame(time_post = x$time_post, death_prog = x$death_prog,
                     switch_rate = x$switch_rate, death_switch = x$death_switch),
          row.names = FALSE)
    if (identical(x$switch_clock, "switch")) {
      cat("Death after a switch: 'death_switch' on the time since the switch.\n")
    }
    if (!is.null(x$switch_af)) {
      cat(sprintf(paste0("Death after a switch: accelerated failure time effect with ",
                         "factor %s ('death_switch' is %s * death_prog).\n"),
                  format(x$switch_af), format(x$switch_af)))
    }
  }
  invisible(x)
}
