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
  cat("Four-state model for one group (clock after progression: \"", x$clock,
      "\")\n", sep = "")
  cat("Before progression (time since randomization):\n")
  print(data.frame(time = x$time, prog = x$prog, death = x$death,
                   switch_prob = x$switch_prob, dropout = x$dropout),
        row.names = FALSE)
  if (any(x$prog > 0)) {
    cat(sprintf("After progression (time since %s):\n",
                if (x$clock == "reset") "progression" else "randomization"))
    print(data.frame(time_post = x$time_post, death_prog = x$death_prog,
                     switch_rate = x$switch_rate, death_switch = x$death_switch),
          row.names = FALSE)
  }
  invisible(x)
}
