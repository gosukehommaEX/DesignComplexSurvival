#' Print an ms_looks object
#'
#' @param x An \code{ms_looks} object.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' print(ms_looks("os", fraction = c(0.5, 1)))
#' print(ms_looks(c("pfs", "os"), fraction = c(1, 1), test = list(c("pfs", "os"), "os")))
#'
#' @export
print.ms_looks <- function(x, ...) {
  if (length(x$trigger) != 1L || !is.null(x$test)) {
    lk <- .looks_table(x)
    cat(sprintf("%d look(s):\n", x$k))
    for (j in seq_len(x$k)) {
      when <- if (lk$trigger[j] == "calendar") {
        sprintf("at calendar time %s", format(lk$time[j]))
      } else {
        sprintf("driven by %s events at fraction %s", toupper(lk$trigger[j]),
                format(lk$fraction[j]))
      }
      tested <- if (is.null(x$test)) {
        "all endpoints"
      } else {
        paste(toupper(x$test[[j]]), collapse = " and ")
      }
      cat(sprintf("  %d: %s; tests %s\n", j, when, tested))
    }
  } else if (x$trigger == "calendar") {
    cat(sprintf("%d look(s) at calendar times: %s\n", x$k,
                paste(format(x$time), collapse = ", ")))
  } else {
    cat(sprintf("%d look(s) driven by %s events at fractions: %s\n", x$k,
                toupper(x$trigger), paste(format(x$fraction), collapse = ", ")))
  }
  invisible(x)
}
