#' Smallest value of a design parameter that gives a target power
#'
#' Searches over a positive parameter x. Without \code{interval}, x is
#' doubled (or halved) from \code{start} until the power brackets the
#' target. The bracket, or \code{interval}, is then divided into a grid and
#' the first grid interval in which the power reaches the target is refined
#' with \code{uniroot()}, so that the smallest solution is returned when the
#' power is not monotone.
#'
#' @param fun Function of x returning the power.
#' @param target Target power.
#' @param start Starting value (positive).
#' @param interval Optional search interval.
#' @param n_grid Number of grid points.
#' @param what Name of the parameter for messages.
#' @return The solution x.
#' @noRd
.solve_power <- function(fun, target, start, interval = NULL, n_grid = 8L,
                         what = "the parameter") {
  if (is.null(interval)) {
    x <- start
    p <- fun(x)
    it <- 0L
    if (p >= target) {
      repeat {
        x_new <- x / 2
        p_new <- fun(x_new)
        it <- it + 1L
        if (p_new < target) {
          interval <- c(x_new, x)
          break
        }
        x <- x_new
        if (it > 30L) {
          stop(sprintf("The target power is reached for any value of %s.", what),
               call. = FALSE)
        }
      }
    } else {
      repeat {
        x_new <- 2 * x
        p_new <- fun(x_new)
        it <- it + 1L
        if (p_new >= target) {
          interval <- c(x, x_new)
          break
        }
        if (it > 30L || (it > 3L && p_new - p < 1e-7)) {
          stop(sprintf("The target power cannot be reached (largest power about %s at %s = %s).",
                       format(signif(p_new, 4)), what, format(signif(x_new, 6))),
               call. = FALSE)
        }
        x <- x_new
        p <- p_new
      }
    }
  }
  if (!is.numeric(interval) || length(interval) != 2L || anyNA(interval) ||
      interval[1L] <= 0 || interval[2L] <= interval[1L]) {
    stop("'interval' must be two increasing positive numbers.", call. = FALSE)
  }
  grid <- seq(interval[1L], interval[2L], length.out = n_grid)
  pw <- vapply(grid, fun, numeric(1))
  i <- which(pw >= target)[1L]
  if (is.na(i)) {
    stop(sprintf("The target power is not reached in the search interval (largest power %s).",
                 format(signif(max(pw), 4))), call. = FALSE)
  }
  if (i == 1L) {
    stop(sprintf("The target power is already reached at the lower end of the search interval (%s = %s).",
                 what, format(signif(grid[1L], 6))), call. = FALSE)
  }
  stats::uniroot(function(x) fun(x) - target, c(grid[i - 1L], grid[i]),
                 f.lower = pw[i - 1L] - target, f.upper = pw[i] - target,
                 tol = 1e-8 * grid[i])$root
}
