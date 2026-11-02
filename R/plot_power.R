#' Plot of the cumulative stopping probabilities of a design
#'
#' @param x An \code{ms_power} object with one endpoint.
#' @param x_axis \code{"events"} or \code{"time"}.
#' @param dots Named list of arguments of \code{plot()}.
#' @return The plotted values (see \code{\link{plot.ms_power}}).
#' @noRd
.plot_power <- function(x, x_axis, dots) {
  lk <- x$looks
  l0 <- x[["looks_null"]]
  if (is.null(l0) || is.null(lk$cum_alpha)) {
    stop("'x' has no looks under the null scenario; compute it again with ms_power().",
         call. = FALSE)
  }
  k <- nrow(lk)
  d <- data.frame(look = lk$look, x = lk[[x_axis]], cum_efficacy = lk$cum_power,
                  cum_futility = cumsum(c(lk$prob_lower[-k], 0)), x_null = l0[[x_axis]],
                  cum_alpha = lk$cum_alpha)
  fut <- x$settings$futility != "none"
  col <- .plot_colors()
  xlab <- if (x_axis == "events") {
    sprintf("Expected %s events", toupper(x$settings$endpoint))
  } else {
    "Calendar time"
  }
  # room for the legend above the values
  top <- min(1, 1.5 * max(d$cum_efficacy, d$cum_futility, d$cum_alpha))
  args <- .plot_args(list(x = range(d$x, d$x_null), y = c(0, top), type = "n", xlab = xlab,
                          ylab = "Cumulative probability"), dots)
  do.call(graphics::plot, args)
  graphics::lines(d$x, d$cum_efficacy, type = "b", pch = 19, col = col[["efficacy"]])
  if (fut) {
    graphics::lines(d$x, d$cum_futility, type = "b", pch = 17, col = col[["futility"]])
  }
  graphics::lines(d$x_null, d$cum_alpha, type = "b", pch = 15, lty = 2, col = col[["alpha"]])
  leg <- c("Stopping for efficacy (design scenario)",
           if (fut) "Stopping for futility (design scenario)",
           "Type I error (null scenario, non-binding)")
  graphics::legend("topleft", legend = leg,
                   col = c(col[["efficacy"]], if (fut) col[["futility"]], col[["alpha"]]),
                   lty = c(1, if (fut) 1, 2), pch = c(19, if (fut) 17, 15), bty = "n")
  d
}
