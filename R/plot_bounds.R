#' Plot of the bounds of a design
#'
#' @param x An \code{ms_power} object with one endpoint.
#' @param scale \code{"z"}, \code{"hr"} or \code{"p"}.
#' @param x_axis \code{"events"} or \code{"time"}.
#' @param dots Named list of arguments of \code{plot()}.
#' @return The plotted values (see \code{\link{plot.ms_power}}).
#' @noRd
.plot_bounds <- function(x, scale, x_axis, dots) {
  lk <- x$looks
  k <- nrow(lk)
  fin <- is.finite(lk$z_lower)
  d <- switch(scale,
              z = data.frame(look = lk$look, x = lk[[x_axis]], upper = lk$z_upper,
                             lower = ifelse(fin, lk$z_lower, NA_real_), mean = lk$z_mean),
              hr = data.frame(look = lk$look, x = lk[[x_axis]], upper = lk$hr_upper,
                              lower = lk$hr_lower),
              p = data.frame(look = lk$look, x = lk[[x_axis]], upper = lk$p_upper,
                             lower = ifelse(fin, stats::pnorm(lk$z_lower, lower.tail = FALSE),
                                            NA_real_)))
  fut <- any(!is.na(d$lower[-k]))
  col <- .plot_colors()
  xlab <- if (x_axis == "events") {
    sprintf("Expected %s events", toupper(x$settings$endpoint))
  } else {
    "Calendar time"
  }
  ylab <- switch(scale, z = "Standardized log-rank statistic", hr = "Hazard ratio",
                 p = "Nominal one-sided p-value")
  # room for the legend above the values (below them for the p-values on the log scale)
  r <- range(c(d$upper, if (fut) d$lower, d$mean), na.rm = TRUE)
  y <- if (scale == "p") c(r[1L] / (r[2L] / r[1L])^0.4, r[2L]) else c(r[1L], r[2L] + 0.4 * diff(r))
  args <- .plot_args(list(x = range(d$x), y = y, type = "n", xlab = xlab, ylab = ylab,
                          log = if (scale == "p") "y" else ""),
                     dots)
  do.call(graphics::plot, args)
  graphics::lines(d$x, d$upper, type = "b", pch = 19, col = col[["efficacy"]])
  if (fut) {
    i <- which(!is.na(d$lower))
    graphics::lines(d$x[i], d$lower[i], col = col[["futility"]])
    j <- i[i < k]
    graphics::points(d$x[j], d$lower[j], pch = 17, col = col[["futility"]])
  }
  if (scale == "z") {
    graphics::lines(d$x, d$mean, lty = 2, col = col[["mean"]])
  }
  leg <- c("Efficacy bound", if (fut) "Futility bound",
           if (scale == "z") "Mean under the design scenario")
  graphics::legend(if (scale == "p") "bottomright" else "topright", legend = leg,
                   col = c(col[["efficacy"]], if (fut) col[["futility"]],
                           if (scale == "z") col[["mean"]]),
                   lty = c(1, if (fut) 1, if (scale == "z") 2),
                   pch = c(19, if (fut) 17, if (scale == "z") NA), bty = "n")
  d
}
