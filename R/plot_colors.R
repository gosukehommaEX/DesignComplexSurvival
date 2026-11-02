#' Colors of the plots
#'
#' @return A named vector of colors of the Okabe-Ito palette
#'   (\code{grDevices::palette.colors()}) for the efficacy and futility bounds,
#'   the mean of the statistic, the type I error and the control and
#'   treatment groups.
#' @noRd
.plot_colors <- function() {
  oi <- unname(grDevices::palette.colors(palette = "Okabe-Ito"))
  c(efficacy = oi[6L], futility = oi[7L], mean = oi[1L], alpha = oi[4L],
    control = oi[1L], treatment = oi[6L])
}
