#' Arguments of a plot with the defaults replaced by those given
#'
#' @param defaults Named list of the default arguments of \code{plot()}.
#' @param dots Named list of the arguments given in \code{...}.
#' @return The list \code{defaults} with the elements of \code{dots} replacing
#'   or added to those of the same name.
#' @noRd
.plot_args <- function(defaults, dots) {
  defaults[names(dots)] <- dots
  defaults
}
