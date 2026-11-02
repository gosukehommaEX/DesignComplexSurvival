#' Plot an ms_power object
#'
#' Draws the bounds of a group sequential design, the cumulative probabilities
#' of stopping at the looks, or the survival functions of the groups of the
#' design scenario.
#'
#' @details
#' \describe{
#'   \item{\code{"bounds"}}{The efficacy and futility bounds at the looks on the
#'     scale \code{scale}: the standardized log-rank statistic (\code{"z"},
#'     with the mean of the statistic under the design scenario as a dashed
#'     line), the hazard ratio (\code{"hr"}, \code{hr_upper} and
#'     \code{hr_lower} of \code{x$looks}) or the nominal one-sided p-value
#'     (\code{"p"}, on a log scale, with \code{pnorm(z_lower, lower.tail =
#'     FALSE)} for the futility bounds). The futility bound of the last look is
#'     the efficacy bound, to which the line of the futility bounds is drawn.}
#'   \item{\code{"power"}}{The cumulative probabilities of stopping for efficacy
#'     (\code{cum_power} of \code{x$looks}) and for futility at the interim
#'     looks under the design scenario, and the cumulative non-binding type I
#'     error under the null scenario (\code{cum_alpha}), drawn at the looks
#'     under the null scenario (\code{x$looks_null}).}
#'   \item{\code{"survival"}}{The survival functions of OS and PFS of the
#'     control and treatment groups (\code{\link{ms_surv}}) from 0 to
#'     \code{max_time}. PFS is drawn when it differs from OS in a group, that
#'     is when the model allows progression. With switching from a calendar
#'     time (\code{switch_start} of \code{\link{ms_scenario}}), OS is that of
#'     the patients who enter after the start of switching, as given by
#'     \code{\link{ms_surv}}, which a note above the plot states.}
#' }
#' The colors are those of \code{grDevices::palette.colors(palette =
#' "Okabe-Ito")}, and the default range of the vertical axis leaves room for
#' the legend above the values of \code{"bounds"} (below them for
#' \code{"p"}) and \code{"power"}. For
#' designs with two endpoints only \code{"survival"} is available.
#'
#' @param x An \code{ms_power} object.
#' @param type \code{"bounds"}, \code{"power"} or \code{"survival"}.
#' @param scale Scale of the bounds for \code{"bounds"}: \code{"z"},
#'   \code{"hr"} or \code{"p"}.
#' @param x_axis Horizontal axis for \code{"bounds"} and \code{"power"}: the
#'   expected events of the tested endpoint (\code{"events"}) or the calendar
#'   time (\code{"time"}).
#' @param max_time Largest time since randomization for \code{"survival"}
#'   (defaults to the calendar time of the final analysis,
#'   \code{x$study_duration}).
#' @param ... Named arguments passed to \code{plot()}, such as \code{main},
#'   \code{xlim}, \code{ylim}, \code{xlab} and \code{ylab}, which replace the
#'   defaults.
#'
#' @return The plotted values, invisibly: for \code{"bounds"} a data frame
#'   with \code{look}, \code{x} (events or calendar time), the efficacy bound
#'   \code{upper}, the futility bound \code{lower} (\code{NA} where there is
#'   none) and, for \code{"z"}, the mean \code{mean}; for \code{"power"} a
#'   data frame with \code{look}, \code{x}, \code{cum_efficacy},
#'   \code{cum_futility}, \code{x_null} (events or calendar time of the looks
#'   under the null scenario) and \code{cum_alpha}; for \code{"survival"} a
#'   data frame with \code{time}, \code{group}, \code{endpoint}, \code{surv}
#'   and the legend label \code{label}, with the note on OS as the attribute
#'   \code{"note"} when switching starts at a calendar time.
#'
#' @seealso \code{\link{ms_power}}, \code{\link{ms_surv}}
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' des <- ms_power(scn, n = 400, looks = ms_looks("os", fraction = c(0.5, 0.75, 1)),
#'                 test = ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2),
#'                 events = 300)
#' plot(des)
#' plot(des, type = "bounds", scale = "hr", x_axis = "time")
#' plot(des, type = "power")
#' plot(des, type = "survival", max_time = 48)
#'
#' @export
plot.ms_power <- function(x, type = c("bounds", "power", "survival"), scale = c("z", "hr", "p"),
                          x_axis = c("events", "time"), max_time = NULL, ...) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  x_axis <- match.arg(x_axis)
  dots <- list(...)
  if (length(dots) > 0L && (is.null(names(dots)) || any(names(dots) == ""))) {
    stop("The arguments in '...' must be named.", call. = FALSE)
  }
  if (type != "survival") {
    if (length(x$settings$endpoint) != 1L) {
      stop(sprintf("type = \"%s\" supports designs with one endpoint only.", type),
           call. = FALSE)
    }
    if (!is.null(max_time)) {
      stop("'max_time' is used only with type = \"survival\".", call. = FALSE)
    }
  }
  res <- switch(type,
                bounds = .plot_bounds(x, scale, x_axis, dots),
                power = .plot_power(x, x_axis, dots),
                survival = .plot_survival(x, max_time, dots))
  invisible(res)
}
