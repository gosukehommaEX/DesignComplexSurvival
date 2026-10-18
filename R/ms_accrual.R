#' Piecewise uniform accrual
#'
#' Defines the accrual pattern either by accrual rates or by the proportions
#' of patients enrolled in intervals of a given accrual duration.
#'
#' @details
#' With \code{rate}, patients are enrolled at \code{rate[j]} patients per
#' unit time from \code{time[j]} on, and the accrual duration follows from
#' the sample size. The last rate applies until the sample size is reached.
#'
#' With \code{prop} and \code{duration}, a proportion \code{prop[j]} of the
#' patients is enrolled uniformly in the interval starting at
#' \code{time[j]}, and the last interval ends at \code{duration}. The
#' accrual rates then scale with the sample size. When
#' \code{\link{ms_size}} solves for the accrual duration, all interval
#' limits are scaled by the same factor.
#'
#' @param time Interval start times (calendar time from the start of the
#'   study). Must start at 0 and be strictly increasing.
#' @param rate Accrual rates in each interval (patients per unit time).
#' @param prop Proportions of patients enrolled in each interval (summing to
#'   1).
#' @param duration Accrual duration (with \code{prop}).
#'
#' @return An object of class \code{ms_accrual}.
#'
#' @seealso \code{\link{ms_scenario}}
#'
#' @examples
#' # 10 patients per month in the first 3 months, then 30 per month
#' ms_accrual(time = c(0, 3), rate = c(10, 30))
#'
#' # 30 percent of the patients in the first 6 months of an 18-month accrual
#' ms_accrual(time = c(0, 6), prop = c(0.3, 0.7), duration = 18)
#'
#' @export
ms_accrual <- function(time = 0, rate = NULL, prop = NULL, duration = NULL) {
  .check_breaks(time, "time")
  if (!is.null(rate)) {
    if (!is.null(prop) || !is.null(duration)) {
      stop("Specify either 'rate', or 'prop' and 'duration'.", call. = FALSE)
    }
    rate <- .check_rates(rate, length(time), "rate")
    if (rate[length(rate)] <= 0) {
      stop("The last accrual rate must be positive.", call. = FALSE)
    }
    out <- list(type = "rate", time = as.numeric(time), rate = rate)
  } else {
    if (is.null(prop) || is.null(duration)) {
      stop("Specify 'rate', or both 'prop' and 'duration'.", call. = FALSE)
    }
    .check_positive(duration, "duration")
    if (time[length(time)] >= duration) {
      stop("All values of 'time' must be smaller than 'duration'.", call. = FALSE)
    }
    prop <- .check_rates(prop, length(time), "prop", upper = 1)
    if (abs(sum(prop) - 1) > 1e-8) {
      stop("'prop' must sum to 1.", call. = FALSE)
    }
    out <- list(type = "prop", time = as.numeric(time), prop = prop,
                duration = as.numeric(duration))
  }
  structure(out, class = "ms_accrual")
}
