#' Test and stopping rules
#'
#' Defines the one-sided log-rank test, the efficacy and futility bounds and
#' the approximation of the distribution of the statistics.
#'
#' @details
#' The efficacy bounds are obtained by alpha spending under the null
#' hypothesis. The correlation of the statistics under the null hypothesis
#' is taken as \code{sqrt(d_j / d_k)}, where d are the expected numbers of
#' events of the tested endpoint at the looks; it is exact under the null
#' hypothesis when the looks are driven by the events of the tested
#' endpoint. The spending time is the event fraction
#' (\code{spending_time = "events"}) or the information fraction under the
#' design scenario (\code{"information"}).
#'
#' Futility bounds may be fixed, on the Z scale (\code{fut_scale = "z"}) or
#' on the hazard ratio scale (\code{fut_scale = "hr"}, converted with
#' \code{z = -log(hr) sqrt(p0 p1 d)}), or obtained by beta spending under
#' the design scenario. With \code{binding = FALSE} the efficacy bounds
#' ignore the futility bounds.
#'
#' @param alpha One-sided type I error.
#' @param endpoints Tested endpoint, \code{"os"} or \code{"pfs"} (one
#'   endpoint in this version).
#' @param sf_eff,param_eff Alpha spending function and its parameter (see
#'   \code{\link{gs_spend}}). For \code{"user"}, \code{param_eff} gives the
#'   cumulative proportions of \code{alpha} at the looks.
#' @param futility \code{"none"}, \code{"fixed"} or \code{"beta"}.
#' @param fut_bound Fixed futility bounds at the interim looks (for
#'   \code{futility = "fixed"}).
#' @param fut_scale Scale of \code{fut_bound}: \code{"z"} or \code{"hr"}.
#' @param sf_fut,param_fut Beta spending function and its parameter.
#' @param beta Total type II error for beta spending in
#'   \code{\link{ms_power}}; \code{\link{ms_size}} uses one minus the target
#'   power instead.
#' @param binding Whether the futility bounds are binding.
#' @param mean Approximation of the mean of the statistic:
#'   \code{"direct"} or \code{"schoenfeld"} (see \code{\link{lr_moments}}).
#' @param variance Approximation of the covariance of the statistics:
#'   \code{"null"}, \code{"alternative"} or \code{"studentized"} (see
#'   \code{\link{lr_cov}}).
#' @param spending_time \code{"events"} or \code{"information"}.
#'
#' @return An object of class \code{ms_test}.
#'
#' @seealso \code{\link{ms_power}}, \code{\link{gs_bounds}}
#'
#' @examples
#' ms_test()
#' ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2)
#' ms_test(futility = "fixed", fut_bound = 1, fut_scale = "hr")
#'
#' @export
ms_test <- function(alpha = 0.025, endpoints = "os", sf_eff = "OF",
                    param_eff = NULL, futility = c("none", "fixed", "beta"),
                    fut_bound = NULL, fut_scale = c("z", "hr"), sf_fut = "OF",
                    param_fut = NULL, beta = 0.1, binding = FALSE,
                    mean = c("direct", "schoenfeld"),
                    variance = c("null", "alternative", "studentized"),
                    spending_time = c("events", "information")) {
  futility <- match.arg(futility)
  fut_scale <- match.arg(fut_scale)
  mean <- match.arg(mean)
  variance <- match.arg(variance)
  spending_time <- match.arg(spending_time)
  .check_prob(alpha, "alpha")
  .check_prob(beta, "beta")
  if (!is.character(endpoints) || length(endpoints) != 1L ||
      !(endpoints %in% c("os", "pfs"))) {
    stop("'endpoints' must be \"os\" or \"pfs\" in this version.", call. = FALSE)
  }
  sf_types <- c("OF", "Pocock", "HSD", "power", "user")
  if (!is.character(sf_eff) || length(sf_eff) != 1L || !(sf_eff %in% sf_types)) {
    stop("'sf_eff' must be one of \"OF\", \"Pocock\", \"HSD\", \"power\", \"user\".",
         call. = FALSE)
  }
  if (!is.character(sf_fut) || length(sf_fut) != 1L || !(sf_fut %in% sf_types)) {
    stop("'sf_fut' must be one of \"OF\", \"Pocock\", \"HSD\", \"power\", \"user\".",
         call. = FALSE)
  }
  if (futility == "fixed") {
    if (!is.numeric(fut_bound) || anyNA(fut_bound) ||
        (fut_scale == "hr" && any(fut_bound <= 0))) {
      stop("'fut_bound' must be numeric (positive on the hazard ratio scale).",
           call. = FALSE)
    }
  }
  if (!is.logical(binding) || length(binding) != 1L || is.na(binding)) {
    stop("'binding' must be TRUE or FALSE.", call. = FALSE)
  }
  structure(list(alpha = alpha, endpoints = endpoints, sf_eff = sf_eff,
                 param_eff = param_eff, futility = futility,
                 fut_bound = fut_bound, fut_scale = fut_scale, sf_fut = sf_fut,
                 param_fut = param_fut, beta = beta, binding = binding,
                 mean = mean, variance = variance,
                 spending_time = spending_time),
            class = "ms_test")
}
