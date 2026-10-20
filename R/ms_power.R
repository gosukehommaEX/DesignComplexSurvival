#' Power of a group sequential log-rank design
#'
#' Computes the analysis times, expected events, bounds, stopping
#' probabilities, power and type I error of a design with a given sample
#' size, without simulation.
#'
#' @details
#' For event-driven looks, the final analysis takes place when the expected
#' number of events of the trigger endpoint reaches \code{events}, or at
#' \code{followup} after the end of accrual. The interim analyses take place
#' when the expected number of events reaches the fractions given in
#' \code{looks}. For calendar-driven looks, \code{events} and
#' \code{followup} are not used.
#'
#' At each look, the standardized log-rank statistics have the mean and
#' covariance of \code{\link{lr_cov}} with the approximations chosen in
#' \code{test}. The bounds are computed by \code{\link{gs_bounds}} and the
#' stopping probabilities by \code{\link{gs_prob}}. The power is the
#' probability of crossing an efficacy bound when the futility bounds are
#' applied. The type I error is computed under \code{\link{ms_null}}
#' applied to the scenario, with event-driven looks at the same numbers of
#' events of the trigger endpoint; the non-binding value ignores the
#' futility bounds and the binding value applies them.
#'
#' The hazard ratios at the bounds are \code{exp(-z / sqrt(p0 p1 d))}, where
#' d is the expected number of events of the tested endpoint.
#'
#' When the looks have different triggers (see \code{\link{ms_looks}}), each
#' trigger endpoint has its own final number of events, and a look driven by
#' it takes place when its expected number of events reaches the fraction of
#' the look times that number. The final numbers are given by \code{events},
#' named by endpoint when the looks are driven by the events of two
#' endpoints, for example \code{c(pfs = 300, os = 250)}. With
#' \code{followup}, or without \code{events} when the final look is
#' calendar-driven, they are the expected numbers of events at the final
#' analysis; the looks driven by an endpoint other than that of the final
#' look then need fractions below 1.
#'
#' With two endpoints (\code{test} from \code{\link{ms_test}} with two
#' endpoints), the statistics of the two endpoints at the looks at which they
#' are tested have the joint distribution of \code{\link{lr_cov}}, with the
#' observed events of the trigger endpoints determining the event-driven
#' looks for \code{variance = "studentized"}. Each endpoint has its own
#' bounds, and the rejection probabilities under the multiple testing
#' procedure are computed by \code{\link{gs_prob_mtp}}. The familywise error
#' rate is computed under \code{\link{ms_null}} applied to the scenario (the
#' global null hypothesis), where the correlations between the endpoints are
#' those of Danzer et al. (2025, Theorem 1), and under each scenario of
#' \code{null_scenarios}, where the hypothesis of one endpoint is true and
#' that of the other is not (a partial null hypothesis). Such scenarios cannot
#' be derived from the multistate model of the design scenario and are given
#' by the user, for example a model with the same OS in both groups and a
#' difference in PFS. Under each null scenario the event-driven looks take
#' place at the same numbers of events of their trigger endpoints.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param n Total sample size.
#' @param looks An \code{ms_looks} object.
#' @param test An \code{ms_test} object.
#' @param events Final number of events of the trigger endpoint (for
#'   event-driven looks), or final numbers of events named by endpoint
#'   (\code{"os"}, \code{"pfs"}).
#' @param followup Time from the end of accrual to the final analysis (for
#'   an event-driven final look, instead of \code{events}).
#' @param null_scenarios For two endpoints, an optional list of
#'   \code{ms_scenario} objects with the allocation and accrual of
#'   \code{scenario}, each named by the endpoint (\code{"os"} or
#'   \code{"pfs"}) whose null hypothesis is true in the scenario, for example
#'   \code{list(os = s_os)}.
#'
#' @return An object of class \code{ms_power}, a list with
#' \describe{
#'   \item{looks}{A data frame with one row per look: \code{time},
#'     \code{enrolled}, the expected events of the tested endpoint
#'     (\code{events}, \code{events_control}, \code{events_treatment}),
#'     \code{info_frac}, the efficacy bound (\code{z_upper}, nominal one-sided
#'     p-value \code{p_upper}, hazard ratio \code{hr_upper}), the futility
#'     bound (\code{z_lower}, \code{hr_lower}), the mean and standard
#'     deviation of the statistic (\code{z_mean}, \code{z_sd}), the
#'     probabilities of stopping for efficacy and futility at the look
#'     (\code{prob_upper}, \code{prob_lower}), the cumulative power
#'     \code{cum_power}, the cumulative alpha spent \code{alpha_spent} and
#'     the cumulative non-binding type I error \code{cum_alpha}. When the
#'     looks are driven by another endpoint, \code{events_trigger} gives the
#'     events of that endpoint.}
#'   \item{power}{Overall power.}
#'   \item{type1}{Type I error, non-binding and binding.}
#'   \item{n, accrual_duration, study_duration, events}{Sample size,
#'     accrual duration, calendar time of the final analysis and expected
#'     events of the tested endpoint at the final analysis.}
#'   \item{settings, scenario, looks_def, test}{The inputs.}
#' }
#' With two endpoints, \code{looks} has one row per endpoint and look at
#' which the endpoint is tested, with the columns above and \code{endpoint},
#' \code{look}, \code{trigger} and \code{events_trigger} (the events of the
#' trigger endpoint at the look, \code{NA} for calendar-driven looks);
#' \code{prob_upper} is the probability of rejecting the hypothesis of the
#' endpoint at the look under the multiple testing procedure,
#' \code{prob_lower} the probability of stopping for futility at the look
#' (first endpoint), \code{cum_power} the cumulative sum of
#' \code{prob_upper}, \code{alpha_spent} the cumulative alpha of the bounds
#' of the endpoint and \code{cum_alpha} the cumulative non-binding
#' probability of rejecting the hypothesis of the endpoint under the global
#' null hypothesis. \code{power} has the powers of the two endpoints (named
#' by endpoint), the probability of rejecting both hypotheses (\code{both})
#' and of rejecting at least one (\code{any}); \code{type1} is the
#' familywise error rate under the global null hypothesis, non-binding and
#' binding; \code{fwer} is a data frame with the familywise error rate under
#' the global null hypothesis and under each scenario of
#' \code{null_scenarios}; \code{events} has the expected events of each
#' endpoint at its last look; \code{final_events} the final numbers of events
#' of the trigger endpoints; and \code{mean} and \code{sigma} are the mean
#' vector and covariance matrix of the statistics in the order of the rows of
#' \code{looks}.
#'
#' @references
#' Danzer MF, Rufibach K, Beyersmann J, Schmidt R (2025). Exhausting the type
#' I error level in event-driven group-sequential designs with a closed
#' testing procedure for progression-free and overall survival.
#' arXiv:2512.08658.
#'
#' @seealso \code{\link{ms_size}}, \code{\link{ms_looks}},
#'   \code{\link{ms_test}}, \code{\link{gs_prob_mtp}}
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' ms_power(scn, n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
#'          test = ms_test(), events = 300)
#'
#' # PFS first, then OS: final PFS analysis at 300 PFS events with an OS
#' # interim analysis, final OS analysis at 220 OS events
#' fl <- ms_scenario(ms_fleischer(4, 12), ms_fleischer(6, 16),
#'                   accrual = ms_accrual(rate = 25))
#' lk <- ms_looks(c("pfs", "os"), fraction = c(1, 1), test = list(c("pfs", "os"), "os"))
#' ms_power(fl, n = 450, looks = lk,
#'          test = ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence"),
#'          events = c(pfs = 300, os = 220))
#'
#' @export
ms_power <- function(scenario, n, looks, test, events = NULL, followup = NULL,
                     null_scenarios = NULL) {
  .check_scenario(scenario)
  .check_positive(n, "n")
  if (!inherits(looks, "ms_looks")) {
    stop("'looks' must be an 'ms_looks' object.", call. = FALSE)
  }
  if (!inherits(test, "ms_test")) {
    stop("'test' must be an 'ms_test' object.", call. = FALSE)
  }
  if (!is.null(events)) {
    .check_events(events)
  }
  if (!is.null(followup) && (!is.numeric(followup) || length(followup) != 1L ||
                             is.na(followup) || followup < 0)) {
    stop("'followup' must be a nonnegative number.", call. = FALSE)
  }
  if (length(test$endpoints) == 1L) {
    if (!is.null(null_scenarios)) {
      stop("'null_scenarios' is used only with two endpoints.", call. = FALSE)
    }
    res <- .design_core(scenario, n, looks, test, events = events,
                        followup = followup)
  } else {
    .check_null_scenarios(null_scenarios, scenario, test$endpoints)
    res <- .design_core2(scenario, n, looks, test, events = events,
                         followup = followup, null_scenarios = null_scenarios)
  }
  structure(res, class = "ms_power")
}
