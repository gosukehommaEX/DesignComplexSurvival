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
#' With two endpoints, each endpoint has its own efficacy bounds, computed as
#' above from the looks at which it is tested (see \code{\link{ms_looks}}).
#' With \code{mtp = "bonferroni"}, the endpoint k is tested at level
#' \code{weights[k] * alpha}; the levels are not recycled after a rejection.
#' With \code{mtp = "fixed_sequence"}, both endpoints are tested at level
#' \code{alpha} and the second endpoint is tested only after the first has
#' been rejected: with \code{hierarchical = "overall"} at the analysis at
#' which the first is rejected and at its later looks until it is rejected
#' (procedure (b) of Glimm et al. 2010), and with \code{"stagewise"} only at
#' the analysis at which the first is rejected, where the trial stops
#' (procedure (a) of Glimm et al. 2010, Tamhane et al. 2010). The bounds of
#' the second endpoint do not depend on the correlation between the
#' endpoints; with them the familywise error rate is at most \code{alpha}
#' (Glimm et al. 2010, Tamhane et al. 2010). Futility bounds are placed on
#' the first endpoint only: the trial stops when its statistic falls below
#' the futility bound at an interim look at which it has not been rejected.
#' See \code{\link{gs_prob_mtp}} for the stopping rules.
#'
#' @param alpha One-sided type I error (with two endpoints, the familywise
#'   error rate).
#' @param endpoints Tested endpoint, \code{"os"} or \code{"pfs"}, or both in
#'   the order of testing (for example \code{c("pfs", "os")}). With two
#'   endpoints the futility bounds are placed on the first, and with
#'   \code{mtp = "fixed_sequence"} the first is tested first.
#' @param sf_eff,param_eff Alpha spending function and its parameter (see
#'   \code{\link{gs_spend}}). For \code{"user"}, \code{param_eff} gives the
#'   cumulative proportions of \code{alpha} at the looks. With two
#'   endpoints, \code{sf_eff} has one value or one value per endpoint, and
#'   \code{param_eff} is a list with one element per endpoint, or a number or
#'   one number per endpoint.
#' @param futility \code{"none"}, \code{"fixed"} or \code{"beta"} (for the
#'   first endpoint).
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
#' @param mtp Multiple testing procedure: \code{"none"} (one endpoint),
#'   \code{"bonferroni"} or \code{"fixed_sequence"} (two endpoints).
#' @param weights Bonferroni weights of the two endpoints (positive, adding
#'   up to 1; defaults to equal weights).
#' @param hierarchical \code{"overall"} or \code{"stagewise"} (for
#'   \code{mtp = "fixed_sequence"}).
#'
#' @return An object of class \code{ms_test}.
#'
#' @references
#' Glimm E, Maurer W, Bretz F (2010). Hierarchical testing of multiple
#' endpoints in group-sequential trials. Statistics in Medicine 29:219-228.
#'
#' Tamhane AC, Mehta CR, Liu L (2010). Testing a primary and a secondary
#' endpoint in a group sequential design. Biometrics 66:1174-1184.
#'
#' @seealso \code{\link{ms_power}}, \code{\link{gs_bounds}},
#'   \code{\link{gs_prob_mtp}}
#'
#' @examples
#' ms_test()
#' ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2)
#' ms_test(futility = "fixed", fut_bound = 1, fut_scale = "hr")
#' # PFS first, then OS, with Pocock-type spending for OS
#' ms_test(endpoints = c("pfs", "os"), sf_eff = c("OF", "Pocock"),
#'         mtp = "fixed_sequence")
#' ms_test(endpoints = c("pfs", "os"), mtp = "bonferroni", weights = c(0.8, 0.2))
#'
#' @export
ms_test <- function(alpha = 0.025, endpoints = "os", sf_eff = "OF",
                    param_eff = NULL, futility = c("none", "fixed", "beta"),
                    fut_bound = NULL, fut_scale = c("z", "hr"), sf_fut = "OF",
                    param_fut = NULL, beta = 0.1, binding = FALSE,
                    mean = c("direct", "schoenfeld"),
                    variance = c("null", "alternative", "studentized"),
                    spending_time = c("events", "information"),
                    mtp = c("none", "bonferroni", "fixed_sequence"), weights = NULL,
                    hierarchical = c("overall", "stagewise")) {
  futility <- match.arg(futility)
  fut_scale <- match.arg(fut_scale)
  mean <- match.arg(mean)
  variance <- match.arg(variance)
  spending_time <- match.arg(spending_time)
  mtp <- match.arg(mtp)
  hierarchical <- match.arg(hierarchical)
  .check_prob(alpha, "alpha")
  .check_prob(beta, "beta")
  if (!is.character(endpoints) || !(length(endpoints) %in% 1:2) || anyNA(endpoints) ||
      !all(endpoints %in% c("os", "pfs")) || anyDuplicated(endpoints) > 0L) {
    stop("'endpoints' must be \"os\", \"pfs\" or both (in the order of testing).",
         call. = FALSE)
  }
  n_ep <- length(endpoints)
  if (n_ep == 2L && mtp == "none") {
    stop("Two 'endpoints' need mtp = \"bonferroni\" or \"fixed_sequence\".", call. = FALSE)
  }
  if (n_ep == 1L && mtp != "none") {
    stop("'mtp' must be \"none\" for one endpoint.", call. = FALSE)
  }
  if (mtp == "bonferroni") {
    if (is.null(weights)) {
      weights <- c(0.5, 0.5)
    }
    if (!is.numeric(weights) || length(weights) != 2L || anyNA(weights) ||
        any(weights <= 0) || abs(sum(weights) - 1) > 1e-10) {
      stop("'weights' must be two positive numbers adding up to 1.", call. = FALSE)
    }
    weights <- as.numeric(weights)
  } else if (!is.null(weights)) {
    stop("'weights' is used only with mtp = \"bonferroni\".", call. = FALSE)
  }
  sf_types <- c("OF", "Pocock", "HSD", "power", "user")
  if (!is.character(sf_eff) || !(length(sf_eff) %in% c(1L, n_ep)) || anyNA(sf_eff) ||
      !all(sf_eff %in% sf_types)) {
    stop(paste0("'sf_eff' must be one of \"OF\", \"Pocock\", \"HSD\", \"power\", \"user\"",
                if (n_ep == 2L) " (one value or one per endpoint)." else "."),
         call. = FALSE)
  }
  if (n_ep == 2L) {
    sf_eff <- rep_len(sf_eff, 2L)
    param_eff <- .param_list(param_eff, sf_eff, "param_eff")
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
                 spending_time = spending_time, mtp = mtp, weights = weights,
                 hierarchical = hierarchical),
            class = "ms_test")
}
