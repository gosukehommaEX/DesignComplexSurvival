#' Two-group trial scenario
#'
#' Combines the models of the control and treatment groups, the allocation
#' ratio and the accrual.
#'
#' @details
#' \code{pfs_switch} defines PFS for patients who switch before progression
#' (\code{switch_pre} in \code{\link{ms_model}}). With \code{"event"}, the
#' switch does not change the definition: PFS is the time to progression or
#' death. With \code{"censor"}, PFS is censored at the switch (for example
#' when the start of a new treatment censors PFS); the log-rank test of PFS
#' then compares the hazard \code{prog + death} of the patients who are
#' progression-free and have not switched.
#'
#' \code{switch_start} is a calendar time (time since the start of accrual)
#' before which nobody switches, for example when switching is allowed only
#' after an interim analysis. A patient enrolled at time e can then switch
#' from the follow-up time \code{max(switch_start - e, 0)} on: progressions
#' before that time do not lead to a switch at progression, and patients who
#' progressed before it can switch with \code{switch_rate} from that time
#' on. The distribution of OS then depends on the entry time, and the
#' computations for OS integrate over the entry time. They cover the OS
#' log-rank statistic alone with \code{mean = "direct"}; two endpoints
#' together, an OS test at looks driven by PFS events, and a PFS test at
#' looks driven by OS events with a variance other than \code{"null"} give
#' an error. PFS is not affected by the switching. Under the null hypothesis
#' the score no longer has independent increments when the switching changes
#' the death hazard: the correlations between the looks are computed from
#' the per-patient influence functions instead of the information (see
#' \code{\link{lr_cov}}). \code{switch_start} cannot be combined with
#' switching before progression (\code{switch_pre}).
#'
#' @param control,treatment \code{ms_model} objects of the control and
#'   treatment groups.
#' @param allocation Allocation ratio, treatment to control (for example 2
#'   for 2:1).
#' @param accrual An \code{ms_accrual} object.
#' @param pfs_switch Definition of PFS at a switch before progression:
#'   \code{"event"} (followed up to progression or death) or \code{"censor"}
#'   (censored at the switch). It applies to both groups.
#' @param switch_start Calendar time (time since the start of accrual) from
#'   which switching after progression is possible, nonnegative. The default
#'   0 means that switching is possible from randomization, as in the
#'   models.
#'
#' @return An object of class \code{ms_scenario}, a list with the two
#'   models, \code{allocation}, the allocation probabilities \code{prob}
#'   (control, treatment), \code{accrual}, \code{pfs_switch} and
#'   \code{switch_start}. The models carry \code{pfs_switch} as well.
#'
#' @seealso \code{\link{ms_model}}, \code{\link{ms_accrual}},
#'   \code{\link{ms_null}}
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#'
#' # Switching after progression in the control group, possible only from
#' # month 15 (for example after an interim analysis)
#' ctl_sw <- ms_model(prog = 0.1, death = 0.01, death_prog = 0.05,
#'                    switch_prob = 0.3, switch_rate = 0.04, death_switch = 0.03)
#' trt_sw <- ms_model(prog = 0.06, death = 0.008, death_prog = 0.045)
#' ms_scenario(ctl_sw, trt_sw, accrual = ms_accrual(rate = 20), switch_start = 15)
#'
#' @export
ms_scenario <- function(control, treatment, allocation = 1, accrual,
                        pfs_switch = c("event", "censor"), switch_start = 0) {
  pfs_switch <- match.arg(pfs_switch)
  if (!inherits(control, "ms_model") || !inherits(treatment, "ms_model")) {
    stop("'control' and 'treatment' must be 'ms_model' objects.", call. = FALSE)
  }
  if (missing(accrual) || !inherits(accrual, "ms_accrual")) {
    stop("'accrual' must be an 'ms_accrual' object.", call. = FALSE)
  }
  .check_positive(allocation, "allocation")
  if (!is.numeric(switch_start) || length(switch_start) != 1L || is.na(switch_start) ||
      !is.finite(switch_start) || switch_start < 0) {
    stop("'switch_start' must be a nonnegative number.", call. = FALSE)
  }
  if (switch_start > 0 && any(c(control$switch_pre, treatment$switch_pre) > 0)) {
    stop(paste0("'switch_start' cannot be combined with switching before progression ",
                "('switch_pre')."), call. = FALSE)
  }
  p1 <- allocation / (1 + allocation)
  control$pfs_switch <- pfs_switch
  treatment$pfs_switch <- pfs_switch
  structure(list(control = control, treatment = treatment,
                 allocation = allocation,
                 prob = c(control = 1 - p1, treatment = p1),
                 accrual = accrual, pfs_switch = pfs_switch,
                 switch_start = as.numeric(switch_start)),
            class = "ms_scenario")
}
