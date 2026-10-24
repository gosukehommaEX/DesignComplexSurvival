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
#' @param control,treatment \code{ms_model} objects of the control and
#'   treatment groups.
#' @param allocation Allocation ratio, treatment to control (for example 2
#'   for 2:1).
#' @param accrual An \code{ms_accrual} object.
#' @param pfs_switch Definition of PFS at a switch before progression:
#'   \code{"event"} (followed up to progression or death) or \code{"censor"}
#'   (censored at the switch). It applies to both groups.
#'
#' @return An object of class \code{ms_scenario}, a list with the two
#'   models, \code{allocation}, the allocation probabilities \code{prob}
#'   (control, treatment), \code{accrual} and \code{pfs_switch}. The models
#'   carry \code{pfs_switch} as well.
#'
#' @seealso \code{\link{ms_model}}, \code{\link{ms_accrual}},
#'   \code{\link{ms_null}}
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#'
#' @export
ms_scenario <- function(control, treatment, allocation = 1, accrual,
                        pfs_switch = c("event", "censor")) {
  pfs_switch <- match.arg(pfs_switch)
  if (!inherits(control, "ms_model") || !inherits(treatment, "ms_model")) {
    stop("'control' and 'treatment' must be 'ms_model' objects.", call. = FALSE)
  }
  if (missing(accrual) || !inherits(accrual, "ms_accrual")) {
    stop("'accrual' must be an 'ms_accrual' object.", call. = FALSE)
  }
  .check_positive(allocation, "allocation")
  p1 <- allocation / (1 + allocation)
  control$pfs_switch <- pfs_switch
  treatment$pfs_switch <- pfs_switch
  structure(list(control = control, treatment = treatment,
                 allocation = allocation,
                 prob = c(control = 1 - p1, treatment = p1),
                 accrual = accrual, pfs_switch = pfs_switch),
            class = "ms_scenario")
}
