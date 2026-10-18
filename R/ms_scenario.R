#' Two-group trial scenario
#'
#' Combines the models of the control and treatment groups, the allocation
#' ratio and the accrual.
#'
#' @param control,treatment \code{ms_model} objects of the control and
#'   treatment groups.
#' @param allocation Allocation ratio, treatment to control (for example 2
#'   for 2:1).
#' @param accrual An \code{ms_accrual} object.
#'
#' @return An object of class \code{ms_scenario}, a list with the two
#'   models, \code{allocation}, the allocation probabilities \code{prob}
#'   (control, treatment) and \code{accrual}.
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
ms_scenario <- function(control, treatment, allocation = 1, accrual) {
  if (!inherits(control, "ms_model") || !inherits(treatment, "ms_model")) {
    stop("'control' and 'treatment' must be 'ms_model' objects.", call. = FALSE)
  }
  if (missing(accrual) || !inherits(accrual, "ms_accrual")) {
    stop("'accrual' must be an 'ms_accrual' object.", call. = FALSE)
  }
  .check_positive(allocation, "allocation")
  p1 <- allocation / (1 + allocation)
  structure(list(control = control, treatment = treatment,
                 allocation = allocation,
                 prob = c(control = 1 - p1, treatment = p1),
                 accrual = accrual),
            class = "ms_scenario")
}
