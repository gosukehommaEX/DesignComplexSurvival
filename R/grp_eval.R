#' Quantities of one group used in the log-rank calculations
#'
#' For PFS, the censoring hazard is \code{dropout + dropout_pfs}, and, with
#' \code{model$pfs_switch = "censor"}, a switch before progression censors
#' PFS as well: its density P1(s) \code{switch_pre}(s) G(s) is added to the
#' density of an observed censoring.
#'
#' @param model An \code{ms_model} object.
#' @param s Follow-up times since randomization.
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @return A list with the survival function \code{S}, density \code{f},
#'   probability of being at risk \code{R} (no event and no censoring),
#'   density of an observed event \code{e} and density of an observed
#'   censoring \code{c} (dropout, and for PFS the censoring of PFS only and a
#'   censoring switch).
#' @noRd
.grp_eval <- function(model, s, endpoint) {
  ep <- .ep_eval(model, s, endpoint)
  rate <- model$dropout
  if (endpoint == "pfs" && any(model$dropout_pfs > 0)) {
    rate <- rate + model$dropout_pfs
  }
  g <- exp(-.pw_cumhaz(s, model$time, rate))
  r <- ep$surv * g
  cr <- .pw_value(pmax(s, 0), model$time, rate)
  if (endpoint == "pfs" && identical(model$pfs_switch, "censor") &&
      any(model$switch_pre > 0)) {
    cr <- cr + .pw_value(pmax(s, 0), model$time, model$switch_pre)
  }
  list(S = ep$surv, f = ep$dens, R = r, e = ep$dens * g, c = r * cr)
}
