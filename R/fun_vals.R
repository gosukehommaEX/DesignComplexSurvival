#' Functions defining the per-patient functionals of the log-rank statistics
#'
#' A per-patient functional is r = delta rho(X) + B(X) with B(x) the integral
#' of beta from 0 to x, where X is the observed follow-up and delta the event
#' indicator. For a patient in group \code{grp} the three functionals are the
#' score (\code{psi}: rho = g, beta = -g lbar), the observed variance
#' (\code{phi}: rho = q (1 - q), beta = (1 - 2 q) g lbar) and the event
#' indicator (\code{dlt}: rho = 1, beta = 0), where q is the proportion of the
#' treatment group among patients at risk, lbar the pooled hazard among
#' patients at risk, and g = 1 - q in the treatment group and g = -q in the
#' control group.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param s Follow-up times.
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @param grp 0 (control group) or 1 (treatment group).
#' @return A list with matrices \code{rho} and \code{beta}, with columns
#'   \code{psi}, \code{phi} and \code{dlt}.
#' @noRd
.fun_vals <- function(scenario, s, endpoint, grp) {
  f <- .alt_funs(scenario, s, endpoint)
  q1 <- f[, "q1"]
  lbar <- f[, "lbar"]
  g <- if (grp == 1L) 1 - q1 else -q1
  ns <- length(s)
  list(rho = cbind(psi = g, phi = q1 * (1 - q1), dlt = rep(1, ns)),
       beta = cbind(psi = -g * lbar, phi = (1 - 2 * q1) * g * lbar,
                    dlt = rep(0, ns)))
}
