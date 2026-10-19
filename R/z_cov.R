#' Covariance matrix of the standardized log-rank statistics
#'
#' @param scenario An \code{ms_scenario} object.
#' @param acc Output of \code{.acc_setup()}.
#' @param t Calendar times (increasing).
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @param core Output of \code{.lr_core()} at \code{t}.
#' @param variance \code{"null"} (unit variances, correlations
#'   sqrt(V_j / V_k)), \code{"alternative"} (Luo et al. 2019) or
#'   \code{"studentized"} (\code{.zeta_cov()}).
#' @param trigger Endpoint whose events determine each time, or \code{NULL}
#'   (used with \code{"studentized"}).
#' @return A covariance matrix.
#' @noRd
.z_cov <- function(scenario, acc, t, endpoint, core, variance, trigger = NULL) {
  info <- core$info
  if (variance == "null") {
    return(sqrt(outer(info, info, pmin) / outer(info, info, pmax)))
  }
  if (variance == "studentized") {
    return(.zeta_cov(scenario, acc, rep(endpoint, length(t)), t, trigger,
                     "studentized")$sigma)
  }
  n <- acc$n
  cv <- .lr_alt_cov(scenario, acc, t, endpoint, core$score_mean / n)
  ipp <- info / n
  cv / sqrt(outer(ipp, ipp))
}
