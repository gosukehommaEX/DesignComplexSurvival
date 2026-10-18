#' Covariance matrix of the standardized log-rank statistics
#'
#' @param scenario An \code{ms_scenario} object.
#' @param acc Output of \code{.acc_setup()}.
#' @param t Calendar times (increasing).
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @param core Output of \code{.lr_core()} at \code{t}.
#' @param variance \code{"null"} (unit variances, correlations
#'   sqrt(V_j / V_k)) or \code{"alternative"} (Luo et al. 2019).
#' @return A covariance matrix.
#' @noRd
.z_cov <- function(scenario, acc, t, endpoint, core, variance) {
  info <- core$info
  if (variance == "null") {
    return(sqrt(outer(info, info, pmin) / outer(info, info, pmax)))
  }
  n <- acc$n
  cv <- .lr_alt_cov(scenario, acc, t, endpoint, core$score_mean / n)
  ipp <- info / n
  cv / sqrt(outer(ipp, ipp))
}
