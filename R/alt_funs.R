#' Functions used for the score covariance under the alternative
#'
#' @param scenario An \code{ms_scenario} object.
#' @param s Follow-up times.
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @return A matrix with columns \code{q1} (proportion of the treatment
#'   group among patients at risk), \code{lbar} (pooled hazard among
#'   patients at risk), the densities of an observed event \code{e0},
#'   \code{e1} and of an observed dropout \code{c0}, \code{c1}, and the
#'   probabilities of being at risk \code{R0}, \code{R1}.
#' @noRd
.alt_funs <- function(scenario, s, endpoint) {
  p0 <- scenario$prob[[1L]]
  p1 <- scenario$prob[[2L]]
  g0 <- .grp_eval(scenario$control, s, endpoint)
  g1 <- .grp_eval(scenario$treatment, s, endpoint)
  ys <- p0 * g0$R + p1 * g1$R
  ok <- ys > 0
  q1 <- rep(p1, length(s))
  lbar <- numeric(length(s))
  q1[ok] <- p1 * g1$R[ok] / ys[ok]
  lbar[ok] <- (p0 * g0$e[ok] + p1 * g1$e[ok]) / ys[ok]
  cbind(q1 = q1, lbar = lbar, e0 = g0$e, e1 = g1$e, c0 = g0$c, c1 = g1$c,
        R0 = g0$R, R1 = g1$R)
}
