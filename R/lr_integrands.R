#' Integrands of the expected events and the log-rank moments
#'
#' Per-patient integrands at follow-up times s, before weighting by the
#' enrolled fraction. With y0 and y1 the probabilities of being at risk in
#' each group multiplied by the allocation probabilities, the score mean has
#' integrand y0 y1 / (y0 + y1) times the hazard difference (treatment minus
#' control) and the information has integrand y0 y1 / (y0 + y1)^2 times
#' (y0 lambda0 + y1 lambda1).
#'
#' @param scenario An \code{ms_scenario} object.
#' @param s Follow-up times.
#' @param endpoint \code{"os"} or \code{"pfs"}.
#' @param schoenfeld Whether to add the event-weighted log hazard ratio.
#' @return A matrix with columns \code{ev0}, \code{ev1} (events),
#'   \code{do0}, \code{do1} (dropouts), \code{u} (score mean), \code{v}
#'   (information) and optionally \code{lw}.
#' @noRd
.lr_integrands <- function(scenario, s, endpoint, schoenfeld = FALSE) {
  p0 <- scenario$prob[[1L]]
  p1 <- scenario$prob[[2L]]
  g0 <- .grp_eval(scenario$control, s, endpoint)
  g1 <- .grp_eval(scenario$treatment, s, endpoint)
  y0 <- p0 * g0$R
  y1 <- p1 * g1$R
  ys <- y0 + y1
  ok <- ys > 0
  u <- numeric(length(s))
  v <- numeric(length(s))
  u[ok] <- (y0[ok] * p1 * g1$e[ok] - y1[ok] * p0 * g0$e[ok]) / ys[ok]
  v[ok] <- y0[ok] * y1[ok] * (p0 * g0$e[ok] + p1 * g1$e[ok]) / ys[ok]^2
  out <- cbind(ev0 = p0 * g0$e, ev1 = p1 * g1$e, do0 = p0 * g0$c,
               do1 = p1 * g1$c, u = u, v = v)
  if (schoenfeld) {
    dd <- p0 * g0$e + p1 * g1$e
    pos <- g0$f > 0 & g1$f > 0 & g0$S > 0 & g1$S > 0
    if (any(!pos & dd > 0)) {
      stop("The hazard ratio is zero or infinite where events occur; ",
           "use mean = \"direct\".", call. = FALSE)
    }
    lhr <- numeric(length(s))
    lhr[pos] <- log((g1$f[pos] * g0$S[pos]) / (g0$f[pos] * g1$S[pos]))
    out <- cbind(out, lw = dd * lhr)
  }
  out
}
