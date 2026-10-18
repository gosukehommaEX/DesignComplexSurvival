#' Auxiliary integral for the switch transition after progression
#'
#' Computes K(x), the integral from 0 to x of r(v) exp(Cd(v) - Cy(v)), where
#' r is the switch hazard, Cd the cumulative hazard of death after a switch
#' and Cy the cumulative hazard of leaving the progressed state (switch or
#' death), all on the post-progression clock. The integral has a closed form
#' because the intensities are piecewise constant.
#'
#' @param model An \code{ms_model} object.
#' @param x Evaluation points on the post-progression clock.
#' @return Numeric vector of K(x).
#' @noRd
.post_kfun <- function(model, x) {
  tp <- model$time_post
  r <- model$switch_rate
  y <- r + model$death_prog
  d <- model$death_switch
  z <- d - y
  nb <- length(tp)
  cd0 <- .pw_cumhaz(tp, tp, d)
  cy0 <- .pw_cumhaz(tp, tp, y)
  scale0 <- r * exp(cd0 - cy0)
  len <- c(diff(tp), 0)
  inc <- scale0 * len * .expm1_ratio(z * len)
  kcum <- c(0, cumsum(inc[-nb]))
  x <- pmax(x, 0)
  idx <- pmax(findInterval(x, tp), 1L)
  h <- x - tp[idx]
  kcum[idx] + scale0[idx] * h * .expm1_ratio(z[idx] * h)
}
