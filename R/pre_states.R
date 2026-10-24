#' Probabilities of the progression-free states
#'
#' Probabilities of being progression-free and alive without a switch
#' (state 1) and after a switch before progression (state 5), at the times
#' \code{s} since randomization, without dropout. With H1 the cumulative
#' hazard of leaving state 1 (progression, death and switch) and A5 that of
#' leaving state 5, P1(s) = exp(-H1(s)) and P5(s) is the integral from 0 to s
#' of P1(w) eta(w) exp(-(A5(s) - A5(w))), where eta is the switch hazard
#' \code{switch_pre}. On each interval of \code{time} the integrand is the
#' exponential of a linear function, so the integral has a closed form. The
#' exponential is evaluated at the end of the interval where it is larger,
#' so that its exponent is nonpositive and the remaining factor lies between
#' 0 and 1.
#'
#' @param model An \code{ms_model} object.
#' @param s Times since randomization.
#' @return A list with \code{p1} and \code{p5}.
#' @noRd
.pre_states <- function(model, s) {
  s <- pmax(s, 0)
  tm <- model$time
  h1 <- model$prog + model$death + model$switch_pre
  p1 <- exp(-.pw_cumhaz(s, tm, h1))
  ns <- length(s)
  if (all(model$switch_pre == 0) || ns == 0L) {
    return(list(p1 = p1, p5 = numeric(ns)))
  }
  h5 <- model$prog_switch + model$death_switch_pre
  nt <- length(tm)
  # interval j of time (columns) for each time s (rows)
  lo <- matrix(tm, ns, nt, byrow = TRUE)
  ln <- pmax(pmin(matrix(c(tm[-1L], Inf), ns, nt, byrow = TRUE), s) - lo, 0)
  a_lo <- -matrix(.pw_cumhaz(tm, tm, h1), ns, nt, byrow = TRUE) -
    (.pw_cumhaz(s, tm, h5) - matrix(.pw_cumhaz(tm, tm, h5), ns, nt, byrow = TRUE))
  slope <- matrix(h5 - h1, ns, nt, byrow = TRUE)
  eta <- matrix(model$switch_pre, ns, nt, byrow = TRUE)
  ok <- ln > 0 & eta > 0
  v <- numeric(ns * nt)
  v[ok] <- eta[ok] * exp(a_lo[ok] + pmax(slope[ok], 0) * ln[ok]) * ln[ok] *
    .expm1_ratio(-abs(slope[ok]) * ln[ok])
  list(p1 = p1, p5 = rowSums(matrix(v, ns, nt)))
}
