#' State occupation probabilities
#'
#' Probabilities of being in each state of the four-state model at given
#' times since randomization, without dropout.
#'
#' @details
#' The probability of being progression-free and alive is
#' \code{exp(-H(t))}, where H is the cumulative hazard of progression and
#' death. The probabilities of the progressed and switched states are
#' integrals over the time of progression of the progression density times
#' the transition probabilities after progression, which have closed forms
#' for piecewise constant intensities. The integrals are computed with
#' Gauss-Legendre rules.
#'
#' @param model An \code{ms_model} object.
#' @param time Times since randomization (nonnegative).
#'
#' @return A data frame with \code{time} and the probabilities
#'   \code{progression_free}, \code{progressed} (without switch),
#'   \code{switched} and \code{dead}.
#'
#' @seealso \code{\link{ms_model}}, \code{\link{ms_surv}}
#'
#' @examples
#' m <- ms_model(time = c(0, 6), prog = c(0.10, 0.08), death = 0.01,
#'               death_prog = 0.05, switch_prob = 0.3, death_switch = 0.03)
#' ms_prob(m, time = c(6, 12, 24))
#'
#' @export
ms_prob <- function(model, time) {
  if (!inherits(model, "ms_model")) {
    stop("'model' must be an 'ms_model' object.", call. = FALSE)
  }
  .check_times(time, "time")
  p1 <- exp(-.pw_cumhaz(time, model$time, model$prog + model$death))
  if (all(model$prog == 0)) {
    p2 <- numeric(length(time))
    p3 <- numeric(length(time))
  } else {
    st <- .os_conv(model, time, what = "states")
    p2 <- st$p12
    p3 <- st$p13
  }
  data.frame(time = time, progression_free = p1, progressed = p2,
             switched = p3, dead = pmax(1 - p1 - p2 - p3, 0))
}
