#' Transition probabilities after progression
#'
#' Transition probabilities of the post-progression process with states
#' progressed without switch (2), switched (3) and dead (4), between the
#' clock values \code{from} and \code{to}. The clock is the time since
#' randomization for \code{clock = "forward"} and the time since progression
#' for \code{clock = "reset"}; the caller supplies the clock values. With
#' \code{switch_clock = "switch"} or an accelerated failure time effect of
#' the switch, the quantities are computed by \code{.post_general()}.
#'
#' @param model An \code{ms_model} object.
#' @param from,to Numeric vectors of equal length with \code{from <= to}.
#' @return A list with \code{r22}, \code{r23}, \code{r33}, the death
#'   hazards \code{c_to} (state 2) and \code{d_to} (after a switch at
#'   progression) at \code{to}, and the death density \code{f23} at
#'   \code{to} of the patients who switched after progression.
#' @noRd
.post_trans <- function(model, from, to) {
  if (identical(model$switch_clock, "switch") ||
      (!is.null(model$switch_af) && model$switch_af != 1)) {
    return(.post_general(model, from, to))
  }
  tp <- model$time_post
  y <- model$switch_rate + model$death_prog
  d <- model$death_switch
  cy_from <- .pw_cumhaz(from, tp, y)
  cy_to <- .pw_cumhaz(to, tp, y)
  cd_from <- .pw_cumhaz(from, tp, d)
  cd_to <- .pw_cumhaz(to, tp, d)
  r22 <- exp(-(cy_to - cy_from))
  r33 <- exp(-(cd_to - cd_from))
  nf <- length(from)
  if (all(model$switch_rate == 0) || nf == 0L) {
    r23 <- numeric(nf)
  } else {
    kf <- .post_kfun(model, c(from, to))
    r23 <- exp(cy_from - cd_to) * (kf[nf + seq_len(nf)] - kf[seq_len(nf)])
    if (any(!is.finite(r23))) {
      stop("Numerical overflow in the post-progression transitions; ",
           "use a larger time unit.", call. = FALSE)
    }
    r23 <- pmax(r23, 0)
  }
  d_to <- .pw_value(to, tp, d)
  list(r22 = r22, r23 = r23, r33 = r33,
       c_to = .pw_value(to, tp, model$death_prog),
       d_to = d_to, f23 = r23 * d_to)
}
