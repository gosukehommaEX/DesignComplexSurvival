#' Transition probabilities after progression
#'
#' Transition probabilities of the post-progression process with states
#' progressed without switch (2), switched (3) and dead (4), between the
#' clock values \code{from} and \code{to}. The clock is the time since
#' randomization for \code{clock = "forward"} and the time since progression
#' for \code{clock = "reset"}; the caller supplies the clock values.
#'
#' @param model An \code{ms_model} object.
#' @param from,to Numeric vectors of equal length with \code{from <= to}.
#' @return A list with \code{r22}, \code{r23}, \code{r33} and the death
#'   hazards \code{c_to} (state 2) and \code{d_to} (state 3) at \code{to}.
#' @noRd
.post_trans <- function(model, from, to) {
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
  list(r22 = r22, r23 = r23, r33 = r33,
       c_to = .pw_value(to, tp, model$death_prog),
       d_to = .pw_value(to, tp, d))
}
