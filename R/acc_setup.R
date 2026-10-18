#' Accrual pattern for a given sample size
#'
#' @param accrual An \code{ms_accrual} object.
#' @param n Total sample size.
#' @return A list with \code{n}, the accrual duration \code{duration}, the
#'   accrual breaks \code{alpha} (from 0 to the accrual duration) and the
#'   accrual densities \code{rho} in each interval, as fractions of \code{n}
#'   per unit time.
#' @noRd
.acc_setup <- function(accrual, n) {
  if (accrual$type == "rate") {
    tm <- accrual$time
    rate <- accrual$rate
    k <- length(tm)
    cum_start <- c(0, cumsum(rate[-k] * diff(tm)))
    cum_end <- c(cum_start[-1L], Inf)
    i <- which(cum_end >= n)[1L]
    dur <- tm[i] + (n - cum_start[i]) / rate[i]
    alpha <- c(tm[tm < dur], dur)
    rho <- rate[seq_len(length(alpha) - 1L)] / n
  } else {
    dur <- accrual$duration
    alpha <- c(accrual$time, dur)
    rho <- accrual$prop / diff(alpha)
  }
  list(n = n, duration = dur, alpha = alpha, rho = rho)
}
