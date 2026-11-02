#' Probabilities of the later looks given the statistics at an interim look
#'
#' Conditions the statistics of the looks after \code{look} on the observed
#' statistic at \code{look} (or on the statistics at the looks 1 to
#' \code{look}) and computes, at each later look, the probabilities of
#' continuing to the look and then crossing the efficacy bound or falling
#' below the futility bound.
#'
#' @param mean,sigma Mean vector and covariance matrix of the statistics at
#'   all looks.
#' @param look Interim look of the observation.
#' @param z Observed statistic at \code{look}, or observed statistics at the
#'   looks 1 to \code{look}.
#' @param lower,upper Futility and efficacy bounds at all looks (the last
#'   futility bound equals the last efficacy bound).
#' @return A list with a data frame \code{looks} (\code{look},
#'   \code{cond_mean}, \code{cond_sd}, \code{prob_upper}, \code{prob_lower},
#'   \code{cum_upper}) of the later looks, the conditional power
#'   \code{cond_power} (the sum of \code{prob_upper}) and the conditional power
#'   without the futility bounds of the later interim looks
#'   \code{cond_power_nonbinding}.
#' @noRd
.cond_probs <- function(mean, sigma, look, z, lower, upper) {
  k_max <- length(upper)
  given <- if (length(z) == 1L) look else seq_len(look)
  future <- seq.int(look + 1L, k_max)
  cm <- .cond_moments(mean, sigma, given, future, z)
  lo <- lower[future]
  up <- upper[future]
  pr <- gs_prob(lo, up, cm$mean, cm$sigma)
  lo_nb <- c(rep(-Inf, length(future) - 1L), lo[length(lo)])
  nb <- gs_prob(lo_nb, up, cm$mean, cm$sigma)
  list(looks = data.frame(look = future, cond_mean = cm$mean,
                          cond_sd = sqrt(diag(cm$sigma)), prob_upper = pr$prob_upper,
                          prob_lower = pr$prob_lower, cum_upper = pr$cum_upper),
       cond_power = sum(pr$prob_upper), cond_power_nonbinding = sum(nb$prob_upper))
}
