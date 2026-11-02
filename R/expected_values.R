#' Expected values of a design when the trial stops at the bounds
#'
#' The trial stops at an interim look when the statistic crosses the efficacy
#' bound or falls below the futility bound, and otherwise at the last look.
#'
#' @param prob_upper,prob_lower Probabilities of continuing to each look and
#'   then crossing the efficacy bound or falling below the futility bound (see
#'   \code{gs_prob()}).
#' @param time,events,enrolled Calendar time, expected events of the tested
#'   endpoint and enrolled patients at each look.
#' @return A data frame with one row: the probabilities of stopping at an
#'   interim look for efficacy (\code{prob_efficacy_early}) and for futility
#'   (\code{prob_futility_early}), and the expected \code{study_duration},
#'   \code{events} and \code{enrolled}, each the sum over the looks of the
#'   probability of stopping at the look times the value at the look. The
#'   probability of stopping at the last look is one minus the probability of
#'   stopping at an interim look.
#' @noRd
.expected_values <- function(prob_upper, prob_lower, time, events, enrolled) {
  early <- seq_len(length(time) - 1L)
  p_stop <- prob_upper[early] + prob_lower[early]
  p_stop <- c(p_stop, max(1 - sum(p_stop), 0))
  data.frame(prob_efficacy_early = sum(prob_upper[early]),
             prob_futility_early = sum(prob_lower[early]),
             study_duration = sum(p_stop * time), events = sum(p_stop * events),
             enrolled = sum(p_stop * enrolled))
}
