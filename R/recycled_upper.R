#' Efficacy bounds of an endpoint at the level after recycling
#'
#' Computes the efficacy bounds of an endpoint at the level at which it is
#' tested after the hypothesis of the other endpoint has been rejected
#' (graphical procedure of Maurer and Bretz 2013), with the spending function
#' of the endpoint over all its looks, and checks that they do not exceed the
#' bounds at the initial level (nominal levels that do not decrease with the
#' level, Maurer and Bretz 2013, Section 3.1).
#'
#' @param info_frac Information fractions of the looks of the endpoint.
#' @param level Level after recycling.
#' @param sf_eff,param_eff Spending function and its parameter.
#' @param corr Correlation matrix of the statistics under the null
#'   hypothesis.
#' @param initial Efficacy bounds at the initial level.
#' @param lower Binding futility bounds at the looks (the last value is not
#'   used), or \code{NULL} for none.
#' @param endpoint Endpoint, for the error messages.
#' @return The efficacy bounds after recycling (set to the initial bounds
#'   where they exceed them by less than 1e-8).
#' @noRd
.recycled_upper <- function(info_frac, level, sf_eff, param_eff, corr, initial,
                            lower = NULL, endpoint) {
  k <- length(info_frac)
  b <- tryCatch({
    if (is.null(lower)) {
      gs_bounds(info_frac, alpha = level, sf_eff = sf_eff, param_eff = param_eff,
                corr = corr)$upper
    } else {
      gs_bounds(info_frac, alpha = level, sf_eff = sf_eff, param_eff = param_eff,
                futility = "fixed", fut_bound = lower[seq_len(k - 1L)], binding = TRUE,
                corr = corr)$upper
    }
  }, error = function(e) {
    stop(sprintf("Cannot compute the efficacy bounds of %s at the level %s after recycling: %s",
                 toupper(endpoint), format(level), conditionMessage(e)), call. = FALSE)
  })
  bad <- which(b > initial + 1e-8)
  if (length(bad) > 0L) {
    stop(sprintf(paste("The efficacy bound of %s at its look %d at the level %s after recycling",
                       "exceeds the bound at the initial level: the nominal levels of the",
                       "spending function must not decrease with the level (Maurer and Bretz",
                       "2013, Section 3.1)."), toupper(endpoint), bad[1L], format(level)),
         call. = FALSE)
  }
  pmin(b, initial)
}
