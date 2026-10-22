#' Multivariate normal probability of a rectangle by the Genz-Bretz algorithm
#'
#' Used by \code{.pmvn()} for rectangles of four or more coordinates whose
#' covariance matrix is not that of a Markov sequence, and by
#' \code{gs_prob_mtp()} for the rectangles that involve both endpoints and
#' have four or more coordinates, for which the Miwa algorithm can be
#' inaccurate (errors of 1e-4 to 3e-3 with 128 grid points were found for
#' such rectangles of designs with two endpoints, and errors of about 1e-3 to
#' above 0.1 for one endpoint with four to eight looks and
#' \code{variance = "alternative"} or \code{"studentized"}).
#' The quasi-Monte Carlo integration of \code{mvtnorm::pmvnorm()} uses a
#' fixed seed, so that the same input always gives the same value, and
#' \code{pmvnorm()} restores the state of the random number generator of the
#' session. The absolute error tolerance can be changed with
#' \code{options(DesignComplexSurvival.genz_abseps = ...)} (default 1e-5).
#'
#' @param lower,upper Limits (may be infinite).
#' @param mean Mean vector.
#' @param sigma Covariance matrix.
#' @return The probability.
#' @noRd
.pmvn_genz <- function(lower, upper, mean, sigma) {
  abseps <- getOption("DesignComplexSurvival.genz_abseps", 1e-5)
  alg <- mvtnorm::GenzBretz(maxpts = 1e6, abseps = abseps, releps = 0)
  p <- mvtnorm::pmvnorm(lower = lower, upper = upper, mean = mean, sigma = sigma,
                        algorithm = alg, seed = 1L)
  min(max(as.numeric(p), 0), 1)
}
