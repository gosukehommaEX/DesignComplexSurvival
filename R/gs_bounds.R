#' Bounds of a group sequential design
#'
#' Computes one-sided efficacy bounds by alpha spending and, optionally,
#' futility bounds that are fixed or obtained by beta spending, for a
#' vector of standardized statistics with a multivariate normal
#' distribution.
#'
#' @details
#' The efficacy bound at look k solves, under the null hypothesis (mean 0,
#' correlation \code{corr}), the equation that the probability of
#' continuing to look k and then crossing the bound equals the alpha spent
#' at look k. With \code{binding = TRUE} the continuation region at the
#' earlier looks includes the futility bounds; with \code{binding = FALSE}
#' the futility bounds are ignored for the efficacy bounds.
#'
#' With \code{futility = "beta"}, the futility bound at look k solves, under
#' the alternative (mean \code{mean}, covariance \code{sigma_alt}), the
#' equation that the probability of continuing to look k and then falling
#' below the bound equals the beta spent at look k. The mean may take any
#' value at each look, so that non-proportional hazards are allowed. A
#' futility bound above the efficacy bound is set to the efficacy bound. At
#' the last look the futility bound equals the efficacy bound. With
#' \code{binding = TRUE}, a futility bound set to the efficacy bound at an
#' interim look stops every trial at that look, so that no alpha can be
#' spent at the later looks; the function then stops with an error that
#' names that look.
#'
#' The probabilities are computed as in \code{\link{gs_prob}}: with the Miwa
#' algorithm for rectangles of at most three coordinates and for covariance
#' matrices of Markov sequences, such as the default \code{corr}, except for
#' rectangles of nine or more coordinates with two finite limits for some
#' coordinate (nine or more looks with futility bounds at the earlier looks),
#' and with the Genz-Bretz algorithm otherwise (for example for
#' \code{sigma_alt} of \code{variance = "alternative"} or \code{"studentized"}
#' with four or more looks).
#'
#' @param info_frac Information fractions used for spending (strictly
#'   increasing and ending at 1).
#' @param alpha One-sided type I error.
#' @param sf_eff,param_eff Alpha spending function and its parameter (see
#'   \code{\link{gs_spend}}).
#' @param futility \code{"none"}, \code{"fixed"} or \code{"beta"}.
#' @param fut_bound Fixed futility bounds on the Z scale at looks 1 to K - 1
#'   (for \code{futility = "fixed"}).
#' @param sf_fut,param_fut Beta spending function and its parameter.
#' @param beta Total type II error for beta spending.
#' @param binding Whether the futility bounds are binding.
#' @param mean Mean vector of the statistics under the alternative
#'   (required for \code{futility = "beta"}).
#' @param corr Correlation matrix of the statistics under the null
#'   hypothesis. Defaults to \code{sqrt(t_j / t_k)} for j <= k with
#'   \code{t = info_frac}.
#' @param sigma_alt Covariance matrix of the statistics under the
#'   alternative. Defaults to \code{corr}.
#'
#' @return A data frame with one row per look: \code{look},
#'   \code{info_frac}, cumulative \code{alpha_spent} and \code{beta_spent},
#'   the efficacy bound \code{upper}, the futility bound \code{lower}
#'   (\code{-Inf} when there is none) and the one-sided nominal p-value of
#'   the efficacy bound \code{nominal_p}.
#'
#' @references
#' Lan KKG, DeMets DL (1983). Discrete sequential boundaries for clinical
#' trials. Biometrika 70:659-663.
#'
#' @seealso \code{\link{gs_spend}}, \code{\link{gs_prob}},
#'   \code{\link{ms_power}}
#'
#' @examples
#' # O'Brien-Fleming type efficacy bounds at 50, 75 and 100 percent
#' gs_bounds(c(0.5, 0.75, 1))
#'
#' # Non-binding futility by beta spending with a delayed effect
#' gs_bounds(c(0.5, 0.75, 1), futility = "beta", sf_fut = "HSD",
#'           param_fut = -2, mean = c(0.8, 1.9, 2.8))
#'
#' @export
gs_bounds <- function(info_frac, alpha = 0.025, sf_eff = "OF", param_eff = NULL,
                      futility = c("none", "fixed", "beta"), fut_bound = NULL,
                      sf_fut = "OF", param_fut = NULL, beta = 0.1,
                      binding = FALSE, mean = NULL, corr = NULL,
                      sigma_alt = NULL) {
  futility <- match.arg(futility)
  .check_info_frac(info_frac)
  .check_prob(alpha, "alpha")
  if (!is.logical(binding) || length(binding) != 1L || is.na(binding)) {
    stop("'binding' must be TRUE or FALSE.", call. = FALSE)
  }
  k_max <- length(info_frac)
  if (is.null(corr)) {
    corr <- sqrt(outer(info_frac, info_frac, pmin) / outer(info_frac, info_frac, pmax))
  }
  .check_sigma(corr, k_max, "corr")
  a_cum <- gs_spend(info_frac, alpha, sf_eff, param_eff)
  a_inc <- diff(c(0, a_cum))
  lower <- rep(-Inf, k_max)
  upper <- rep(Inf, k_max)
  b_cum <- rep(NA_real_, k_max)
  if (futility == "fixed") {
    if (!is.numeric(fut_bound) || length(fut_bound) != k_max - 1L || anyNA(fut_bound)) {
      stop("'fut_bound' must have one value per interim look.", call. = FALSE)
    }
    lower[seq_len(k_max - 1L)] <- fut_bound
  }
  if (futility == "beta") {
    if (!is.numeric(mean) || length(mean) != k_max || anyNA(mean)) {
      stop("'mean' must have one value per look when futility = \"beta\".",
           call. = FALSE)
    }
    if (is.null(sigma_alt)) {
      sigma_alt <- corr
    }
    .check_sigma(sigma_alt, k_max, "sigma_alt")
    .check_prob(beta, "beta")
    b_cum <- gs_spend(info_frac, beta, sf_fut, param_fut)
  }
  b_inc <- diff(c(0, b_cum))
  for (k in seq_len(k_max)) {
    prev <- seq_len(k - 1L)
    lp <- if (binding) lower[prev] else rep(-Inf, k - 1L)
    upper[k] <- .solve_upper(k, lp, upper[prev], a_inc[k], corr)
    if (k < k_max) {
      if (futility == "beta") {
        lower[k] <- .solve_lower(k, lower[prev], upper[seq_len(k)], b_inc[k],
                                 mean, sigma_alt)
      }
      if (lower[k] > upper[k]) {
        lower[k] <- upper[k]
      }
    }
  }
  lower[k_max] <- upper[k_max]
  data.frame(look = seq_len(k_max), info_frac = info_frac, alpha_spent = a_cum,
             beta_spent = b_cum, upper = upper, lower = lower,
             nominal_p = stats::pnorm(upper, lower.tail = FALSE))
}
