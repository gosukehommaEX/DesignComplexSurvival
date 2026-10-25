#' Rejection probabilities of a group sequential design with two endpoints
#'
#' For the standardized statistics of two endpoints with a joint multivariate
#' normal distribution, computes under a multiple testing procedure the
#' probability of rejecting each null hypothesis at each look, the
#' probability of rejecting both and the probability of rejecting at least
#' one.
#'
#' @details
#' Each endpoint is tested at its own looks, and \code{analysis} places the
#' looks of both endpoints on one sequence of analyses; an analysis may
#' include looks of one endpoint or of both. A hypothesis is rejected at the
#' first look at which its statistic is at or above its efficacy bound and is
#' not tested again. The trial stops for futility at an interim look of the
#' first endpoint at which its hypothesis has not been rejected and its
#' statistic is below the futility bound. The futility bound at the last look
#' of the first endpoint is not used.
#'
#' With \code{mtp = "bonferroni"}, each endpoint is tested at all its looks
#' until the trial stops, with bounds computed at levels that add up to
#' alpha. The trial stops when both hypotheses are rejected or for futility.
#'
#' With \code{mtp = "bonferroni"} and \code{upper_recycled}, the level of a
#' rejected hypothesis is passed on to the other hypothesis (the graphical
#' procedure of Maurer and Bretz 2013, Algorithm 1, for two hypotheses).
#' Each hypothesis is tested with its bounds in \code{upper} until the other
#' hypothesis is rejected, and with its bounds in \code{upper_recycled} from
#' the analysis at which the other is rejected on: a hypothesis rejected with
#' its bound in \code{upper} at an analysis makes the other be tested with its
#' bound in \code{upper_recycled} at the same analysis. The bounds in
#' \code{upper_recycled} are those of the level after recycling, for example
#' \code{(w1 + g21 w2) alpha} for the first hypothesis with the weights w and
#' the fraction g21 of the level of the second hypothesis that is passed on
#' to the first. They must not exceed the bounds in \code{upper}, which holds
#' when the nominal levels of the spending function do not decrease with the
#' level (Maurer and Bretz 2013, Section 3.1). The futility stop of the first
#' endpoint also applies after the second hypothesis has been rejected.
#'
#' With \code{mtp = "fixed_sequence"}, the second hypothesis is tested only
#' after the first has been rejected. With \code{hierarchical = "overall"} it
#' is tested at the analysis at which the first hypothesis is rejected and at
#' the later looks of the second endpoint, until it is rejected (strategy (b)
#' of Glimm et al. 2010). With \code{hierarchical = "stagewise"} the trial
#' stops at the analysis at which the first hypothesis is rejected, and the
#' second hypothesis is tested only at that analysis (strategy (a) of Glimm
#' et al. 2010, the procedure of Tamhane et al. 2010).
#'
#' The probabilities are sums of rectangle probabilities of the multivariate
#' normal distribution, one for each combination of the outcome of the first
#' endpoint (rejection or stop for futility at a look, or no rejection) and
#' the look at which the second hypothesis is rejected. With
#' \code{upper_recycled}, the probabilities are summed over the outcomes at
#' each analysis: at an analysis at which both hypotheses are tested and
#' neither has been rejected, with the bounds A1, A2 in \code{upper} and B1,
#' B2 in \code{upper_recycled}, both hypotheses are rejected when Z1 >= A1 and
#' Z2 >= B2 or when B1 <= Z1 < A1 and Z2 >= A2, only the first when Z1 >= A1
#' and Z2 < B2, only the second when Z1 < B1 and Z2 >= A2, and neither when
#' Z1 < A1 and Z2 < A2; at an analysis at which one hypothesis is tested, it
#' is rejected or not. The futility stop splits the outcomes in which the
#' first hypothesis is tested and not rejected at an interim look. The
#' rectangle probabilities are computed
#' with \code{mvtnorm::pmvnorm()} (at most 20 looks in total): with the
#' Genz-Bretz algorithm for the rectangles that involve both endpoints and
#' have four or more coordinates and, as in \code{\link{gs_prob}}, for the
#' other rectangles of four or more coordinates whose covariance matrix is
#' not that of a Markov sequence and of nine or more coordinates with two
#' finite limits for some coordinate, for which the Miwa algorithm can be
#' inaccurate, and with the Miwa algorithm otherwise. The Genz-Bretz
#' algorithm uses a fixed seed, so that the results are reproducible and the
#' random number generator of the session is not changed, and an absolute
#' error tolerance of 1e-5 for each rectangle, which can be changed with
#' \code{options(DesignComplexSurvival.genz_abseps = ...)}.
#'
#' @param upper List of two numeric vectors: the efficacy bounds of the first
#'   and of the second endpoint at their looks (\code{Inf} for no efficacy
#'   test at a look).
#' @param lower Futility bounds of the first endpoint at its looks
#'   (\code{-Inf} for none), or \code{NULL} for no futility bounds.
#' @param mean Mean vector of the statistics, first the looks of the first
#'   endpoint and then those of the second (recycled if of length 1).
#' @param sigma Covariance matrix of the statistics in the same order, for
#'   example \code{sigma} of \code{\link{lr_cov}} with two endpoints.
#' @param analysis List of two integer vectors: the analysis at which each
#'   look of the first and of the second endpoint takes place (strictly
#'   increasing within an endpoint; equal numbers mean the same analysis).
#'   \code{NULL} means that both endpoints are tested at analyses 1 to K.
#' @param mtp \code{"fixed_sequence"} or \code{"bonferroni"}.
#' @param hierarchical \code{"overall"} or \code{"stagewise"} (used with
#'   \code{mtp = "fixed_sequence"}).
#' @param upper_recycled For \code{mtp = "bonferroni"}, an optional list of
#'   two numeric vectors: the efficacy bounds of the first and of the second
#'   endpoint at their looks at the level after the other hypothesis has been
#'   rejected (not above \code{upper}). \code{NULL} means no recycling.
#'
#' @return A list with
#' \describe{
#'   \item{\code{looks}}{A data frame with one row per look: \code{endpoint}
#'     (1 or 2), \code{look}, \code{analysis}, the probability of rejecting
#'     the hypothesis of the endpoint at the look (\code{prob_upper}), its
#'     cumulative sum (\code{cum_upper}) and the probability of stopping for
#'     futility at the look (\code{prob_lower}, 0 for the second endpoint and
#'     for the last look of the first).}
#'   \item{\code{reject}}{The probabilities of rejecting the first hypothesis
#'     (\code{first}), the second (\code{second}), both (\code{both}) and at
#'     least one (\code{any}).}
#' }
#'
#' @references
#' Glimm E, Maurer W, Bretz F (2010). Hierarchical testing of multiple
#' endpoints in group-sequential trials. Statistics in Medicine 29:219-228.
#'
#' Maurer W, Bretz F (2013). Multiple testing in group sequential trials
#' using graphical approaches. Statistics in Biopharmaceutical Research
#' 5:311-320.
#'
#' Tamhane AC, Mehta CR, Liu L (2010). Testing a primary and a secondary
#' endpoint in a group sequential design. Biometrics 66:1174-1184.
#'
#' @seealso \code{\link{gs_prob}}, \code{\link{gs_bounds}},
#'   \code{\link{lr_cov}}
#'
#' @examples
#' # Two looks at information fractions 0.5 and 1 for both endpoints and
#' # correlation 0.5 between the endpoints
#' tt <- c(0.5, 1)
#' r <- sqrt(outer(tt, tt, pmin) / outer(tt, tt, pmax))
#' sigma <- rbind(cbind(r, 0.5 * r), cbind(0.5 * r, r))
#' b1 <- gs_bounds(tt, sf_eff = "OF")$upper
#' b2 <- gs_bounds(tt, sf_eff = "Pocock")$upper
#' gs_prob_mtp(list(b1, b2), mean = c(2, 2.9, 1.5, 2.1), sigma = sigma)$reject
#' # Probability of rejecting the second hypothesis when it is true
#' gs_prob_mtp(list(b1, b2), mean = c(2, 2.9, 0, 0), sigma = sigma)$reject
#' # Bonferroni with the level alpha / 2 for each endpoint and the level
#' # alpha after the other hypothesis has been rejected (full recycling)
#' h1 <- gs_bounds(tt, alpha = 0.0125, sf_eff = "OF")$upper
#' h2 <- gs_bounds(tt, alpha = 0.0125, sf_eff = "Pocock")$upper
#' gs_prob_mtp(list(h1, h2), mean = c(2, 2.9, 1.5, 2.1), sigma = sigma,
#'             mtp = "bonferroni", upper_recycled = list(b1, b2))$reject
#'
#' @export
gs_prob_mtp <- function(upper, lower = NULL, mean = 0, sigma, analysis = NULL,
                        mtp = c("fixed_sequence", "bonferroni"),
                        hierarchical = c("overall", "stagewise"),
                        upper_recycled = NULL) {
  mtp <- match.arg(mtp)
  hierarchical <- match.arg(hierarchical)
  if (!is.list(upper) || length(upper) != 2L ||
      !all(vapply(upper, function(u) is.numeric(u) && length(u) >= 1L && !anyNA(u),
                  logical(1)))) {
    stop("'upper' must be a list of two numeric vectors of efficacy bounds.",
         call. = FALSE)
  }
  k1 <- length(upper[[1L]])
  k2 <- length(upper[[2L]])
  k <- k1 + k2
  if (k > 20L) {
    stop("The two endpoints may have at most 20 looks in total.", call. = FALSE)
  }
  u1 <- upper[[1L]]
  u2 <- upper[[2L]]
  if (is.null(lower)) {
    lower <- rep(-Inf, k1)
  }
  if (!is.numeric(lower) || length(lower) != k1 || anyNA(lower)) {
    stop("'lower' must be NULL or have one value per look of the first endpoint.",
         call. = FALSE)
  }
  lower[k1] <- -Inf
  if (any(lower > u1)) {
    stop("'lower' must not exceed the efficacy bounds of the first endpoint.",
         call. = FALSE)
  }
  if (length(mean) == 1L) {
    mean <- rep(mean, k)
  }
  if (!is.numeric(mean) || length(mean) != k || anyNA(mean)) {
    stop("'mean' must have length 1 or one value per look of both endpoints.",
         call. = FALSE)
  }
  .check_sigma(sigma, k, "sigma")
  if (is.null(analysis)) {
    if (k1 != k2) {
      stop("'analysis' must be given when the endpoints have different numbers of looks.",
           call. = FALSE)
    }
    analysis <- list(seq_len(k1), seq_len(k2))
  }
  ok <- is.list(analysis) && length(analysis) == 2L &&
    all(vapply(seq_len(2L), function(e) {
      a <- analysis[[e]]
      is.numeric(a) && length(a) == length(upper[[e]]) && !anyNA(a) &&
        all(a >= 1) && all(a == round(a)) && !any(diff(a) <= 0)
    }, logical(1)))
  if (!ok) {
    stop(paste0("'analysis' must be a list of two strictly increasing vectors of ",
                "positive integers, one value per look of each endpoint."), call. = FALSE)
  }
  a1 <- analysis[[1L]]
  a2 <- analysis[[2L]]
  if (!is.null(upper_recycled)) {
    if (mtp != "bonferroni") {
      stop("'upper_recycled' is used only with mtp = \"bonferroni\".", call. = FALSE)
    }
    ok <- is.list(upper_recycled) && length(upper_recycled) == 2L &&
      all(vapply(seq_len(2L), function(e) {
        u <- upper_recycled[[e]]
        is.numeric(u) && length(u) == length(upper[[e]]) && !anyNA(u)
      }, logical(1)))
    if (!ok) {
      stop(paste0("'upper_recycled' must be NULL or a list of two numeric vectors with ",
                  "one value per look of each endpoint."), call. = FALSE)
    }
    if (any(upper_recycled[[1L]] > u1) || any(upper_recycled[[2L]] > u2)) {
      stop("'upper_recycled' must not exceed the efficacy bounds in 'upper'.", call. = FALSE)
    }
  }

  # probability of a rectangle on the coordinates idx; coordinates without
  # finite limits are left out, and rectangles over both endpoints with four
  # or more coordinates use the Genz-Bretz algorithm (.pmvn() also uses it
  # for the other rectangles of four or more coordinates whose covariance is
  # not Markov)
  box <- function(idx, lo, hi) {
    if (any(lo >= hi)) {
      return(0)
    }
    keep <- is.finite(lo) | is.finite(hi)
    if (!any(keep)) {
      return(1)
    }
    i <- idx[keep]
    if (length(i) >= 4L && any(i <= k1) && any(i > k1)) {
      return(.pmvn_genz(lo[keep], hi[keep], mean[i], sigma[i, i, drop = FALSE]))
    }
    .pmvn(lo[keep], hi[keep], mean[i], sigma[i, i, drop = FALSE])
  }

  p1 <- numeric(k1)
  pf <- numeric(k1)
  p2 <- numeric(k2)
  both <- 0
  pats <- if (is.null(upper_recycled)) {
    .mtp_patterns(a1, a2, is.finite(lower), mtp, hierarchical)
  } else {
    list()
  }
  for (pat in pats) {
    if (pat$type == "none" && !length(pat$eligible)) {
      next
    }
    j <- pat$look
    prev <- seq_len(j - 1L)
    lo1 <- c(lower[prev], switch(pat$type, reject = u1[j], -Inf))
    hi1 <- c(u1[prev], switch(pat$type, reject = Inf, futility = lower[j], none = u1[j]))
    idx1 <- seq_len(j)
    pr <- box(idx1, lo1, hi1)
    if (pat$type == "reject") {
      p1[j] <- pr
    } else if (pat$type == "futility") {
      pf[j] <- pr
    }
    elig <- pat$eligible
    for (r in seq_along(elig)) {
      i <- elig[r]
      if (!is.finite(u2[i]) || pr == 0) {
        next
      }
      before <- elig[seq_len(r - 1L)]
      q <- box(c(idx1, k1 + before, k1 + i), c(lo1, rep(-Inf, r - 1L), u2[i]),
               c(hi1, u2[before], Inf))
      p2[i] <- p2[i] + q
      if (pat$type == "reject") {
        both <- both + q
      }
    }
  }
  if (!is.null(upper_recycled)) {
    gr <- .mtp_graph(list(u1, u2), upper_recycled, lower, list(a1, a2), box)
    p1 <- gr$p1
    pf <- gr$pf
    p2 <- gr$p2
    both <- gr$both
  }
  looks <- data.frame(endpoint = rep(1:2, c(k1, k2)), look = c(seq_len(k1), seq_len(k2)),
                      analysis = c(a1, a2), prob_upper = c(p1, p2),
                      cum_upper = c(cumsum(p1), cumsum(p2)), prob_lower = c(pf, numeric(k2)))
  first <- sum(p1)
  second <- sum(p2)
  list(looks = looks,
       reject = c(first = first, second = second, both = both,
                  any = first + second - both))
}
