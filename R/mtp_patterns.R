#' Outcomes of the first endpoint and looks of the second
#'
#' Lists the outcomes of the first endpoint (rejection at a look, stop for
#' futility at an interim look, or no rejection) and, for each outcome, the
#' looks at which the second hypothesis is tested until it is rejected.
#'
#' @param a1,a2 Analyses of the looks of the first and of the second
#'   endpoint.
#' @param fut Logical vector, one value per look of the first endpoint:
#'   whether it has a futility bound (the last value is not used).
#' @param mtp \code{"fixed_sequence"} or \code{"bonferroni"}.
#' @param hierarchical \code{"overall"} or \code{"stagewise"}.
#' @return A list with one element per outcome: \code{type}
#'   (\code{"reject"}, \code{"futility"} or \code{"none"}), \code{look} (the
#'   look of the first endpoint at which the outcome occurs; the last look
#'   for \code{"none"}) and \code{eligible} (looks of the second endpoint).
#' @noRd
.mtp_patterns <- function(a1, a2, fut, mtp, hierarchical) {
  k1 <- length(a1)
  all2 <- seq_along(a2)
  out <- list()
  for (j in seq_len(k1)) {
    elig <- if (mtp == "bonferroni") {
      all2
    } else if (hierarchical == "overall") {
      all2[a2 >= a1[j]]
    } else {
      all2[a2 == a1[j]]
    }
    out[[length(out) + 1L]] <- list(type = "reject", look = j, eligible = elig)
  }
  for (j in seq_len(k1 - 1L)) {
    if (fut[j]) {
      elig <- if (mtp == "bonferroni") all2[a2 <= a1[j]] else integer(0)
      out[[length(out) + 1L]] <- list(type = "futility", look = j, eligible = elig)
    }
  }
  elig <- if (mtp == "bonferroni") all2 else integer(0)
  out[[length(out) + 1L]] <- list(type = "none", look = k1, eligible = elig)
  out
}
