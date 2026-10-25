#' Rejection probabilities of the graphical procedure with two endpoints
#'
#' Sums the rectangle probabilities of the outcomes of the Bonferroni-based
#' graphical procedure with recycling of the level (Maurer and Bretz 2013,
#' Algorithm 1) for two hypotheses, analysis by analysis.
#'
#' A hypothesis is tested with its bound at the initial level until the
#' other hypothesis is rejected, and with its bound at the level after
#' recycling afterwards. At an analysis at which both hypotheses are tested
#' and neither has been rejected, with the bounds A1, A2 at the initial level
#' and B1 <= A1, B2 <= A2 after recycling, the outcomes are the rectangles
#' \{Z1 >= A1, Z2 >= B2\} and \{B1 <= Z1 < A1, Z2 >= A2\} (both rejected),
#' \{Z1 >= A1, Z2 < B2\} (first only), \{Z1 < B1, Z2 >= A2\} (second only)
#' and \{Z1 < A1, Z2 < A2\} (neither). At an analysis at which one hypothesis
#' is tested, the outcomes are rejection and no rejection. The futility stop
#' of the first endpoint splits the outcomes in which its hypothesis is
#' tested and not rejected at an interim look. The trial stops when both
#' hypotheses are rejected or for futility.
#'
#' @param upper List of the bounds of the two endpoints at the initial level.
#' @param recycled List of the bounds after recycling (not above
#'   \code{upper}).
#' @param lower Futility bounds of the first endpoint (\code{-Inf} for none).
#' @param analysis List of the analyses of the looks of the two endpoints.
#' @param box Function of the coordinates and the lower and upper limits that
#'   returns the probability of a rectangle (from \code{gs_prob_mtp()}).
#' @return A list with the probabilities of rejecting the first hypothesis
#'   (\code{p1}) and of stopping for futility (\code{pf}) at each look of the
#'   first endpoint, of rejecting the second hypothesis at each look of the
#'   second endpoint (\code{p2}), and of rejecting both (\code{both}).
#' @noRd
.mtp_graph <- function(upper, recycled, lower, analysis, box) {
  k1 <- length(upper[[1L]])
  an <- sort(unique(c(analysis[[1L]], analysis[[2L]])))
  lower[k1] <- -Inf
  prob <- list(numeric(k1), numeric(length(upper[[2L]])))
  fut <- numeric(k1)
  both <- 0
  coord <- function(e, j) if (e == 1L) j else k1 + j
  # probability of the rectangle of the limits fixed so far (cons) and new
  # limits on the coordinates idx
  rect <- function(cons, idx, lo, hi) {
    i <- c(cons$idx, idx)
    o <- order(i)
    box(i[o], c(cons$lo, lo)[o], c(cons$hi, hi)[o])
  }
  extend <- function(cons, idx, lo, hi) {
    list(idx = c(cons$idx, idx), lo = c(cons$lo, lo), hi = c(cons$hi, hi))
  }
  # outcomes from analysis ai on, given the hypotheses rejected before (rej)
  step <- function(ai, rej, cons) {
    if (ai > length(an) || all(rej)) {
      return(invisible(NULL))
    }
    lv <- c(match(an[ai], analysis[[1L]]), match(an[ai], analysis[[2L]]))
    tested <- which(!is.na(lv) & !rej)
    if (length(tested) == 0L) {
      return(step(ai + 1L, rej, cons))
    }
    if (length(tested) == 1L) {
      e <- tested
      j <- lv[e]
      b <- if (rej[3L - e]) recycled[[e]][j] else upper[[e]][j]
      ce <- coord(e, j)
      pr <- rect(cons, ce, b, Inf)
      prob[[e]][j] <<- prob[[e]][j] + pr
      if (rej[3L - e]) {
        both <<- both + pr
      }
      if (b < Inf) {
        nr <- rej
        nr[e] <- TRUE
        step(ai + 1L, nr, extend(cons, ce, b, Inf))
      }
      return(go_on(ai, rej, cons, ce, -Inf, b, if (e == 1L) j else NA_integer_))
    }
    j1 <- lv[1L]
    j2 <- lv[2L]
    a1 <- upper[[1L]][j1]
    b1 <- recycled[[1L]][j1]
    a2 <- upper[[2L]][j2]
    b2 <- recycled[[2L]][j2]
    idx <- c(coord(1L, j1), coord(2L, j2))
    pr <- rect(cons, idx, c(a1, b2), c(Inf, Inf)) + rect(cons, idx, c(b1, a2), c(a1, Inf))
    prob[[1L]][j1] <<- prob[[1L]][j1] + pr
    prob[[2L]][j2] <<- prob[[2L]][j2] + pr
    both <<- both + pr
    pr <- rect(cons, idx, c(a1, -Inf), c(Inf, b2))
    prob[[1L]][j1] <<- prob[[1L]][j1] + pr
    if (a1 < Inf) {
      step(ai + 1L, c(TRUE, FALSE), extend(cons, idx, c(a1, -Inf), c(Inf, b2)))
    }
    pr <- rect(cons, idx, c(-Inf, a2), c(b1, Inf))
    prob[[2L]][j2] <<- prob[[2L]][j2] + pr
    if (a2 < Inf) {
      go_on(ai, c(FALSE, TRUE), cons, idx, c(-Inf, a2), c(b1, Inf), j1)
    }
    go_on(ai, c(FALSE, FALSE), cons, idx, c(-Inf, -Inf), c(a1, a2), j1)
  }
  # continue to the next analysis with the limits lo, hi on the coordinates
  # idx; when the first hypothesis is tested at its look j1 and not rejected,
  # the trial stops below its futility bound
  go_on <- function(ai, rej, cons, idx, lo, hi, j1) {
    if (!is.na(j1) && is.finite(lower[j1])) {
      p <- which(idx == coord(1L, j1))
      hs <- hi
      hs[p] <- min(hi[p], lower[j1])
      fut[j1] <<- fut[j1] + rect(cons, idx, lo, hs)
      lo[p] <- max(lo[p], lower[j1])
    }
    if (all(lo < hi)) {
      step(ai + 1L, rej, extend(cons, idx, lo, hi))
    }
    invisible(NULL)
  }
  step(1L, c(FALSE, FALSE), list(idx = integer(0), lo = numeric(0), hi = numeric(0)))
  list(p1 = prob[[1L]], pf = fut, p2 = prob[[2L]], both = both)
}
