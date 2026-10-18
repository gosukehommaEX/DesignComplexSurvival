#' Composite Gauss-Legendre rule on an interval
#'
#' Splits the interval from \code{a} to \code{b} at the given break points
#' and into pieces no longer than \code{hmax}, and places \code{n}
#' Gauss-Legendre nodes in each piece.
#'
#' @param a,b Interval limits.
#' @param breaks Break points (those outside the open interval are ignored).
#' @param hmax Maximum piece length.
#' @param n Number of nodes per piece.
#' @return A list with nodes \code{x}, weights \code{w}, the piece index of
#'   each node \code{piece}, and the piece limits \code{lo} and \code{hi}.
#' @noRd
.gl_pieces <- function(a, b, breaks = numeric(0), hmax = Inf, n = 12L) {
  if (!(b > a)) {
    return(list(x = numeric(0), w = numeric(0), piece = integer(0),
                lo = numeric(0), hi = numeric(0)))
  }
  pts <- sort(unique(c(a, breaks[breaks > a & breaks < b], b)))
  if (is.finite(hmax)) {
    len <- diff(pts)
    m <- pmax(1L, as.integer(ceiling(len / hmax - 1e-9)))
    if (any(m > 1L)) {
      sub <- lapply(seq_along(len), function(i) {
        pts[i] + (seq_len(m[i]) - 1L) * len[i] / m[i]
      })
      pts <- c(unlist(sub, use.names = FALSE), b)
    }
  }
  rule <- .gl_rule(n)
  lo <- pts[-length(pts)]
  hi <- pts[-1L]
  half <- (hi - lo) / 2
  mid <- (hi + lo) / 2
  list(x = as.vector(outer(rule$x, half)) + rep(mid, each = n),
       w = as.vector(outer(rule$w, half)),
       piece = rep(seq_along(lo), each = n),
       lo = lo, hi = hi)
}
