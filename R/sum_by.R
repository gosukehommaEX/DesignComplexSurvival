#' Column sums by group index
#'
#' @param x Numeric vector or matrix.
#' @param idx Integer group index of each row of \code{x}.
#' @param n Number of groups.
#' @return An \code{n} by \code{ncol(x)} matrix (groups without rows are 0).
#' @noRd
.sum_by <- function(x, idx, n) {
  x <- as.matrix(x)
  out <- matrix(0, n, ncol(x))
  if (length(idx) > 0L) {
    r <- rowsum(x, idx, reorder = TRUE)
    out[as.integer(rownames(r)), ] <- r
  }
  out
}
