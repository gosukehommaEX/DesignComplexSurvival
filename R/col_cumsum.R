#' Cumulative sums of each column of a matrix
#'
#' @param x Numeric matrix.
#' @return Matrix of the same dimension.
#' @noRd
.col_cumsum <- function(x) {
  x <- as.matrix(x)
  for (j in seq_len(ncol(x))) {
    x[, j] <- cumsum(x[, j])
  }
  x
}
