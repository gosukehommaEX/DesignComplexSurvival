#' Gauss-Legendre nodes and weights on the interval from -1 to 1
#'
#' Uses the Golub-Welsch algorithm and caches the result.
#'
#' @param n Number of nodes.
#' @return A list with nodes \code{x} (increasing) and weights \code{w}.
#' @noRd
.gl_rule <- function(n) {
  key <- paste0("gl", n)
  if (!is.null(.dcs_env[[key]])) {
    return(.dcs_env[[key]])
  }
  if (n == 1L) {
    out <- list(x = 0, w = 2)
  } else {
    k <- seq_len(n - 1L)
    b <- k / sqrt(4 * k^2 - 1)
    jac <- matrix(0, n, n)
    jac[cbind(k, k + 1L)] <- b
    jac[cbind(k + 1L, k)] <- b
    e <- eigen(jac, symmetric = TRUE)
    ord <- order(e$values)
    out <- list(x = e$values[ord], w = 2 * e$vectors[1L, ord]^2)
  }
  assign(key, out, envir = .dcs_env)
  out
}
