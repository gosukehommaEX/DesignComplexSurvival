#' The function (exp(z) - 1) / z, equal to 1 at z = 0
#'
#' @param z Numeric vector.
#' @return Numeric vector.
#' @noRd
.expm1_ratio <- function(z) {
  out <- rep(1, length(z))
  big <- abs(z) >= 1e-8
  out[big] <- expm1(z[big]) / z[big]
  out[!big] <- 1 + z[!big] / 2
  out
}
