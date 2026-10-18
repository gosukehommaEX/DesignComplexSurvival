#' Maximum quadrature piece length for a scenario
#'
#' Pieces are short enough that the integrands change by a factor of about
#' exp(2) at most within a piece, with at most about 400 pieces over the
#' range of integration.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param span Length of the range of integration.
#' @return Maximum piece length.
#' @noRd
.scn_hmax <- function(scenario, span) {
  lam <- unlist(lapply(list(scenario$control, scenario$treatment), function(x) {
    c(x$prog + x$death + x$dropout, x$switch_rate + x$death_prog, x$death_switch)
  }))
  max(2 / max(lam, 1e-12), span / 400)
}
