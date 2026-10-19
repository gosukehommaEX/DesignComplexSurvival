#' Products of the PFS and OS functionals
#'
#' Per-patient E[r_P(t_p) r_O(t_o)] summed over the groups with the
#' allocation probabilities (see \code{.cross_group()}).
#'
#' @param scenario An \code{ms_scenario} object.
#' @param acc Output of \code{.acc_setup()}.
#' @param t_p,t_o Calendar times of the PFS and OS functionals.
#' @return A 3 x 3 matrix with rows for the PFS functionals and columns for
#'   the OS functionals (\code{psi}, \code{phi}, \code{dlt}).
#' @noRd
.cross_moments <- function(scenario, acc, t_p, t_o) {
  keys <- c("psi", "phi", "dlt")
  out <- matrix(0, 3L, 3L, dimnames = list(keys, keys))
  if (!(t_p > 0 && t_o > 0)) {
    return(out)
  }
  brk <- .scn_breaks(scenario)
  hmax <- .scn_hmax(scenario, max(t_p, t_o, 1))
  p <- scenario$prob
  for (grp in 0:1) {
    out <- out + p[[grp + 1L]] * .cross_group(scenario, grp, acc, t_p, t_o, brk, hmax)
  }
  out
}
