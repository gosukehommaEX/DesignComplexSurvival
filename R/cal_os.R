#' OS of one group given the follow-up time from which switching is possible
#'
#' \code{.os_avail()} for a group that can switch, and \code{.ep_eval()} (no
#' dependence on \code{o}) for a group that cannot.
#'
#' @param model An \code{ms_model} object.
#' @param s,o Numeric vectors of equal length.
#' @return A list with \code{surv} and \code{dens}.
#' @noRd
.cal_os <- function(model, s, o) {
  if (.has_switch(model)) {
    return(.os_avail(model, s, o))
  }
  ep <- .ep_eval(model, s, "os")
  list(surv = ep$surv, dens = ep$dens)
}
