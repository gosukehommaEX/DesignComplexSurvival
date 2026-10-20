#' Type I error of a design under the null scenario
#'
#' Applies the bounds of the design to the statistics under the null
#' scenario. Event-driven looks take place when the same numbers of events
#' of the trigger endpoint are expected under the null scenario;
#' calendar-driven looks keep their calendar times.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param acc Output of \code{.acc_setup()}.
#' @param looks An \code{ms_looks} object.
#' @param targets Event targets of the trigger endpoints of the looks
#'   (\code{NA} for calendar-driven looks).
#' @param endpoint Tested endpoint.
#' @param bnd Output of \code{gs_bounds()}.
#' @return A list with \code{type1} (non-binding and binding) and the
#'   cumulative non-binding type I error \code{cum_nonbinding}.
#' @noRd
.type1_core <- function(scenario, acc, looks, targets, endpoint, bnd) {
  s0 <- ms_null(scenario)
  times0 <- .look_times(s0, acc, .looks_table(looks), targets = targets,
                        what = "the null scenario")$time
  core0 <- .lr_core(s0, acc, times0, endpoint, "direct")
  v0 <- core0$info
  corr <- sqrt(outer(v0, v0, pmin) / outer(v0, v0, pmax))
  k_max <- length(times0)
  nb <- gs_prob(rep(-Inf, k_max), bnd$upper, 0, corr)
  bd <- gs_prob(bnd$lower, bnd$upper, 0, corr)
  list(type1 = c(nonbinding = sum(nb$prob_upper), binding = sum(bd$prob_upper)),
       cum_nonbinding = nb$cum_upper)
}
