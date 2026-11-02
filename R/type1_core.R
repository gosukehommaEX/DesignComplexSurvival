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
#' @return A list with \code{type1} (non-binding and binding), the
#'   cumulative non-binding type I error \code{cum_nonbinding}, a data frame
#'   \code{looks} with the calendar \code{time}, the \code{enrolled} patients
#'   and the \code{events} of the endpoint at the looks under the null
#'   scenario and the probabilities of stopping for efficacy and for futility
#'   at the looks with the futility bounds applied (\code{prob_upper},
#'   \code{prob_lower}), and the correlation matrix \code{sigma} of the
#'   statistics under the null scenario.
#' @details The statistics under the null scenario have unit variances and
#'   the correlations sqrt(V_j / V_k), except for OS with switching from a
#'   calendar time (\code{.is_calendar()}), where the score has no
#'   independent increments and the correlations are those of
#'   \code{.cal_zeta()}; the type I error then differs from the level of the
#'   bounds.
#' @noRd
.type1_core <- function(scenario, acc, looks, targets, endpoint, bnd) {
  s0 <- ms_null(scenario)
  times0 <- .look_times(s0, acc, .looks_table(looks), targets = targets,
                        what = "the null scenario")$time
  core0 <- .lr_core(s0, acc, times0, endpoint, "direct")
  v0 <- core0$info
  corr <- if (endpoint == "os" && .is_calendar(s0)) {
    .cal_zeta(s0, acc, times0, NULL, "null")$sigma
  } else {
    sqrt(outer(v0, v0, pmin) / outer(v0, v0, pmax))
  }
  k_max <- length(times0)
  nb <- gs_prob(rep(-Inf, k_max), bnd$upper, 0, corr)
  bd <- gs_prob(bnd$lower, bnd$upper, 0, corr)
  list(type1 = c(nonbinding = sum(nb$prob_upper), binding = sum(bd$prob_upper)),
       cum_nonbinding = nb$cum_upper,
       looks = data.frame(look = seq_len(k_max), time = times0,
                          enrolled = acc$n * .acc_frac(acc, times0),
                          events = core0$events, prob_upper = bd$prob_upper,
                          prob_lower = bd$prob_lower),
       sigma = corr)
}
