#' Distribution of the statistics for the conditional power
#'
#' The mean vector and covariance matrix of the statistics at the looks of a
#' design with one endpoint under the distribution chosen in
#' \code{ms_cond_power()}, with the calendar times and the expected events of
#' the looks.
#'
#' @param x An \code{ms_power} object with one endpoint.
#' @param look Interim look of the observation.
#' @param z Observed statistic at \code{look}.
#' @param under \code{"design"}, \code{"null"}, \code{"trend"}, \code{"hr"} or
#'   \code{"scenario"}.
#' @param hr Hazard ratio (for \code{"hr"}).
#' @param scenario An \code{ms_scenario} object with the allocation and
#'   accrual of the design (for \code{"scenario"}).
#' @return A list with \code{mean}, \code{sigma}, \code{time}, \code{events}
#'   and the hazard ratio \code{hr} of \code{"trend"} and \code{"hr"}
#'   (\code{NA} otherwise).
#' @details With \code{"trend"} and \code{"hr"} the hazard ratio is constant:
#'   the mean of the statistic at look k is \code{-log(hr) sqrt(p0 p1 d_k)},
#'   where d are the expected events of the design scenario, and the
#'   correlations are \code{sqrt(d_j / d_k)}, those used for the bounds. The
#'   hazard ratio of \code{"trend"} is \code{exp(-z / sqrt(p0 p1 d_look))}, so
#'   that the mean at the look is \code{z}. With \code{"scenario"} the
#'   event-driven looks take place when the expected events of their trigger
#'   endpoint reach the targets of the design (\code{x$targets}) and the
#'   calendar-driven looks keep their times, and the moments follow the
#'   approximations of \code{x$test}.
#' @noRd
.cp_moments <- function(x, look, z, under, hr = NULL, scenario = NULL) {
  lk <- x$looks
  d <- lk$events
  p <- x$scenario$prob
  corr_d <- sqrt(outer(d, d, pmin) / outer(d, d, pmax))
  switch(under,
    design = list(mean = x$mean, sigma = x$sigma, time = lk$time, events = d,
                  hr = NA_real_),
    null = list(mean = rep(0, length(d)), sigma = x$sigma_null,
                time = x$looks_null$time, events = x$looks_null$events, hr = NA_real_),
    trend = list(mean = z * sqrt(d / d[look]), sigma = corr_d, time = lk$time, events = d,
                 hr = exp(-z / sqrt(p[[1L]] * p[[2L]] * d[look]))),
    hr = list(mean = -log(hr) * sqrt(p[[1L]] * p[[2L]] * d), sigma = corr_d, time = lk$time,
              events = d, hr = hr),
    scenario = {
      ep <- x$settings$endpoint
      acc <- .acc_setup(scenario$accrual, x$n)
      lt <- .look_times(scenario, acc, .looks_table(x$looks_def), targets = x$targets,
                        what = "the given scenario")
      core <- .lr_core(scenario, acc, lt$time, ep, x$test$mean)
      list(mean = core$z_mean,
           sigma = .z_cov(scenario, acc, lt$time, ep, core, x$test$variance, lt$trigger),
           time = lt$time, events = core$events, hr = NA_real_)
    })
}
