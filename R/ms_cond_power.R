#' Conditional power at an interim look
#'
#' Computes, for a design with one endpoint, the probability that the trial
#' crosses an efficacy bound at a later look, given the standardized log-rank
#' statistic observed at an interim look, under a chosen distribution of the
#' statistics.
#'
#' @details
#' At the interim look j with the observed statistic z, the statistics of the
#' later looks given \eqn{Z_j = z} have the conditional normal distribution
#' of the joint normal distribution of the statistics chosen by \code{under}.
#' The conditional power is the probability that, at some later look k, the
#' statistic crosses the efficacy bound after the statistics of the looks
#' between j and k have stayed between the futility and efficacy bounds, as
#' in \code{gsCP()} of gsDesign; \code{cond_power_nonbinding} leaves out the
#' futility bounds of the later interim looks (a trial that continues despite
#' non-binding futility bounds). The bounds are those of the design
#' (\code{x$looks}). When the events observed at the interim look differ from
#' those planned, the bounds of the design no longer apply: compute the design
#' again with \code{\link{ms_power}} with looks that give the observed events
#' and the remaining planned events, and use the new design.
#'
#' With \code{z} of length \code{look}, the statistics of the looks 1 to
#' \code{look} are all conditioned on. For statistics with independent
#' increments (\code{variance = "null"} of \code{\link{ms_test}}, except OS
#' with switching from a calendar time) the result is the same as with the
#' statistic of the look alone; with \code{"alternative"} and
#' \code{"studentized"} it can differ. In both cases the conditioning does not
#' use the fact that the trial did not stop at the earlier looks.
#'
#' The distributions (\code{under}):
#' \describe{
#'   \item{\code{"design"}}{The mean and covariance of the statistics under
#'     the design scenario (\code{x$mean}, \code{x$sigma}).}
#'   \item{\code{"null"}}{Mean zero and the correlations under the null
#'     scenario (\code{x$sigma_null}), at the looks of the type I error
#'     (\code{x$looks_null}); the result is the conditional type I error.}
#'   \item{\code{"trend"}}{The current trend: a constant hazard ratio
#'     \code{exp(-z / sqrt(p0 p1 d_j))} estimated from z, where p0 and p1 are
#'     the allocation probabilities and d the expected events of the tested
#'     endpoint at the looks of the design, so that the mean of the statistic
#'     at look k is \code{z sqrt(d_k / d_j)}, with the correlations
#'     \code{sqrt(d_j / d_k)} used for the bounds.}
#'   \item{\code{"hr"}}{A constant hazard ratio \code{hr}: the mean at look k
#'     is \code{-log(hr) sqrt(p0 p1 d_k)}, with the correlations
#'     \code{sqrt(d_j / d_k)}.}
#'   \item{\code{"scenario"}}{The mean and covariance under \code{scenario},
#'     with the approximations of the test of the design; the event-driven
#'     looks take place when the expected events of their trigger endpoint
#'     reach the event targets of the design (\code{x$targets}), also when the
#'     final look of the design is fixed by \code{followup}, and the
#'     calendar-driven looks keep their times.}
#' }
#' \code{"trend"} and \code{"hr"} assume proportional hazards (the Schoenfeld
#' approximation) and do not use the shape of the hazard ratio of the design
#' scenario, such as a delayed effect; \code{"scenario"} does.
#'
#' @param x An \code{ms_power} object with one endpoint (from
#'   \code{\link{ms_power}} or \code{\link{ms_size}}).
#' @param look Interim look at which \code{z} is observed (from 1 to the
#'   number of looks minus 1).
#' @param z Observed standardized log-rank statistic at \code{look}
#'   (positive when the treatment is better), or the observed statistics at
#'   the looks 1 to \code{look}. It must be below the efficacy bounds and, with
#'   binding futility bounds, not below the futility bounds.
#' @param under Distribution of the statistics: \code{"design"},
#'   \code{"null"}, \code{"trend"}, \code{"hr"} or \code{"scenario"} (see
#'   Details).
#' @param hr Hazard ratio (treatment versus control) for \code{under = "hr"}.
#' @param scenario An \code{ms_scenario} object with the allocation and
#'   accrual of the design scenario, for \code{under = "scenario"}.
#'
#' @return An object of class \code{ms_cond_power}, a list with
#' \describe{
#'   \item{cond_power}{The conditional power.}
#'   \item{cond_power_nonbinding}{The conditional power without the futility
#'     bounds of the later interim looks.}
#'   \item{looks}{A data frame with one row per later look: \code{look},
#'     \code{time} and \code{events} (calendar time and expected events of the
#'     tested endpoint under the distribution of \code{under}), the bounds
#'     \code{z_upper} and \code{z_lower}, the conditional mean and standard
#'     deviation of the statistic (\code{cond_mean}, \code{cond_sd}), the
#'     conditional probabilities of stopping for efficacy and for futility at
#'     the look (\code{prob_upper}, \code{prob_lower}) and the cumulative
#'     conditional power \code{cum_upper}.}
#'   \item{look, z, under}{The inputs.}
#'   \item{hr}{The hazard ratio of \code{"trend"} and \code{"hr"}
#'     (\code{NA} otherwise).}
#' }
#'
#' @references
#' Ortega-Villa AM, Grieco MC, Rubenstein K, Wang J, Proschan MA (2025).
#' Futility monitoring in clinical trials. Statistics in Medicine, 44,
#' e70157.
#'
#' @seealso \code{\link{ms_power}}, \code{\link{gs_prob}}
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' des <- ms_power(scn, n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
#'                 test = ms_test(), events = 300)
#' ms_cond_power(des, look = 1, z = 1.2)
#' ms_cond_power(des, look = 1, z = 1.2, under = "trend")
#' ms_cond_power(des, look = 1, z = 1.2, under = "hr", hr = 0.8)
#'
#' @export
ms_cond_power <- function(x, look, z, under = c("design", "null", "trend", "hr", "scenario"),
                          hr = NULL, scenario = NULL) {
  if (!inherits(x, "ms_power")) {
    stop("'x' must be an 'ms_power' object (see ms_power()).", call. = FALSE)
  }
  if (length(x$settings$endpoint) != 1L) {
    stop("ms_cond_power() supports designs with one endpoint only.", call. = FALSE)
  }
  if (is.null(x[["sigma"]]) || is.null(x[["sigma_null"]])) {
    stop("'x' has no covariance matrices of the statistics; compute it again with ms_power().",
         call. = FALSE)
  }
  under <- match.arg(under)
  lk <- x$looks
  k_max <- nrow(lk)
  if (k_max < 2L) {
    stop("The design has no interim look.", call. = FALSE)
  }
  if (!is.numeric(look) || length(look) != 1L || is.na(look) || look != round(look) ||
      look < 1 || look > k_max - 1L) {
    stop(sprintf("'look' must be an interim look, an integer from 1 to %d.", k_max - 1L),
         call. = FALSE)
  }
  look <- as.integer(look)
  if (!is.numeric(z) || !(length(z) %in% c(1L, look)) || any(!is.finite(z))) {
    stop(sprintf("'z' must be the statistic at look %d or the statistics at looks 1 to %d.",
                 look, look), call. = FALSE)
  }
  obs <- if (length(z) == 1L) look else seq_len(look)
  if (any(z >= lk$z_upper[obs])) {
    stop("'z' crosses the efficacy bound, at which the trial stops.", call. = FALSE)
  }
  if (x$settings$binding && any(z < lk$z_lower[obs])) {
    stop("'z' falls below the binding futility bound, at which the trial stops.", call. = FALSE)
  }
  if (under == "hr") {
    if (!is.numeric(hr) || length(hr) != 1L || !is.finite(hr) || !(hr > 0)) {
      stop("'hr' must be a positive number for under = \"hr\".", call. = FALSE)
    }
  } else if (!is.null(hr)) {
    stop("'hr' is used only with under = \"hr\".", call. = FALSE)
  }
  if (under == "scenario") {
    .check_scenario(scenario)
    if (!isTRUE(all.equal(scenario$prob, x$scenario$prob)) ||
        !isTRUE(all.equal(scenario$accrual, x$scenario$accrual))) {
      stop("'scenario' must have the allocation and accrual of the design scenario.",
           call. = FALSE)
    }
    .check_calendar(scenario, x$settings$endpoint, x$test$mean, x$test$variance,
                    .looks_table(x$looks_def)$trigger)
  } else if (!is.null(scenario)) {
    stop("'scenario' is used only with under = \"scenario\".", call. = FALSE)
  }
  mo <- .cp_moments(x, look, z[length(z)], under, hr = hr, scenario = scenario)
  cp <- .cond_probs(mo$mean, mo$sigma, look, z, lk$z_lower, lk$z_upper)
  fut <- seq.int(look + 1L, k_max)
  looks <- data.frame(look = fut, time = mo$time[fut], events = mo$events[fut],
                      z_upper = lk$z_upper[fut], z_lower = lk$z_lower[fut],
                      cp$looks[-1L])
  structure(list(cond_power = cp$cond_power,
                 cond_power_nonbinding = cp$cond_power_nonbinding, looks = looks,
                 look = look, z = z, under = under, hr = mo$hr),
            class = "ms_cond_power")
}
