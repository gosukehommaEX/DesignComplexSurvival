#' Operating characteristics of a design
#'
#' Determines the analysis times, the moments of the tested statistic, the
#' bounds and the stopping probabilities under the design scenario, and
#' optionally the type I error under the null scenario (one endpoint).
#'
#' @param scenario An \code{ms_scenario} object.
#' @param n Total sample size.
#' @param looks An \code{ms_looks} object.
#' @param test An \code{ms_test} object with one endpoint.
#' @param events Final number of events of the trigger endpoint (named by
#'   endpoint when the looks are driven by the events of two endpoints).
#' @param followup Time from the end of accrual to the final analysis.
#' @param beta Total type II error for beta spending (defaults to
#'   \code{test$beta}).
#' @param type1 Whether to compute the type I error.
#' @return A list (see \code{\link{ms_power}}).
#' @noRd
.design_core <- function(scenario, n, looks, test, events = NULL,
                         followup = NULL, beta = NULL, type1 = TRUE) {
  acc <- .acc_setup(scenario$accrual, n)
  k_max <- looks$k
  ep <- test$endpoints
  .tested_looks(looks, ep)
  lk <- .looks_table(looks)
  lt <- .look_times(scenario, acc, lk, events = events, followup = followup,
                    what = "the design scenario")
  times <- lt$time
  targets <- lt$target
  core <- .lr_core(scenario, acc, times, ep, test$mean)
  sigma <- .z_cov(scenario, acc, times, ep, core, test$variance, lt$trigger)
  d <- core$events
  t_frac <- if (test$spending_time == "events") d / d[k_max] else core$info / core$info[k_max]
  corr0 <- sqrt(outer(d, d, pmin) / outer(d, d, pmax))
  p <- scenario$prob
  root_info <- sqrt(p[[1L]] * p[[2L]] * d)
  fut <- NULL
  if (test$futility == "fixed") {
    fb <- test$fut_bound
    if (length(fb) != k_max - 1L) {
      stop("'fut_bound' must have one value per interim look (", k_max - 1L, ").",
           call. = FALSE)
    }
    fut <- if (test$fut_scale == "z") fb else -log(fb) * root_info[seq_len(k_max - 1L)]
  }
  if (test$sf_eff == "user" && length(test$param_eff) != k_max) {
    stop("'param_eff' must have one value per look for sf_eff = \"user\".", call. = FALSE)
  }
  if (test$futility == "beta" && test$sf_fut == "user" && length(test$param_fut) != k_max) {
    stop("'param_fut' must have one value per look for sf_fut = \"user\".", call. = FALSE)
  }
  bnd <- gs_bounds(t_frac, alpha = test$alpha, sf_eff = test$sf_eff,
                   param_eff = test$param_eff, futility = test$futility,
                   fut_bound = fut, sf_fut = test$sf_fut,
                   param_fut = test$param_fut,
                   beta = if (is.null(beta)) test$beta else beta,
                   binding = test$binding, mean = core$z_mean, corr = corr0,
                   sigma_alt = sigma)
  pr <- gs_prob(bnd$lower, bnd$upper, core$z_mean, sigma)
  hr_lower <- exp(-bnd$lower / root_info)
  hr_lower[is.infinite(bnd$lower)] <- NA_real_
  out <- data.frame(look = seq_len(k_max), time = times,
                    enrolled = n * .acc_frac(acc, times),
                    events = d, events_control = core$events_control,
                    events_treatment = core$events_treatment,
                    info_frac = t_frac, z_upper = bnd$upper,
                    p_upper = bnd$nominal_p,
                    hr_upper = exp(-bnd$upper / root_info),
                    z_lower = bnd$lower, hr_lower = hr_lower,
                    z_mean = core$z_mean, z_sd = sqrt(diag(sigma)),
                    prob_upper = pr$prob_upper, prob_lower = pr$prob_lower,
                    cum_power = cumsum(pr$prob_upper),
                    alpha_spent = bnd$alpha_spent)
  if (any(!is.na(targets) & lk$trigger != ep)) {
    out$events_trigger <- targets
  }
  if (length(unique(lk$trigger)) > 1L) {
    out <- cbind(out[1L], trigger = lk$trigger, out[-1L], stringsAsFactors = FALSE)
  }
  res <- list(looks = out, power = sum(pr$prob_upper), n = n,
              accrual_duration = acc$duration, study_duration = times[k_max],
              events = d[k_max],
              type1 = c(nonbinding = NA_real_, binding = NA_real_),
              settings = list(endpoint = ep, trigger = looks$trigger,
                              mean = test$mean, variance = test$variance,
                              futility = test$futility, binding = test$binding,
                              spending_time = test$spending_time,
                              alpha = test$alpha),
              scenario = scenario, looks_def = looks, test = test)
  if (type1) {
    t1 <- .type1_core(scenario, acc, looks, targets, ep, bnd)
    res$type1 <- t1$type1
    res$looks$cum_alpha <- t1$cum_nonbinding
  }
  res
}
