#' Operating characteristics of a design with two endpoints
#'
#' Determines the analysis times, the joint distribution of the statistics of
#' the two endpoints at the looks at which they are tested, the bounds of each
#' endpoint and the rejection probabilities under the multiple testing
#' procedure, and optionally the familywise error rate under the global null
#' scenario and under \code{null_scenarios}.
#'
#' The bounds of each endpoint are computed by \code{gs_bounds()} with the
#' information fractions of the looks at which the endpoint is tested and the
#' null correlations sqrt(d_j / d_k) of its expected events, at level alpha
#' (fixed sequence) or w_k alpha (Bonferroni). Futility bounds are computed for
#' the first endpoint only. With recycling of the level
#' (\code{test$transition}), the bounds at the level after recycling are
#' computed by \code{.recycled_upper()}, with the binding futility bounds for
#' the first endpoint when they are binding. The rejection probabilities are
#' computed by \code{gs_prob_mtp()} with the mean vector and covariance matrix
#' of \code{.zeta_cov()}.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param n Total sample size.
#' @param looks An \code{ms_looks} object.
#' @param test An \code{ms_test} object with two endpoints.
#' @param events Final numbers of events of the trigger endpoints.
#' @param followup Time from the end of accrual to the final analysis.
#' @param beta Total type II error for beta spending (defaults to
#'   \code{test$beta}).
#' @param type1 Whether to compute the familywise error rate.
#' @param null_scenarios Checked list of null scenarios (see
#'   \code{.fwer_core()}), or \code{NULL}.
#' @return A list (see \code{\link{ms_power}}).
#' @noRd
.design_core2 <- function(scenario, n, looks, test, events = NULL, followup = NULL,
                          beta = NULL, type1 = TRUE, null_scenarios = NULL) {
  acc <- .acc_setup(scenario$accrual, n)
  eps <- test$endpoints
  tl <- .tested_looks(looks, eps)
  lk <- .looks_table(looks)
  lt <- .look_times(scenario, acc, lk, events = events, followup = followup,
                    what = "the design scenario")
  times <- lt$time
  item_ep <- rep(eps, lengths(tl))
  item_lk <- unlist(tl, use.names = FALSE)
  idx <- split(seq_along(item_ep), factor(item_ep, levels = eps))
  cores <- lapply(eps, function(e) .lr_core(scenario, acc, times[tl[[e]]], e, test$mean))
  sigma <- .zeta_cov(scenario, acc, item_ep, times[item_lk], lt$trigger[item_lk],
                     test$variance)$sigma
  mean <- unlist(lapply(cores, `[[`, "z_mean"), use.names = FALSE)
  p <- scenario$prob
  alpha_e <- if (test$mtp == "bonferroni") test$alpha * test$weights else rep(test$alpha, 2L)
  bnd <- vector("list", 2L)
  for (j in 1:2) {
    core <- cores[[j]]
    kj <- length(tl[[j]])
    d <- core$events
    t_frac <- if (test$spending_time == "events") d / d[kj] else core$info / core$info[kj]
    corr0 <- sqrt(outer(d, d, pmin) / outer(d, d, pmax))
    if (test$sf_eff[j] == "user" && length(test$param_eff[[j]]) != kj) {
      stop(sprintf("'param_eff' must have one value per look of %s for sf_eff = \"user\".",
                   toupper(eps[j])), call. = FALSE)
    }
    if (j == 2L) {
      bnd[[j]] <- gs_bounds(t_frac, alpha = alpha_e[j], sf_eff = test$sf_eff[j],
                            param_eff = test$param_eff[[j]], corr = corr0)
      next
    }
    fut <- NULL
    if (test$futility == "fixed") {
      fb <- test$fut_bound
      if (length(fb) != kj - 1L) {
        stop("'fut_bound' must have one value per interim look of ", toupper(eps[1L]),
             " (", kj - 1L, ").", call. = FALSE)
      }
      fut <- if (test$fut_scale == "z") {
        fb
      } else {
        -log(fb) * sqrt(p[[1L]] * p[[2L]] * d[seq_len(kj - 1L)])
      }
    }
    if (test$futility == "beta" && test$sf_fut == "user" && length(test$param_fut) != kj) {
      stop("'param_fut' must have one value per look of ", toupper(eps[1L]),
           " for sf_fut = \"user\".", call. = FALSE)
    }
    bnd[[j]] <- gs_bounds(t_frac, alpha = alpha_e[j], sf_eff = test$sf_eff[j],
                          param_eff = test$param_eff[[j]], futility = test$futility,
                          fut_bound = fut, sf_fut = test$sf_fut,
                          param_fut = test$param_fut,
                          beta = if (is.null(beta)) test$beta else beta,
                          binding = test$binding, mean = core$z_mean, corr = corr0,
                          sigma_alt = sigma[idx[[1L]], idx[[1L]], drop = FALSE])
  }
  up_r <- NULL
  tr <- test$transition
  if (test$mtp == "bonferroni" && !is.null(tr) && any(tr > 0)) {
    lev <- test$alpha * (test$weights + rev(tr * test$weights))
    up_r <- lapply(1:2, function(j) {
      if (!(lev[j] > alpha_e[j])) {
        return(bnd[[j]]$upper)
      }
      d <- cores[[j]]$events
      lw <- if (j == 1L && test$binding && test$futility != "none") bnd[[1L]]$lower
      .recycled_upper(bnd[[j]]$info_frac, lev[j], test$sf_eff[j], test$param_eff[[j]],
                      sqrt(outer(d, d, pmin) / outer(d, d, pmax)), bnd[[j]]$upper, lw, eps[j])
    })
  }
  lower <- if (test$futility == "none") NULL else bnd[[1L]]$lower
  pr <- gs_prob_mtp(list(bnd[[1L]]$upper, bnd[[2L]]$upper), lower = lower, mean = mean,
                    sigma = sigma, analysis = unname(tl), mtp = test$mtp,
                    hierarchical = test$hierarchical, upper_recycled = up_r)
  out <- do.call(rbind, lapply(1:2, function(j) {
    core <- cores[[j]]
    b <- bnd[[j]]
    lks <- tl[[j]]
    root_info <- sqrt(p[[1L]] * p[[2L]] * core$events)
    hr_lower <- exp(-b$lower / root_info)
    hr_lower[is.infinite(b$lower)] <- NA_real_
    pj <- pr$looks[pr$looks$endpoint == j, ]
    df <- data.frame(endpoint = eps[j], look = lks, trigger = lk$trigger[lks],
               time = times[lks], enrolled = n * .acc_frac(acc, times[lks]),
               events = core$events, events_control = core$events_control,
               events_treatment = core$events_treatment, info_frac = b$info_frac,
               z_upper = b$upper, p_upper = b$nominal_p,
               hr_upper = exp(-b$upper / root_info), z_lower = b$lower,
               hr_lower = hr_lower, z_mean = core$z_mean,
               z_sd = sqrt(diag(sigma)[idx[[j]]]), prob_upper = pj$prob_upper,
               prob_lower = pj$prob_lower, cum_power = pj$cum_upper,
               alpha_spent = b$alpha_spent, events_trigger = lt$target[lks],
               stringsAsFactors = FALSE)
    if (!is.null(up_r)) {
      pos <- match("hr_upper", names(df))
      df <- cbind(df[seq_len(pos)], z_upper_recycled = up_r[[j]],
                  p_upper_recycled = stats::pnorm(up_r[[j]], lower.tail = FALSE),
                  df[-seq_len(pos)])
    }
    df
  }))
  rownames(out) <- NULL
  rj <- pr$reject
  power <- c(rj[["first"]], rj[["second"]], rj[["both"]], rj[["any"]])
  names(power) <- c(eps, "both", "any")
  ev_last <- vapply(1:2, function(j) cores[[j]]$events[length(tl[[j]])], numeric(1))
  names(ev_last) <- eps
  res <- list(looks = out, power = power, n = n, accrual_duration = acc$duration,
              study_duration = times[length(times)], events = ev_last,
              final_events = lt$final,
              type1 = c(nonbinding = NA_real_, binding = NA_real_), fwer = NULL,
              mean = mean, sigma = sigma,
              settings = list(endpoint = eps, trigger = looks$trigger, mtp = test$mtp,
                              weights = test$weights, hierarchical = test$hierarchical,
                              mean = test$mean, variance = test$variance,
                              futility = test$futility, binding = test$binding,
                              spending_time = test$spending_time, alpha = test$alpha,
                              transition = test$transition),
              scenario = scenario, looks_def = looks, test = test)
  if (type1) {
    fw <- .fwer_core(scenario, acc, lk, lt$target, tl, test, bnd, null_scenarios,
                     upper_recycled = up_r)
    res$type1 <- fw$type1
    res$fwer <- fw$fwer
    res$looks$cum_alpha <- fw$cum_alpha
  }
  res
}
