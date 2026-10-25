#' Familywise error rate of a design with two endpoints
#'
#' Applies the bounds and the multiple testing procedure of the design to the
#' statistics under the global null scenario (\code{ms_null()}) and under
#' each scenario of \code{null_scenarios}. Event-driven looks take place when
#' the same numbers of events of their trigger endpoints are expected under
#' the null scenario; calendar-driven looks keep their calendar times.
#'
#' Under the global null scenario the statistics have mean 0, the statistics
#' of one endpoint have the correlations sqrt(v_j / v_k) and the
#' correlations between the endpoints are those of \code{.zeta_cov()} (the
#' covariance of Danzer et al. 2025, Theorem 1). Under a scenario of
#' \code{null_scenarios}, the mean vector and covariance matrix are computed
#' with the approximations of the test, as for the design scenario, and the
#' familywise error rate is the probability of rejecting at least one of the
#' hypotheses of the endpoints named by the scenario.
#'
#' @param scenario An \code{ms_scenario} object (design scenario).
#' @param acc Output of \code{.acc_setup()}.
#' @param lk Output of \code{.looks_table()}.
#' @param targets Event targets of the looks (\code{NA} for calendar-driven
#'   looks).
#' @param tl Output of \code{.tested_looks()}.
#' @param test An \code{ms_test} object with two endpoints.
#' @param bnd List of the outputs of \code{gs_bounds()} for the two endpoints.
#' @param null_scenarios \code{NULL}, or a list of \code{ms_scenario} objects
#'   named by the endpoint (\code{"os"} or \code{"pfs"}) whose hypothesis is
#'   true in the scenario.
#' @param upper_recycled Bounds of the two endpoints at the level after
#'   recycling (see \code{gs_prob_mtp()}), or \code{NULL} for none.
#' @return A list with \code{type1} (familywise error rate under the global
#'   null scenario, non-binding and binding), \code{fwer} (a data frame with
#'   one row per null scenario) and \code{cum_alpha} (cumulative non-binding
#'   probability of rejecting the hypothesis of each endpoint under the global
#'   null scenario, in the order of the tested looks).
#' @noRd
.fwer_core <- function(scenario, acc, lk, targets, tl, test, bnd, null_scenarios = NULL,
                       upper_recycled = NULL) {
  eps <- test$endpoints
  upper <- list(bnd[[1L]]$upper, bnd[[2L]]$upper)
  lower <- if (test$futility == "none") NULL else bnd[[1L]]$lower
  an <- unname(tl)
  item_ep <- rep(eps, lengths(tl))
  item_lk <- unlist(tl, use.names = FALSE)
  trig <- ifelse(lk$trigger == "calendar", NA_character_, lk$trigger)
  probs <- function(mean0, sig) {
    nb <- gs_prob_mtp(upper, NULL, mean0, sig, an, test$mtp, test$hierarchical,
                      upper_recycled = upper_recycled)
    bd <- if (is.null(lower)) nb else gs_prob_mtp(upper, lower, mean0, sig, an,
                                                   test$mtp, test$hierarchical,
                                                   upper_recycled = upper_recycled)
    list(nb = nb, bd = bd)
  }
  s0 <- ms_null(scenario)
  t0 <- .look_times(s0, acc, lk, targets = targets, what = "the global null scenario")$time
  pr0 <- probs(0, .zeta_cov(s0, acc, item_ep, t0[item_lk], NULL, "null")$sigma)
  type1 <- c(nonbinding = pr0$nb$reject[["any"]], binding = pr0$bd$reject[["any"]])
  fwer <- data.frame(scenario = "global", null_endpoints = paste(toupper(eps), collapse = ", "),
                     nonbinding = type1[["nonbinding"]], binding = type1[["binding"]],
                     stringsAsFactors = FALSE)
  for (q in seq_along(null_scenarios)) {
    sq <- null_scenarios[[q]]
    lab <- sprintf("null_scenarios[[%d]]", q)
    tq <- .look_times(sq, acc, lk, targets = targets, what = lab)$time
    mq <- unlist(lapply(eps, function(e) {
      .lr_core(sq, acc, tq[tl[[e]]], e, test$mean)$z_mean
    }), use.names = FALSE)
    sig <- .zeta_cov(sq, acc, item_ep, tq[item_lk], trig[item_lk], test$variance)$sigma
    prq <- probs(mq, sig)
    key <- if (names(null_scenarios)[q] == eps[1L]) "first" else "second"
    fwer <- rbind(fwer, data.frame(scenario = lab,
                                   null_endpoints = toupper(names(null_scenarios)[q]),
                                   nonbinding = prq$nb$reject[[key]],
                                   binding = prq$bd$reject[[key]], stringsAsFactors = FALSE))
  }
  list(type1 = type1, fwer = fwer, cum_alpha = pr0$nb$looks$cum_upper)
}
