# Expected values: inst/validation/python/results/cond_power_stage4.csv (cases M5_*; conditional
# moments from the precision matrix and probabilities of the later looks by recursive conditioning
# with adaptive quadrature in expected_stage4.py, m5_inputs())

m5_inputs <- function() {
  tt <- 1:5
  sd <- sqrt(c(1, 1.1, 0.95, 1.05, 1))
  list(mean = c(0.8, 1.3, 1.7, 2.0, 2.3),
       sigma = (0.9 * sqrt(outer(tt, tt, pmin) / outer(tt, tt, pmax)) + 0.1) * outer(sd, sd),
       lower = c(-0.5, 0.2, 0.8, 1.4, 2.0), upper = c(3.5, 2.9, 2.5, 2.2, 2.0))
}

test_that(".cond_probs matches the reference for a covariance without independent increments", {
  x <- m5_inputs()
  # given the statistic of look 1, the last later look is a rectangle of four coordinates whose
  # covariance is not that of a Markov sequence, computed by the Genz-Bretz algorithm (error
  # tolerance 1e-6 with with_abseps()); given look 2 there are at most three coordinates (Miwa)
  cases <- list(list("M5_l1", 1L, 0.9, 1e-5), list("M5_l2_hist", 2L, c(0.9, 1.4), 1e-6),
                list("M5_l2_cur", 2L, 1.4, 1e-6))
  for (cs in cases) {
    lab <- cs[[1]]
    tol <- cs[[4]]
    r <- with_abseps(.cond_probs(x$mean, x$sigma, cs[[2]], cs[[3]], x$lower, x$upper))
    expect_identical(r$looks$look, seq.int(cs[[2]] + 1L, 5L))
    expect_equal(r$looks$cond_mean, stage4_cp(lab, "cond_mean"), tolerance = 1e-10, label = lab)
    expect_equal(r$looks$cond_sd, stage4_cp(lab, "cond_sd"), tolerance = 1e-10, label = lab)
    expect_lt(max(abs(r$looks$prob_upper - stage4_cp(lab, "prob_upper"))), tol, label = lab)
    expect_lt(max(abs(r$looks$prob_lower - stage4_cp(lab, "prob_lower"))), tol, label = lab)
    expect_lt(abs(r$cond_power - stage4_cp(lab, "cond_power")), tol, label = lab)
    expect_lt(abs(r$cond_power_nonbinding - stage4_cp(lab, "cond_power_nonbinding")), tol,
              label = lab)
    expect_equal(r$looks$cum_upper, cumsum(r$looks$prob_upper))
  }
  expect_false(.is_markov(.cond_moments(x$mean, x$sigma, 1, 2:5, 0.9)$sigma))
})

test_that(".cond_probs uses the earlier statistics and drops the later futility bounds", {
  x <- m5_inputs()
  h <- .cond_probs(x$mean, x$sigma, 2L, c(0.9, 1.4), x$lower, x$upper)
  cu <- .cond_probs(x$mean, x$sigma, 2L, 1.4, x$lower, x$upper)
  # the covariance is not that of a Markov sequence, so that the statistic of look 1 matters
  expect_gt(abs(h$cond_power - cu$cond_power), 1e-4)
  lo <- x$lower
  lo[3:4] <- -Inf
  nb <- .cond_probs(x$mean, x$sigma, 2L, 1.4, lo, x$upper)
  expect_equal(cu$cond_power_nonbinding, nb$cond_power)
  expect_equal(nb$cond_power_nonbinding, nb$cond_power)
  expect_gt(cu$cond_power_nonbinding, cu$cond_power)
})
