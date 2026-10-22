# .pmvn() uses the Miwa algorithm for at most three coordinates and for
# covariance matrices of Markov sequences, and the Genz-Bretz algorithm
# (.pmvn_genz()) for four or more coordinates otherwise and for nine or more
# coordinates with two finite limits for some coordinate. Inputs and expected
# values: helper-pmvn.R (design V2_K5_T2_alternative, and the rectangles of
# nine and ten coordinates).

test_that(".pmvn uses the Genz-Bretz algorithm for four coordinates that are not Markov", {
  x <- pmvn_v2k5
  i <- 1:4
  lo <- c(x$lower[1:3], x$upper[4])
  hi <- c(x$upper[1:3], Inf)
  expect_identical(.pmvn(lo, hi, x$mean[i], x$sigma[i, i]),
                   .pmvn_genz(lo, hi, x$mean[i], x$sigma[i, i]))
  p <- with_abseps(.pmvn(lo, hi, x$mean[i], x$sigma[i, i]), 1e-6)
  expect_lt(abs(p - x$gs_prob_upper[4]), 2e-6)
})

test_that(".pmvn keeps the Miwa algorithm for three coordinates and for Markov matrices", {
  x <- pmvn_v2k5
  miwa <- function(lo, hi, m, s) {
    max(as.numeric(mvtnorm::pmvnorm(lower = lo, upper = hi, mean = m, sigma = s,
                                    algorithm = mvtnorm::Miwa(steps = 128))), 0)
  }
  # three coordinates with one-sided limits (no replacement of infinite
  # limits in .pmvn())
  lo3 <- c(-Inf, -Inf, x$upper[3])
  hi3 <- c(x$upper[1:2], Inf)
  expect_identical(.pmvn(lo3, hi3, x$mean[1:3], x$sigma[1:3, 1:3]),
                   miwa(lo3, hi3, x$mean[1:3], x$sigma[1:3, 1:3]))
  # four coordinates with the correlation of independent increments
  tt <- x$info_frac[1:4]
  corr <- sqrt(outer(tt, tt, pmin) / outer(tt, tt, pmax))
  lo4 <- c(-Inf, -Inf, -Inf, x$upper[4])
  hi4 <- c(x$upper[1:3], Inf)
  expect_identical(.pmvn(lo4, hi4, x$mean[1:4], corr), miwa(lo4, hi4, x$mean[1:4], corr))
})

test_that(".pmvn uses the Genz-Bretz algorithm for nine coordinates with two finite limits", {
  # the Miwa algorithm of mvtnorm leaves out the terms of nine or more
  # coordinates of the inclusion-exclusion and gave 0.08981538 (version
  # 0.2.1) against the exact value 0.08906131
  x <- pmvn_k9
  expect_identical(.pmvn(x$lower, x$upper, x$mean, x$corr),
                   .pmvn_genz(x$lower, x$upper, x$mean, x$corr))
  p <- with_abseps(.pmvn(x$lower, x$upper, x$mean, x$corr), 1e-6)
  expect_lt(abs(p - x$exact), 2e-6)
})

test_that(".pmvn keeps the Miwa algorithm for ten coordinates with one-sided limits", {
  x <- pmvn_k10
  p <- .pmvn(x$lower, x$upper, x$mean, x$corr)
  expect_identical(p, max(as.numeric(mvtnorm::pmvnorm(lower = x$lower, upper = x$upper,
                                                      mean = x$mean, sigma = x$corr,
                                                      algorithm = mvtnorm::Miwa(steps = 128))),
                          0))
  expect_lt(abs(p - x$exact), 1e-6)
})
