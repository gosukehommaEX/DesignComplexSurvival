# .pmvn() uses the Miwa algorithm for at most three coordinates and for
# covariance matrices of Markov sequences, and the Genz-Bretz algorithm
# (.pmvn_genz()) for four or more coordinates otherwise. Inputs and expected
# values: helper-pmvn.R (design V2_K5_T2_alternative).

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
