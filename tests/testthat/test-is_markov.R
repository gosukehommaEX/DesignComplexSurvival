# Expected values: inst/validation/python/results/pmvn.csv (case V2K5 of
# helper-pmvn.R: the largest element of the inverse correlation matrix
# outside the tridiagonal band is 0.0039 of the largest element for the five
# looks and 0.0010 for the first three)

test_that(".is_markov recognizes the covariance of independent increments", {
  tt <- c(0.2, 0.4, 0.6, 0.8, 1)
  corr <- sqrt(outer(tt, tt, pmin) / outer(tt, tt, pmax))
  expect_true(.is_markov(corr))
  sd <- c(1.2, 0.9, 1.1, 1, 0.95)
  expect_true(.is_markov(outer(sd, sd) * corr))
  expect_true(.is_markov(matrix(c(1, 0.3, 0.3, 1), 2)))
  expect_true(.is_markov(matrix(2)))
})

test_that(".is_markov rejects the covariance of the alternative variance", {
  x <- pmvn_v2k5
  expect_false(.is_markov(x$sigma))
  expect_false(.is_markov(x$sigma[1:3, 1:3]))
  expect_true(.is_markov(x$sigma[1:2, 1:2]))
  expect_false(.is_markov(matrix(1, 3, 3)))
})
