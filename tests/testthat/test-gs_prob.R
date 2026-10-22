# Expected values: inst/validation/python/results/gs_prob.csv

test_that("gs_prob matches the reference stopping probabilities", {
  lower <- c(-1.128537662917749, 0.19207987721354994, 2.0140836682405148)
  upper <- c(2.962588042727569, 2.3590177072405347, 2.0140836682405148)
  tt <- c(0.5, 0.75, 1)
  corr <- sqrt(outer(tt, tt, pmin) / outer(tt, tt, pmax))
  h1 <- gs_prob(lower, upper, mean = c(0.8, 1.9, 2.8), sigma = corr)
  expect_lt(max(abs(h1$prob_upper - c(0.015286439374457172, 0.30789497811439753,
                                      0.4569427960945711))), 1e-6)
  expect_lt(max(abs(h1$prob_lower - c(0.026894142136999408, 0.027600434470659366,
                                      0.16538120980891566))), 1e-6)
  expect_equal(h1$cum_upper[3] + h1$cum_lower[3], 1, tolerance = 1e-6)
  h0 <- gs_prob(lower, upper, mean = 0, sigma = corr)
  expect_lt(max(abs(h0$prob_upper - c(0.0015253227579888895, 0.00812400205057425,
                                      0.015343484178349576))), 1e-6)
  expect_lt(max(abs(h0$prob_lower - c(0.12954646011832432, 0.44776793447248414,
                                      0.39769279642227906))), 1e-6)
})

test_that("gs_prob is exact with five looks and a covariance that is not Markov", {
  # design V2_K5_T2_alternative (helper-pmvn.R); with the Miwa algorithm
  # (version 0.2.0) the probabilities of looks 4 and 5 were off by up to
  # 4.0e-3
  x <- pmvn_v2k5
  p <- with_abseps(gs_prob(x$lower, x$upper, x$mean, x$sigma), 1e-6)
  expect_lt(max(abs(p$prob_upper - x$gs_prob_upper)), 2e-6)
  expect_lt(max(abs(p$prob_lower - x$gs_prob_lower)), 2e-6)
})

test_that("gs_prob without futility spends exactly alpha under the null", {
  b <- gs_bounds(c(0.5, 0.75, 1))
  tt <- c(0.5, 0.75, 1)
  corr <- sqrt(outer(tt, tt, pmin) / outer(tt, tt, pmax))
  p0 <- gs_prob(rep(-Inf, 3), b$upper, 0, corr)
  expect_lt(max(abs(p0$cum_upper - b$alpha_spent)), 1e-8)
  expect_error(gs_prob(0, c(1, 2), 0, diag(2)), "same length")
  expect_error(gs_prob(c(0, 0), c(1, 2), 0, diag(3)), "sigma")
})
