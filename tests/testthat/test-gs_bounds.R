# Expected values: inst/validation/python/results/gs_bounds.csv (exact
# recursive integration with adaptive quadrature and root finding)

tf3 <- c(0.5, 0.75, 1)

test_that("efficacy bounds match the reference values", {
  a <- gs_bounds(tf3)
  expect_lt(max(abs(a$upper - c(2.962588042727569, 2.3590177072405347,
                                2.0140836682405148))), 1e-5)
  expect_equal(a$lower[1:2], c(-Inf, -Inf))
  expect_equal(a$lower[3], a$upper[3])
  expect_equal(a$nominal_p, pnorm(a$upper, lower.tail = FALSE))
  b <- gs_bounds(c(1 / 3, 2 / 3, 1), sf_eff = "Pocock")
  expect_lt(max(abs(b$upper - c(2.279428238917005, 2.294911134901885,
                                2.29593957933839))), 1e-5)
})

test_that("beta-spending futility bounds match the reference values", {
  mu <- c(0.8, 1.9, 2.8)
  c_nb <- gs_bounds(tf3, futility = "beta", sf_fut = "HSD", param_fut = -2, beta = 0.1,
                    mean = mu)
  expect_lt(max(abs(c_nb$upper - c(2.962588042727569, 2.3590177072405347,
                                   2.0140836682405148))), 1e-5)
  expect_lt(max(abs(c_nb$lower - c(-1.128537662917749, 0.19207987721354994,
                                   2.0140836682405148))), 1e-5)
  d_b <- gs_bounds(tf3, futility = "beta", sf_fut = "HSD", param_fut = -2, beta = 0.1,
                   mean = mu, binding = TRUE)
  expect_lt(max(abs(d_b$upper - c(2.962588042727569, 2.3590177011481708,
                                  2.013928655942471))), 1e-5)
  expect_lt(max(abs(d_b$lower[1:2] - c(-1.128537662917749, 0.19207987721354994))), 1e-5)
  e_f <- gs_bounds(tf3, futility = "fixed", fut_bound = c(0, 0.5), binding = TRUE)
  expect_lt(max(abs(e_f$upper - c(2.962588042727569, 2.358961190543406,
                                  2.0080612734559042))), 1e-5)
  sd <- c(1.05, 1.03, 1.02)
  rc <- matrix(c(1, 0.80, 0.68, 0.80, 1, 0.85, 0.68, 0.85, 1), 3)
  f_a <- gs_bounds(tf3, futility = "beta", sf_fut = "OF", beta = 0.15,
                   mean = c(1.2, 2.0, 2.7), sigma_alt = outer(sd, sd) * rc)
  expect_lt(max(abs(f_a$lower - c(-0.6170313923930154, 0.5750385493546711,
                                  2.0140836682405148))), 1e-5)
})

test_that("gs_bounds agrees with gsDesign::gsDesign()", {
  skip_if_not_installed("gsDesign")
  x1 <- gsDesign::gsDesign(k = 3, test.type = 1, alpha = 0.025,
                           sfu = gsDesign::sfLDOF, timing = tf3)
  expect_lt(max(abs(gs_bounds(tf3)$upper - x1$upper$bound)), 1e-4)
  for (tt in c(4, 3)) {
    x <- gsDesign::gsDesign(k = 3, test.type = tt, alpha = 0.025, beta = 0.1,
                            sfu = gsDesign::sfLDOF, sfl = gsDesign::sfHSD, sflpar = -2,
                            timing = tf3)
    b <- gs_bounds(x$timing, futility = "beta", sf_fut = "HSD", param_fut = -2,
                   beta = 0.1, mean = x$delta * sqrt(x$n.I), binding = (tt == 3))
    expect_lt(max(abs(b$upper - x$upper$bound)), 1e-4)
    expect_lt(max(abs(b$lower - x$lower$bound)), 1e-4)
  }
})

test_that("gs_bounds agrees with rpact::getDesignGroupSequential()", {
  skip_if_not_installed("rpact")
  d <- rpact::getDesignGroupSequential(kMax = 3, alpha = 0.025, sided = 1,
                                       typeOfDesign = "asOF", informationRates = tf3)
  expect_lt(max(abs(gs_bounds(tf3)$upper - d$criticalValues)), 1e-4)
  dp <- rpact::getDesignGroupSequential(kMax = 3, alpha = 0.025, sided = 1,
                                        typeOfDesign = "asP", informationRates = tf3)
  expect_lt(max(abs(gs_bounds(tf3, sf_eff = "Pocock")$upper - dp$criticalValues)), 1e-4)
})

test_that("gs_bounds agrees with gsDesign2 under a non-constant effect", {
  skip_if_not_installed("gsDesign2")
  skip_if_not_installed("gsDesign")
  info <- 100 * tf3
  mu <- c(0.8, 1.9, 2.8)
  g2 <- as.data.frame(gsDesign2::gs_power_npe(
    theta = mu / sqrt(info), info = info,
    upper = gsDesign2::gs_spending_bound,
    upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL),
    lower = gsDesign2::gs_spending_bound,
    lpar = list(sf = gsDesign::sfHSD, total_spend = 0.1, param = -2, timing = NULL),
    binding = FALSE))
  up2 <- g2$z[g2$bound == "upper"][order(g2$analysis[g2$bound == "upper"])]
  lo2 <- g2$z[g2$bound == "lower"][order(g2$analysis[g2$bound == "lower"])]
  b <- gs_bounds(tf3, futility = "beta", sf_fut = "HSD", param_fut = -2, beta = 0.1,
                 mean = mu)
  expect_lt(max(abs(b$upper - up2)), 1e-4)
  expect_lt(max(abs(b$lower[1:2] - lo2[1:2])), 1e-4)
})

test_that("bounds and probabilities give no warning with mixed infinite limits", {
  # boxes with two finite limits for some looks and one for others made
  # mvtnorm (Miwa) warn "Approximating +/-Inf by +/-1000" before .pmvn()
  # replaced the infinite limits itself
  mu <- c(0.8, 1.9, 2.8)
  expect_warning(b <- gs_bounds(tf3, futility = "beta", sf_fut = "HSD",
                                param_fut = -2, mean = mu), NA)
  corr <- sqrt(outer(tf3, tf3, pmin) / outer(tf3, tf3, pmax))
  expect_warning(gs_prob(b$lower, b$upper, mu, corr), NA)
})

test_that("gs_bounds handles zero spending and checks its input", {
  u <- gs_bounds(tf3, sf_eff = "user", param_eff = c(0, 0.5, 1))
  expect_equal(u$upper[1], Inf)
  expect_error(gs_bounds(c(0.5, 0.9)), "end at 1")
  expect_error(gs_bounds(tf3, futility = "beta"), "mean")
  expect_error(gs_bounds(tf3, futility = "fixed", fut_bound = 0), "fut_bound")
  expect_error(gs_bounds(tf3, corr = diag(2)), "corr")
})
