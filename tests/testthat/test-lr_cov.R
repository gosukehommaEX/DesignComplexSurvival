# Expected values: inst/validation/python/results/alt_cov.csv (per-patient
# U, I and the covariance CV(t, t') of Luo et al. 2019, Section 3.1, by
# nested adaptive quadrature; the influence-function derivation used in the
# R code was checked against these formulas in Python to about 1e-14)

test_that("lr_cov with the null variance uses the information ratios", {
  lc <- lr_cov(scn_b(), n = 300, calendar_time = c(15, 25))
  expect_equal(diag(lc$sigma), c(1, 1))
  expect_equal(lc$sigma[1, 2], sqrt(lc$info[1] / lc$info[2]))
  expect_equal(lc$info, c(18.817550458247293, 29.418218596536327), tolerance = 1e-10)
  expect_equal(lc$mean, c(0.3652132322828265, 1.174471710657412), tolerance = 1e-10)
})

test_that("lr_cov with the alternative variance matches Luo et al. (2019)", {
  lc <- lr_cov(scn_b(), n = 300, calendar_time = c(15, 25), variance = "alternative")
  i15 <- 0.06272516819415765
  i25 <- 0.09806072865512112
  expect_equal(lc$info / 300, c(i15, i25), tolerance = 1e-10)
  expect_equal(diag(lc$sigma), c(0.06317798711200127 / i15, 0.09823341210532265 / i25),
               tolerance = 1e-8)
  expect_equal(lc$sigma[1, 2], 0.06113467590571379 / sqrt(i15 * i25), tolerance = 1e-8)
  expect_equal(lc$sigma[1, 2], lc$sigma[2, 1])
})

test_that("lr_cov with the alternative variance reduces to the null under the null", {
  lc <- lr_cov(ms_null(scn_b()), n = 300, calendar_time = c(15, 25),
               variance = "alternative")
  ln <- lr_cov(ms_null(scn_b()), n = 300, calendar_time = c(15, 25))
  expect_lt(max(abs(lc$sigma - ln$sigma)), 1e-8)
  expect_error(lr_cov(scn_b(), n = 300, calendar_time = c(25, 15)), "increasing")
})
