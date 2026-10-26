test_that("ms_scenario stores the allocation probabilities", {
  s <- scn_b()
  expect_s3_class(s, "ms_scenario")
  expect_equal(unname(s$prob), c(1 / 3, 2 / 3))
  expect_equal(names(s$prob), c("control", "treatment"))
})

test_that("ms_scenario checks its input", {
  m <- ms_model(death = 0.1)
  expect_error(ms_scenario(m, list(), accrual = ms_accrual(rate = 1)), "ms_model")
  expect_error(ms_scenario(m, m), "ms_accrual")
  expect_error(ms_scenario(m, m, allocation = 0, accrual = ms_accrual(rate = 1)),
               "positive")
})

test_that("ms_scenario stores pfs_switch in the scenario and in both models", {
  s <- scn_s3("PRC")
  expect_identical(s$pfs_switch, "censor")
  expect_identical(c(s$control$pfs_switch, s$treatment$pfs_switch), c("censor", "censor"))
  expect_identical(scn_c()$pfs_switch, "event")
  expect_identical(scn_c()$control$pfs_switch, "event")
  expect_identical(ms_null(s)$treatment$pfs_switch, "censor")
  m <- ms_model(death = 0.1)
  expect_error(ms_scenario(m, m, accrual = ms_accrual(rate = 1), pfs_switch = "x"),
               "should be one of")
})

test_that("ms_scenario stores switch_start and checks it", {
  expect_identical(scn_c()$switch_start, 0)
  expect_identical(scn_cal()$switch_start, 15)
  expect_identical(ms_null(scn_cal())$switch_start, 15)
  m <- models_c()
  acc <- ms_accrual(rate = 20)
  for (bad in list(-1, NA_real_, Inf, c(1, 2), "15")) {
    expect_error(ms_scenario(m$control, m$treatment, accrual = acc, switch_start = bad),
                 "'switch_start' must be a nonnegative number")
  }
  expect_error(ms_scenario(models_s3("P"), models_s3("T"), accrual = acc, switch_start = 10),
               "cannot be combined with switching before progression")
  expect_identical(ms_scenario(models_s3("P"), models_s3("T"), accrual = acc)$switch_start, 0)
})
