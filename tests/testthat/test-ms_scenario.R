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
