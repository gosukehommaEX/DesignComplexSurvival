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
