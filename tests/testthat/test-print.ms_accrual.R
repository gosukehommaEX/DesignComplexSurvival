test_that("print.ms_accrual shows the accrual pattern", {
  expect_output(print(ms_accrual(time = c(0, 3), rate = c(10, 30))), "by rates")
  expect_output(print(ms_accrual(prop = 1, duration = 12)), "by proportions")
})
