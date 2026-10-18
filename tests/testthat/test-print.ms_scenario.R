test_that("print.ms_scenario shows the groups and the accrual", {
  out <- capture.output(print(scn_b()))
  expect_true(any(grepl("allocation 2:1", out)))
  expect_true(any(grepl("Treatment group", out)))
  expect_true(any(grepl("accrual", out)))
})
