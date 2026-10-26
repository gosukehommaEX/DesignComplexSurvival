test_that("print.ms_scenario shows the groups and the accrual", {
  out <- capture.output(print(scn_b()))
  expect_true(any(grepl("allocation 2:1", out)))
  expect_true(any(grepl("Treatment group", out)))
  expect_true(any(grepl("accrual", out)))
})

test_that("print.ms_scenario shows the definition of PFS at a switch before progression", {
  expect_true(any(grepl("PFS at a switch before progression: censored",
                        capture.output(print(scn_s3("PRC"))))))
  expect_true(any(grepl("PFS at a switch before progression: followed up",
                        capture.output(print(scn_s3("PR"))))))
  expect_false(any(grepl("PFS at a switch", capture.output(print(scn_c())))))
})

test_that("print.ms_scenario shows the calendar time from which switching is possible", {
  expect_true(any(grepl("Switching possible from calendar time 15 ",
                        capture.output(print(scn_cal())))))
  expect_false(any(grepl("Switching possible", capture.output(print(scn_c())))))
})
