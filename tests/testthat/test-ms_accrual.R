test_that("accrual by rates reaches the sample size at the right time", {
  ctl <- ms_model(death = 0.05)
  scn <- ms_scenario(ctl, ctl, accrual = ms_accrual(time = c(0, 3), rate = c(10, 30)))
  # 30 patients by month 3, then 30 per month: 300 patients at month 12
  e <- ms_events(scn, n = 300, calendar_time = c(2, 6, 12, 20))
  expect_equal(e$enrolled, c(20, 120, 300, 300), tolerance = 1e-12)
})

test_that("accrual by proportions is uniform within intervals", {
  ctl <- ms_model(death = 0.05)
  acc <- ms_accrual(time = c(0, 6), prop = c(0.3, 0.7), duration = 18)
  scn <- ms_scenario(ctl, ctl, accrual = acc)
  e <- ms_events(scn, n = 500, calendar_time = c(3, 6, 12, 18, 24))
  expect_equal(e$enrolled, c(75, 150, 325, 500, 500), tolerance = 1e-12)
})

test_that("ms_accrual rejects invalid input", {
  expect_error(ms_accrual(rate = 0), "last accrual rate")
  expect_error(ms_accrual(prop = 1), "prop")
  expect_error(ms_accrual(time = c(0, 6), prop = c(0.5, 0.4), duration = 12), "sum to 1")
  expect_error(ms_accrual(time = c(0, 12), prop = c(0.5, 0.5), duration = 12), "smaller")
  expect_error(ms_accrual(rate = 10, prop = 1, duration = 5), "either")
})
