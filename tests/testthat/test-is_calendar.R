test_that(".is_calendar needs a positive switch_start and a group that switches", {
  expect_true(.is_calendar(scn_cal()))
  expect_false(.is_calendar(scn_cal(tc = 0)))
  expect_false(.is_calendar(scn_c()))
  m <- models_c()
  s <- ms_scenario(m$treatment, m$treatment, accrual = ms_accrual(rate = 20), switch_start = 10)
  expect_false(.is_calendar(s))
  old <- scn_cal()
  old$switch_start <- NULL
  expect_false(.is_calendar(old))
})
