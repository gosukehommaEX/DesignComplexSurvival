test_that("ms_null copies the control model and gives a zero score mean", {
  s0 <- ms_null(scn_b())
  expect_identical(s0$treatment, s0$control)
  expect_equal(s0$allocation, 2)
  mo <- lr_moments(s0, n = 300, calendar_time = c(15, 25))
  expect_lt(max(abs(mo$score_mean)), 1e-10)
  # under the null hypothesis the information is p0 p1 times the events
  expect_lt(rel_diff(mo$info, mo$events * 2 / 9), 1e-10)
  expect_error(ms_null(list()), "ms_scenario")
})
