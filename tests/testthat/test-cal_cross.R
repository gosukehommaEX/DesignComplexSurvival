test_that(".cal_cross gives the meeting points inside the range of the availability", {
  acc <- .acc_setup(scn_cal()$accrual, 500)
  ln <- .cal_lines(models_c("forward")$control)
  # t = 36, switching from month 15, accrual breaks 0, 6 and 18: o = 3 and 4 meet o = s at
  # s = 3, 4 and o = s - 21 at s = 24, 25; o = 0 meets o = s - 21 at s = 21; the edges o = 15
  # and 9 meet o = s at s = 15, 9, and o = 9 meets o = s - 21 at s = 30 (o = -3 at s = 18 is
  # outside the range)
  expect_equal(sort(unique(.cal_cross(ln, 36, 15, acc))), c(3, 4, 9, 15, 21, 24, 25, 30))
  expect_equal(sort(unique(.cal_cross(ln[0L, , drop = FALSE], 36, 15, acc))), c(9, 15, 21, 30))
})
