test_that(".acc_cum agrees with .acc_frac", {
  acc <- .acc_setup(ms_accrual(time = c(0, 3), rate = c(10, 30)), 300)
  x <- c(-1, 0, 1.5, 3, 7.2, acc$duration, acc$duration + 5)
  expect_equal(.acc_cum(acc, x), .acc_frac(acc, x), tolerance = 1e-14)
  expect_equal(.acc_cum(acc, c(-1, acc$duration + 1)), c(0, 1), tolerance = 1e-14)
  expect_identical(.acc_cum(acc, numeric(0)), numeric(0))
})
