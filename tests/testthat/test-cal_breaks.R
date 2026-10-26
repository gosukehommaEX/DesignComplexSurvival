test_that(".cal_breaks keeps the points inside the range without near duplicates", {
  x <- c(5, 2, 2 + 1e-12, 0, 1e-12, 7, 10 - 1e-12, 12, 3, 3)
  expect_equal(.cal_breaks(x, 0, 10), c(2, 3, 5, 7))
  expect_equal(.cal_breaks(c(1, 1 + 2e-8), 0, 10), c(1, 1 + 2e-8))
  expect_equal(.cal_breaks(c(4, 4 + 1e-10), 0, 10), 4)
  expect_equal(.cal_breaks(numeric(0), 0, 10), numeric(0))
  expect_equal(.cal_breaks(c(-1, 11), 0, 10), numeric(0))
})
