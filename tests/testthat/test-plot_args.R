test_that(".plot_args replaces and adds arguments", {
  a <- .plot_args(list(x = 1:2, y = 3:4, xlab = "a", log = ""), list(xlab = "b", main = "m"))
  expect_identical(a, list(x = 1:2, y = 3:4, xlab = "b", log = "", main = "m"))
  expect_identical(.plot_args(list(x = 1), list()), list(x = 1))
})
