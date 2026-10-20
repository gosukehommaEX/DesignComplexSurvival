test_that(".param_list gives one parameter per endpoint", {
  expect_identical(.param_list(NULL, c("OF", "OF"), "p"), list(NULL, NULL))
  expect_identical(.param_list(-2, c("HSD", "HSD"), "p"), list(-2, -2))
  expect_identical(.param_list(c(-4, -2), c("HSD", "HSD"), "p"), list(-4, -2))
  expect_identical(.param_list(list(c(0.2, 1), NULL), c("user", "OF"), "p"),
                   list(c(0.2, 1), NULL))
  expect_error(.param_list(c(0.2, 1), c("user", "OF"), "p"), "'p'")
  expect_error(.param_list(list(1, 2, 3), c("HSD", "HSD"), "p"), "'p'")
})
