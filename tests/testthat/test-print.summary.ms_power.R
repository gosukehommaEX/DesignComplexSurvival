test_that("print.summary.ms_power prints all columns", {
  res <- ms_power(scn_a(), n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
                  test = ms_test(), events = 300)
  expect_output(print(summary(res)), "Overview")
  expect_output(print(summary(res)), "cum_alpha")
})
