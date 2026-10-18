test_that("print.ms_looks describes the looks", {
  expect_output(print(ms_looks("os", fraction = c(0.5, 1))), "driven by OS events")
  expect_output(print(ms_looks("calendar", time = c(24, 36))), "calendar times")
})
