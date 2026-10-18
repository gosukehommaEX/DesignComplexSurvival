test_that("print.ms_model shows both sets of intensities", {
  m <- models_c()$control
  expect_output(print(m), "Four-state model")
  expect_output(print(m), "After progression")
  expect_output(print(ms_model(death = 0.1)), "Before progression")
  expect_output(expect_invisible(print(m)), "Four-state model")
})
