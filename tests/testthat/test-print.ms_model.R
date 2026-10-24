test_that("print.ms_model shows both sets of intensities", {
  m <- models_c()$control
  expect_output(print(m), "Four-state model")
  expect_output(print(m), "After progression")
  expect_output(print(ms_model(death = 0.1)), "Before progression")
  expect_output(expect_invisible(print(m)), "Four-state model")
})

test_that("print.ms_model shows the options of version 0.3.0", {
  expect_output(print(models_s3("P")), "Five-state model")
  expect_output(print(models_s3("P")), "Switch before progression")
  expect_output(print(models_s3("W")), "time since the switch")
  expect_output(print(models_s3("A06")), "accelerated failure time effect with factor 0.6")
  expect_output(print(models_s3("T0", dropout_pfs = 0.01)), "dropout_pfs")
  expect_output(print(models_s3("Z")), "After progression")
  out <- capture.output(print(models_c()$control))
  expect_false(any(grepl("dropout_pfs|Switch before|accelerated", out)))
})
