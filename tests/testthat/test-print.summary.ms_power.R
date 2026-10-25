test_that("print.summary.ms_power prints all columns", {
  res <- ms_power(scn_a(), n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
                  test = ms_test(), events = 300)
  expect_output(print(summary(res)), "Overview")
  expect_output(print(summary(res)), "cum_alpha")
})

test_that("print.summary.ms_power prints a design with two endpoints", {
  res <- run_design("D3")
  expect_output(print(summary(res)), "fixed sequence \\(stagewise\\)")
  expect_output(print(summary(res)), "Familywise error rate")
})

test_that("print.summary.ms_power shows the transition weights", {
  res <- run_design("D5")
  expect_output(print(summary(res)), "Bonferroni \\(weights 0.8, 0.2; transition weights 1, 1\\)")
  expect_output(print(summary(res)), "z_upper_recycled")
})
