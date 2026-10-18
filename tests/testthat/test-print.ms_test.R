test_that("print.ms_test describes the test", {
  expect_output(print(ms_test()), "One-sided log-rank test of OS")
  expect_output(print(ms_test(futility = "beta", binding = TRUE)), "binding")
  expect_output(print(ms_test(futility = "fixed", fut_bound = 1, fut_scale = "hr")),
                "hazard ratio")
})
