test_that("print.ms_test describes the test", {
  expect_output(print(ms_test()), "One-sided log-rank test of OS")
  expect_output(print(ms_test(futility = "beta", binding = TRUE)), "binding")
  expect_output(print(ms_test(futility = "fixed", fut_bound = 1, fut_scale = "hr")),
                "hazard ratio")
})

test_that("print.ms_test describes tests of two endpoints", {
  expect_output(print(ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence")),
                "Fixed sequence: PFS, then OS \\(overall\\)")
  expect_output(print(ms_test(endpoints = c("pfs", "os"), mtp = "bonferroni",
                              weights = c(0.8, 0.2), futility = "beta")),
                "PFS at level 0.02, OS at level 0.005")
  expect_output(print(ms_test(endpoints = c("pfs", "os"), mtp = "bonferroni",
                              futility = "beta")), "Futility \\(PFS\\)")
})
