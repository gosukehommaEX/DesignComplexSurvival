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

test_that("print.ms_test shows the recycling of the level", {
  te <- ms_test(endpoints = c("pfs", "os"), mtp = "bonferroni", weights = c(0.8, 0.2),
                transition = c(1, 0.5))
  expect_output(print(te), "Recycling: fraction 1 from PFS to OS and 0.5 from OS to PFS")
  expect_output(print(te), paste("OS at level 0.025 after PFS is rejected, PFS at level 0.0225",
                                 "after OS is rejected"))
  out <- capture.output(print(ms_test(endpoints = c("pfs", "os"), mtp = "bonferroni")))
  expect_false(any(grepl("Recycling", out)))
})
