test_that("ms_test stores the settings with defaults", {
  te <- ms_test()
  expect_s3_class(te, "ms_test")
  expect_equal(te$alpha, 0.025)
  expect_identical(te$endpoints, "os")
  expect_identical(te$futility, "none")
  expect_identical(te$mean, "direct")
  expect_identical(te$variance, "null")
  expect_identical(te$spending_time, "events")
  expect_false(te$binding)
  tb <- ms_test(endpoints = "pfs", futility = "beta", sf_fut = "HSD", param_fut = -2,
                binding = TRUE, mean = "schoenfeld", variance = "alternative")
  expect_identical(tb$sf_fut, "HSD")
  expect_true(tb$binding)
})

test_that("ms_test checks its input", {
  expect_error(ms_test(alpha = 0.6 + 0.5), "alpha")
  expect_error(ms_test(endpoints = c("pfs", "os")), "endpoints")
  expect_error(ms_test(sf_eff = "OBF"), "sf_eff")
  expect_error(ms_test(futility = "fixed"), "fut_bound")
  expect_error(ms_test(futility = "fixed", fut_bound = -1, fut_scale = "hr"), "fut_bound")
  expect_error(ms_test(binding = NA), "binding")
})

test_that("ms_test accepts the studentized variance", {
  expect_identical(ms_test(variance = "studentized")$variance, "studentized")
})

test_that("ms_test defines tests of two endpoints", {
  te <- ms_test(endpoints = c("pfs", "os"), sf_eff = c("OF", "Pocock"), mtp = "fixed_sequence")
  expect_identical(te$endpoints, c("pfs", "os"))
  expect_identical(te$sf_eff, c("OF", "Pocock"))
  expect_identical(te$param_eff, list(NULL, NULL))
  expect_identical(te$mtp, "fixed_sequence")
  expect_identical(te$hierarchical, "overall")
  expect_null(te$weights)
  tb <- ms_test(endpoints = c("os", "pfs"), mtp = "bonferroni", sf_eff = "HSD",
                param_eff = c(-4, -2))
  expect_equal(tb$weights, c(0.5, 0.5))
  expect_identical(tb$sf_eff, c("HSD", "HSD"))
  expect_identical(tb$param_eff, list(-4, -2))
  tw <- ms_test(endpoints = c("pfs", "os"), mtp = "bonferroni", weights = c(0.8, 0.2),
                sf_eff = c("user", "OF"), param_eff = list(c(0.3, 1), NULL))
  expect_equal(tw$weights, c(0.8, 0.2))
  expect_identical(tw$param_eff, list(c(0.3, 1), NULL))
  expect_identical(ms_test()$mtp, "none")
  expect_identical(ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence",
                           hierarchical = "stagewise")$hierarchical, "stagewise")
})

test_that("ms_test checks the input for two endpoints", {
  expect_error(ms_test(endpoints = c("pfs", "pfs"), mtp = "bonferroni"), "endpoints")
  expect_error(ms_test(endpoints = c("pfs", "os", "os")), "endpoints")
  expect_error(ms_test(mtp = "bonferroni"), "one endpoint")
  expect_error(ms_test(endpoints = c("pfs", "os"), mtp = "bonferroni", weights = c(0.8, 0.3)),
               "weights")
  expect_error(ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence",
                       weights = c(0.5, 0.5)), "weights")
  expect_error(ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence",
                       sf_eff = c("OF", "OF", "OF")), "sf_eff")
  expect_error(ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence", sf_eff = "user",
                       param_eff = c(0.5, 1)), "param_eff")
  expect_error(ms_test(sf_eff = c("OF", "Pocock")), "sf_eff")
})
