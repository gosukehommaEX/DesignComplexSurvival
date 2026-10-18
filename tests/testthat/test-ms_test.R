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
