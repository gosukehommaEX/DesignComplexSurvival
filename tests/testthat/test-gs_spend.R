# Expected values: inst/validation/python/results/gs_spend.csv

test_that("gs_spend matches the reference values", {
  tt <- c(0.25, 0.5, 0.75, 1)
  expect_equal(gs_spend(tt, 0.025, "OF"),
               c(7.36680843586938e-06, 0.0015253227579888986, 0.009649324953511989,
                 0.024999999999999894), tolerance = 1e-12)
  expect_equal(gs_spend(tt, 0.025, "Pocock"),
               c(0.00893435048771971, 0.015502862673956938, 0.020699723481071745, 0.025),
               tolerance = 1e-12)
  expect_equal(gs_spend(tt, 0.025, "HSD", -4),
               c(0.0008014650820021247, 0.0029800730505529396, 0.008902143502800693, 0.025),
               tolerance = 1e-12)
  expect_equal(gs_spend(tt, 0.025, "HSD", 0), 0.025 * tt, tolerance = 1e-14)
  expect_equal(gs_spend(tt, 0.025, "power", 3), 0.025 * tt^3, tolerance = 1e-14)
  expect_equal(gs_spend(tt, 0.025, "user", c(0.1, 0.3, 0.6, 1)),
               0.025 * c(0.1, 0.3, 0.6, 1))
  expect_equal(gs_spend(0, 0.025, "OF"), 0)
})

test_that("gs_spend checks its input", {
  expect_error(gs_spend(0.5, 0.025, "HSD"), "param")
  expect_error(gs_spend(1.5, 0.025, "OF"), "between 0 and 1")
  expect_error(gs_spend(c(0.5, 1), 0.025, "user", c(0.6, 0.5)), "nondecreasing")
  expect_error(gs_spend(0.5, 1.2, "OF"), "alpha")
})
