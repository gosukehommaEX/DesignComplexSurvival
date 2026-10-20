test_that(".target_power selects the targeted power", {
  res <- list(power = c(pfs = 0.9, os = 0.7, both = 0.65, any = 0.95))
  expect_equal(.target_power(res, "first"), 0.9)
  expect_equal(.target_power(res, "each"), 0.7)
  expect_equal(.target_power(res, "both"), 0.65)
  expect_equal(.target_power(res, "any"), 0.95)
  expect_equal(.target_power(list(power = 0.8), "first"), 0.8)
})
