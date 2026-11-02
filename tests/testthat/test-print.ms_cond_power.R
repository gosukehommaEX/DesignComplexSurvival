test_that("print.ms_cond_power shows the conditional power and the later looks", {
  te <- ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2, beta = 0.1)
  x <- ms_power(scn_b(), n = 300, looks = ms_looks("os", fraction = c(0.5, 0.75, 1)), test = te,
                events = 120)
  cp <- ms_cond_power(x, 1, 0.5, under = "trend")
  expect_output(print(cp), "Conditional power of a group sequential log-rank design")
  expect_output(print(cp), "Observed at look 1: z = 0.5")
  expect_output(print(cp), sprintf("current trend (hazard ratio %s)", format(signif(cp$hr, 4))),
                fixed = TRUE)
  expect_output(print(cp), sprintf("Conditional power: %s", format(signif(cp$cond_power, 4))),
                fixed = TRUE)
  expect_output(print(cp), "Without the futility bounds of the later interim looks")
  expect_output(print(cp), "cond_mean")
  expect_output(expect_invisible(print(cp)), "cum_upper")
  h <- ms_cond_power(x, 2, c(0.5, 1))
  expect_output(print(h), "Observed at looks 1 to 2: z = 0.5, 1")
  expect_output(print(h), "Distribution: design scenario")
  expect_output(print(ms_cond_power(x, 1, 0.5, under = "hr", hr = 0.75)), "hazard ratio 0.75")
  expect_output(print(ms_cond_power(x, 1, 0.5, under = "null")), "null hypothesis")
  s <- ms_cond_power(x, 1, 0.5, under = "scenario", scenario = x$scenario)
  expect_output(print(s), "given scenario")
})
