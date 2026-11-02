test_that(".cp_moments gives the distributions of the statistics of the looks", {
  te <- ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2, beta = 0.1)
  x <- ms_power(scn_b(), n = 300, looks = ms_looks("os", fraction = c(0.5, 0.75, 1)),
                test = te, events = 120)
  lk <- x$looks
  d <- lk$events
  p <- x$scenario$prob
  m <- .cp_moments(x, 1L, 0.5, "design")
  expect_identical(m[c("mean", "sigma", "time", "events")],
                   list(mean = x$mean, sigma = x$sigma, time = lk$time, events = d))
  expect_identical(m$hr, NA_real_)
  m <- .cp_moments(x, 1L, 0.5, "null")
  expect_identical(m[c("mean", "sigma", "time", "events")],
                   list(mean = rep(0, 3), sigma = x$sigma_null, time = x$looks_null$time,
                        events = x$looks_null$events))
  corr_d <- sqrt(outer(d, d, pmin) / outer(d, d, pmax))
  m <- .cp_moments(x, 2L, 1, "trend")
  expect_equal(m$mean, sqrt(d / d[2]))
  expect_equal(m$sigma, corr_d)
  expect_equal(m$hr, exp(-1 / sqrt(p[[1]] * p[[2]] * d[2])))
  # the hazard ratio of the trend gives the same mean
  mh <- .cp_moments(x, 2L, 1, "hr", hr = m$hr)
  expect_equal(mh$mean, m$mean)
  expect_identical(mh$sigma, m$sigma)
  expect_identical(mh$hr, m$hr)
  expect_identical(mh$time, lk$time)
  # the design scenario given as the scenario reproduces the design
  ms <- .cp_moments(x, 1L, 0.5, "scenario", scenario = x$scenario)
  expect_equal(ms$mean, x$mean, tolerance = 1e-10)
  expect_equal(ms$sigma, x$sigma, tolerance = 1e-10)
  expect_equal(ms$time, lk$time, tolerance = 1e-10)
  expect_equal(ms$events, d, tolerance = 1e-10)
  # the null scenario given as the scenario gives the looks and correlations of the type I error
  m0 <- .cp_moments(x, 1L, 0.5, "scenario", scenario = ms_null(x$scenario))
  expect_equal(m0$mean, rep(0, 3), tolerance = 1e-12)
  expect_equal(m0$sigma, x$sigma_null, tolerance = 1e-10)
  expect_equal(m0$time, x$looks_null$time, tolerance = 1e-10)
})
