# Expected values: inst/validation/python/results/size.csv

test_that("ms_size reproduces the Schoenfeld number of events for a fixed design", {
  res <- ms_size(scn_a(), looks = ms_looks("os"), test = ms_test(mean = "schoenfeld"),
                 power = 0.9, solve_for = "events", n = 500, round_up = FALSE)
  # (z_0.975 + z_0.9)^2 / (0.25 log(0.7)^2)
  d_ref <- (qnorm(0.975) + qnorm(0.9))^2 / (0.25 * log(0.7)^2)
  expect_equal(d_ref, 330.3779139638097, tolerance = 1e-12)
  expect_lt(abs(res$solution - d_ref) / d_ref, 1e-6)
  expect_lt(abs(res$power - 0.9), 1e-6)
  expect_s3_class(res, "ms_power")
})

test_that("ms_size matches the reference two-look design and rounds up", {
  lo <- ms_looks("os", fraction = c(0.5, 1))
  res <- ms_size(scn_a(), looks = lo, test = ms_test(), power = 0.9,
                 solve_for = "events", n = 500)
  expect_lt(abs(res$solution - 332.2830834262169), 0.01)
  expect_equal(res$events, ceiling(res$solution))
  expect_gte(res$power, 0.9)
})

test_that("solving for the follow-up gives the same final analysis time", {
  lo <- ms_looks("os", fraction = c(0.5, 1))
  res <- ms_size(scn_a(), looks = lo, test = ms_test(), power = 0.9,
                 solve_for = "followup", n = 500)
  # accrual of 500 patients at 20 per month ends at month 25
  expect_lt(abs(res$solution + 25 - 38.73277391908485), 1e-4)
  expect_lt(abs(res$power - 0.9), 1e-6)
})

test_that("ms_size solves for the sample size", {
  lo <- ms_looks("os", fraction = c(0.5, 1))
  te <- ms_test()
  res <- ms_size(scn_a(), looks = lo, test = te, power = 0.9, solve_for = "n",
                 followup = 12)
  expect_equal(res$n, ceiling(res$solution))
  expect_gte(res$power, 0.9)
  p_sol <- ms_power(scn_a(), n = res$solution, looks = lo, test = te, followup = 12)$power
  expect_lt(abs(p_sol - 0.9), 1e-6)
  p_less <- ms_power(scn_a(), n = res$n - 1, looks = lo, test = te, followup = 12)$power
  expect_lt(p_less, 0.9)
})

test_that("ms_size solves for the accrual duration", {
  s <- scn_a()
  s$accrual <- ms_accrual(prop = 1, duration = 20)
  lo <- ms_looks("os", fraction = c(0.5, 1))
  res <- ms_size(s, looks = lo, test = ms_test(), power = 0.85,
                 solve_for = "accrual_duration", n = 400, followup = 6)
  expect_lt(abs(res$power - 0.85), 1e-6)
  expect_equal(res$accrual_duration, res$solution)
})

test_that("ms_size checks its input", {
  lc <- ms_looks("calendar", time = c(24, 36))
  expect_error(ms_size(scn_a(), lc, ms_test(), solve_for = "events", n = 400),
               "event-driven")
  expect_error(ms_size(scn_a(), ms_looks("os"), ms_test(), solve_for = "n"), "followup")
  expect_error(ms_size(scn_a(), ms_looks("os"), ms_test(), solve_for = "accrual_duration",
                       n = 400, followup = 6), "prop")
  expect_error(ms_size(scn_a(), ms_looks("os"), ms_test(), solve_for = "events"), "'n'")
})
