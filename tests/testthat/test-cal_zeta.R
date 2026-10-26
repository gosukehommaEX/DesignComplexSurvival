# Expected values: inst/validation/python/results/calendar_zeta.csv (ref_calendar.py: moments of
# the per-patient functionals given the entry time, equation (5) of the covariance note, by
# Gauss-Legendre rules with the OS from matrix exponentials; time derivatives by central
# differences)

test_that(".cal_zeta matches the reference under the alternative", {
  scn <- scn_cal()
  acc <- .acc_setup(scn$accrual, 500)
  z <- .cal_zeta(scn, acc, c(24, 36), NULL, "studentized")
  expect_lt(max(abs(z$sigma - expected_zeta_cal("alt", "none", "studentized"))), 1e-7)
  ex <- expected_csv("calendar_looks.csv")
  expect_equal(z$mu, ex$mu[ex$design == "alt"], tolerance = 1e-9)
  expect_equal(z$v, ex$v[ex$design == "alt"], tolerance = 1e-9)
  z <- .cal_zeta(scn, acc, c(24, 36), c("os", "os"), "studentized")
  expect_lt(max(abs(z$sigma - expected_zeta_cal("alt", "os", "studentized"))), 1e-7)
  z <- .cal_zeta(scn, acc, c(24, 36), NULL, "alternative")
  expect_lt(max(abs(z$sigma - expected_zeta_cal("alt", "none", "alternative"))), 1e-7)
})

test_that(".cal_zeta matches the reference with the null variance", {
  scn <- scn_cal()
  acc <- .acc_setup(scn$accrual, 500)
  z <- .cal_zeta(scn, acc, c(24, 36), c("os", "os"), "null")
  expect_lt(max(abs(z$sigma - expected_zeta_cal("alt", "os", "null"))), 1e-7)
  expect_equal(diag(z$sigma), c(1, 1))
  # under the null hypothesis the three variances agree (correlation of calendar_delta.csv)
  s0 <- scn_cal(null = TRUE)
  z0 <- .cal_zeta(s0, acc, c(24, 36), NULL, "studentized")
  expect_lt(max(abs(z0$sigma - expected_zeta_cal("null", "none", "studentized"))), 1e-7)
  expect_lt(max(abs(z0$sigma - expected_zeta_cal("null", "none", "null"))), 1e-7)
  expect_equal(z0$sigma[1, 2], expected_csv("calendar_delta.csv")$corr, tolerance = 1e-7)
})

test_that(".cal_zeta agrees with .zeta_cov when everybody can switch", {
  acc <- .acc_setup(scn_cal()$accrual, 500)
  z <- .cal_zeta(scn_cal(tc = 1e-8), acc, c(24, 36), c("os", NA), "studentized")
  ref <- .zeta_cov(scn_cal(tc = 0), acc, c("os", "os"), c(24, 36), c("os", NA),
                   "studentized")$sigma
  expect_lt(max(abs(z$sigma - ref)), 1e-9)
})

test_that(".cal_zeta agrees with .zeta_cov when nobody switches by the last look", {
  acc <- .acc_setup(scn_cal()$accrual, 500)
  m <- models_c("forward")
  m0 <- m$control
  m0$switch_prob[] <- 0
  m0$switch_rate[] <- 0
  no_sw <- ms_scenario(m0, m$treatment, accrual = scn_cal()$accrual)
  z <- .cal_zeta(scn_cal(tc = 40), acc, c(24, 36), c(NA, "os"), "studentized")
  ref <- .zeta_cov(no_sw, acc, c("os", "os"), c(24, 36), c(NA, "os"), "studentized")$sigma
  expect_lt(max(abs(z$sigma - ref)), 1e-7)
  expect_error(.cal_zeta(scn_cal(), acc, 24, "pfs"), "Internal error")
})
