# Functionals of .cal_zeta() at one look (psi, phi and dlt as columns)
cal_funs <- function(lk) {
  lapply(1:2, function(j) {
    gj <- (j - 1) - lk$q1
    rho <- cbind(psi = gj, phi = lk$q1 * (1 - lk$q1), dlt = 1)
    beta <- cbind(psi = -gj * lk$lbar, phi = (1 - 2 * lk$q1) * gj * lk$lbar, dlt = 0)
    list(rho = rho, beta = beta, pr = .pw_poly(lk$g, rho, 12L), pb = .pw_poly(lk$g, beta, 12L))
  })
}

test_that(".cal_pair gives the null variance and the closed form of the null covariance", {
  # Expected value: inst/validation/python/results/calendar_delta.csv (delta of equation (3) of
  # the calendar-time switching note, in units of events; times p0 p1 = 1/4 for psi)
  scn <- scn_cal(null = TRUE)
  acc <- .acc_setup(scn$accrual, 500)
  l1 <- .cal_look(scn, acc, 24)
  l2 <- .cal_look(scn, acc, 36)
  f1 <- cal_funs(l1)
  p11 <- .cal_pair(scn, acc, l1, l1, f1, f1)
  expect_equal(p11, base::t(p11), tolerance = 1e-12)
  expect_equal(p11[["psi", "psi"]], l1$v, tolerance = 1e-9)
  expect_equal(p11[["dlt", "dlt"]], sum(l1$ev), tolerance = 1e-10)
  p12 <- .cal_pair(scn, acc, l1, l2, f1, cal_funs(l2))
  dl <- expected_csv("calendar_delta.csv")
  expect_equal(p12[["psi", "psi"]], l1$v + dl$delta / 4, tolerance = 1e-9)
  # an event by the first look is an event by the second look
  expect_equal(p12[["dlt", "dlt"]], sum(l1$ev), tolerance = 1e-10)
})

test_that(".cal_pair gives the probability of an event under the alternative", {
  scn <- scn_cal()
  acc <- .acc_setup(scn$accrual, 500)
  l1 <- .cal_look(scn, acc, 24)
  l2 <- .cal_look(scn, acc, 36)
  p12 <- .cal_pair(scn, acc, l1, l2, cal_funs(l1), cal_funs(l2))
  expect_equal(p12[["dlt", "dlt"]], sum(l1$ev), tolerance = 1e-10)
  expect_equal(dimnames(p12), list(c("psi", "phi", "dlt"), c("psi", "phi", "dlt")))
})
