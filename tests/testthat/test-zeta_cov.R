# Expected values: inst/validation/python/results/zeta_cov.csv (scenario F,
# items PFS and OS at months 12 and 20; ref_zeta.zeta_cov with derivatives by
# central differences)

test_that(".zeta_cov matches the reference for calendar and event-driven times", {
  scn <- scn_f()
  acc <- .acc_setup(scn$accrual, 400)
  ex <- expected_csv("zeta_cov.csv")
  ep <- c("pfs", "pfs", "os", "os")
  tt <- c(12, 20, 12, 20)
  for (tr in c("calendar", "os")) {
    zc <- .zeta_cov(scn, acc, ep, tt, if (tr == "calendar") NULL else tr, "studentized")
    e <- ex[ex$trigger == tr, ]
    expect_equal(zc$sigma[cbind(e$row, e$col)], e$value, tolerance = 1e-8)
  }
})

test_that(".zeta_cov reduces to the null covariance under the null", {
  scn <- ms_null(scn_f())
  acc <- .acc_setup(scn$accrual, 400)
  ep <- c("pfs", "os")
  tt <- c(20, 20)
  zs <- .zeta_cov(scn, acc, ep, tt, "pfs", "studentized")
  zn <- .zeta_cov(scn, acc, ep, tt, "pfs", "null")
  za <- .zeta_cov(scn, acc, ep, tt, "pfs", "alternative")
  expect_equal(diag(zs$sigma), c(1, 1), tolerance = 1e-8)
  expect_equal(zs$sigma, zn$sigma, tolerance = 1e-8)
  expect_equal(zs$sigma, za$sigma, tolerance = 1e-8)
  expect_equal(zs$mu, c(0, 0), tolerance = 1e-12)
})
