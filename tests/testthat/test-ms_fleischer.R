test_that("ms_fleischer gives exponential PFS and OS with the given medians", {
  m <- ms_fleischer(median_pfs = 6, median_os = 15)
  tt <- c(1, 6, 15, 30)
  pfs <- ms_surv(m, tt, endpoint = "pfs")
  os <- ms_surv(m, tt, endpoint = "os")
  # closed forms: PFS hazard log(2) / 6, OS hazard log(2) / 15
  expect_lt(rel_diff(pfs$surv, exp(-log(2) * tt / 6)), 1e-12)
  expect_lt(rel_diff(os$surv, exp(-log(2) * tt / 15)), 1e-10)
  expect_lt(rel_diff(os$haz, rep(log(2) / 15, 4)), 1e-9)
  expect_equal(ms_surv(m, 6, "pfs")$surv, 0.5, tolerance = 1e-12)
  expect_equal(ms_surv(m, 15, "os")$surv, 0.5, tolerance = 1e-10)
})

test_that("ms_fleischer checks the medians", {
  expect_error(ms_fleischer(median_pfs = 10, median_os = 5), "median_os")
  expect_error(ms_fleischer(median_pfs = -1, median_os = 5), "positive")
  expect_equal(ms_fleischer(median_pfs = 5, median_os = 5)$prog, 0)
})
