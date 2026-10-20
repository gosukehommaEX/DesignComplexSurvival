test_that(".design_core2 orders the statistics by endpoint and look", {
  d <- design_def("D3")
  res <- .design_core2(scn_f(), 400, d$looks, d$test, events = 220, type1 = FALSE)
  expect_identical(res$looks$endpoint, c("pfs", "pfs", "os", "os"))
  expect_equal(res$mean, res$looks$z_mean)
  expect_equal(sqrt(diag(res$sigma)), res$looks$z_sd)
  expect_true(all(is.na(res$type1)))
  expect_null(res$fwer)
  expect_null(res$looks$cum_alpha)
  lc <- lr_cov(scn_f(), n = 400, calendar_time = unique(res$looks$time),
               endpoint = c("pfs", "os"), variance = "alternative")
  expect_equal(res$sigma, lc$sigma, tolerance = 1e-12)
  expect_equal(res$mean, lc$mean, tolerance = 1e-12)
})
