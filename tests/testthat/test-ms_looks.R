test_that("ms_looks defines event-driven and calendar-driven looks", {
  l <- ms_looks("os", fraction = c(0.5, 0.75, 1))
  expect_s3_class(l, "ms_looks")
  expect_equal(l$k, 3L)
  expect_identical(l$trigger, "os")
  expect_equal(ms_looks("pfs")$fraction, 1)
  lc <- ms_looks("calendar", time = c(24, 36))
  expect_equal(lc$time, c(24, 36))
  expect_null(lc$fraction)
})

test_that("ms_looks checks its input", {
  expect_error(ms_looks("os", fraction = c(0.5, 0.9)), "end at 1")
  expect_error(ms_looks("os", fraction = c(0.7, 0.5, 1)), "increasing")
  expect_error(ms_looks("calendar"), "time")
  expect_error(ms_looks("calendar", time = c(36, 24)), "increasing")
  expect_error(ms_looks("os", time = 24), "fraction")
  expect_error(ms_looks("os", fraction = (1:25) / 25), "20 looks")
})
