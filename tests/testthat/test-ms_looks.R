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

test_that("ms_looks defines one trigger per look and the tested endpoints", {
  l <- ms_looks(c("pfs", "pfs", "os"), fraction = c(0.6, 1, 1),
                test = list("pfs", c("pfs", "os"), "os"))
  expect_s3_class(l, "ms_looks")
  expect_equal(l$k, 3L)
  expect_identical(l$trigger, c("pfs", "pfs", "os"))
  expect_equal(l$fraction, c(0.6, 1, 1))
  expect_equal(l$time, c(NA_real_, NA_real_, NA_real_))
  expect_identical(l$test, list("pfs", c("pfs", "os"), "os"))
  lc <- ms_looks(c("calendar", "os"), fraction = c(NA, 1), time = c(24, NA))
  expect_equal(lc$fraction, c(NA, 1))
  expect_equal(lc$time, c(24, NA))
  expect_identical(ms_looks("os", fraction = c(0.5, 1), test = c("os", "pfs"))$test,
                   list("os", "pfs"))
  expect_null(ms_looks("os")$test)
})

test_that("ms_looks checks one trigger per look", {
  expect_error(ms_looks(c("os", "dfs")), "trigger")
  expect_error(ms_looks(c("pfs", "os"), fraction = c(1, NA)), "fraction")
  expect_error(ms_looks(c("calendar", "os"), fraction = c(0.5, 1), time = c(24, NA)),
               "fraction")
  expect_error(ms_looks(c("calendar", "os"), fraction = c(NA, 1), time = c(24, 30)), "time")
  expect_error(ms_looks(c("pfs", "pfs", "os"), fraction = c(1, 0.6, 1)), "PFS")
  expect_error(ms_looks(c("pfs", "os"), fraction = c(1, 0.9)), "final look")
  expect_error(ms_looks(c("os", "pfs"), fraction = c(1, 1.2)), "\\(0, 1\\]")
  expect_error(ms_looks(c("calendar", "calendar"), time = c(36, 24)), "increasing")
  expect_error(ms_looks(rep(c("pfs", "os"), 11), fraction = rep(1, 22)), "20 looks")
  expect_error(ms_looks("os", fraction = c(0.5, 1), test = list("os")), "test")
  expect_error(ms_looks("os", fraction = c(0.5, 1), test = list("os", "dfs")), "test")
  expect_error(ms_looks("os", fraction = c(0.5, 1), test = list("os", c("os", "os"))),
               "test")
})
