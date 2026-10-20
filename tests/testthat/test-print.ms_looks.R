test_that("print.ms_looks describes the looks", {
  expect_output(print(ms_looks("os", fraction = c(0.5, 1))), "driven by OS events")
  expect_output(print(ms_looks("calendar", time = c(24, 36))), "calendar times")
})

test_that("print.ms_looks describes one trigger per look", {
  l <- ms_looks(c("calendar", "os"), fraction = c(NA, 1), time = c(15, NA),
                test = list("pfs", c("pfs", "os")))
  expect_output(print(l), "at calendar time 15; tests PFS")
  expect_output(print(l), "driven by OS events at fraction 1; tests PFS and OS")
  expect_output(print(ms_looks(c("pfs", "os"), fraction = c(1, 1))), "tests all endpoints")
})
