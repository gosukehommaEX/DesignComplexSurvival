test_that("summary.ms_power returns the overview and the looks", {
  res <- ms_power(scn_a(), n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
                  test = ms_test(), events = 300)
  s <- summary(res)
  expect_s3_class(s, "summary.ms_power")
  expect_equal(s$overview[["power"]], res$power)
  expect_equal(s$overview[["events"]], res$events)
  expect_identical(s$looks, res$looks)
})

test_that("summary.ms_power returns the overview of a design with two endpoints", {
  res <- run_design("D3")
  s <- summary(res)
  expect_s3_class(s, "summary.ms_power")
  expect_equal(s$overview[["power_pfs"]], res$power[["pfs"]])
  expect_equal(s$overview[["power_any"]], res$power[["any"]])
  expect_equal(s$overview[["events_os"]], res$events[["os"]])
  expect_equal(s$overview[["fwer_nonbinding"]], res$type1[["nonbinding"]])
  expect_identical(s$fwer, res$fwer)
})
