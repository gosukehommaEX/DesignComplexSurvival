test_that("print.ms_size shows the solved quantity", {
  res <- ms_size(scn_a(), looks = ms_looks("os"), test = ms_test(), power = 0.9,
                 solve_for = "events", n = 500)
  expect_output(print(res), "Solved for the final number of events")
  expect_output(print(res), "Power:")
})

test_that("print.ms_size shows the final numbers of events of two endpoints", {
  x <- run_design("D3")
  x$solve_for <- "events"
  x$solution <- c(pfs = 250.4, os = 180.2)
  x$target_power <- 0.8
  x$target <- "first"
  class(x) <- c("ms_size", "ms_power")
  expect_output(print(x),
                "final numbers of events: PFS 250.4, OS 180.2 \\(target power of PFS 0.8\\)")
  expect_output(print(x), "two endpoints")
})
