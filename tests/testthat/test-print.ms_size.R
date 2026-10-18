test_that("print.ms_size shows the solved quantity", {
  res <- ms_size(scn_a(), looks = ms_looks("os"), test = ms_test(), power = 0.9,
                 solve_for = "events", n = 500)
  expect_output(print(res), "Solved for the final number of events")
  expect_output(print(res), "Power:")
})
