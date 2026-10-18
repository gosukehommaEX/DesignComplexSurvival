test_that("as.data.frame.ms_power returns the per-look table", {
  res <- ms_power(scn_a(), n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
                  test = ms_test(), events = 300)
  df <- as.data.frame(res)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 2L)
  expect_equal(df[["cum_power"]][2], res$power)
})
