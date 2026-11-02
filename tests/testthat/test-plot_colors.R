test_that(".plot_colors takes the colors from the Okabe-Ito palette", {
  col <- .plot_colors()
  oi <- unname(grDevices::palette.colors(palette = "Okabe-Ito"))
  expect_identical(names(col), c("efficacy", "futility", "mean", "alpha", "control",
                                 "treatment"))
  expect_true(all(col %in% oi))
  expect_identical(col[["efficacy"]], "#0072B2")
  expect_identical(col[["futility"]], "#D55E00")
  expect_false(col[["control"]] == col[["treatment"]])
})
