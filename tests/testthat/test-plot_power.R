test_that(".plot_power draws the cumulative stopping probabilities", {
  x <- plot_design_b()
  lk <- x$looks
  d <- with_pdf(.plot_power(x, "events", list()))
  expect_gt(attr(d, "pdf_size"), 0)
  expect_identical(d$x, lk$events)
  expect_identical(d$cum_efficacy, lk$cum_power)
  expect_equal(d$cum_efficacy, cumsum(lk$prob_upper))
  expect_equal(d$cum_futility, c(lk$prob_lower[1], sum(lk$prob_lower[1:2]),
                                 sum(lk$prob_lower[1:2])))
  expect_identical(d$x_null, x$looks_null$events)
  expect_identical(d$cum_alpha, lk$cum_alpha)
  expect_equal(d$cum_alpha[3], x$type1[["nonbinding"]])
  d <- with_pdf(.plot_power(x, "time", list(main = "Stopping")))
  expect_identical(d$x, lk$time)
  expect_identical(d$x_null, x$looks_null$time)
})

test_that(".plot_power needs the looks under the null scenario", {
  x <- plot_design_b()
  x$looks_null <- NULL
  expect_error(with_pdf(.plot_power(x, "events", list())), "compute it again")
})
