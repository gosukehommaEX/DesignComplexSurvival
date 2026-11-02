test_that("plot.ms_power returns the plotted values invisibly", {
  x <- plot_design_b()
  d <- with_pdf(expect_invisible(plot(x)))
  expect_gt(attr(d, "pdf_size"), 0)
  expect_identical(d$upper, x$looks$z_upper)
  expect_identical(d$x, x$looks$events)
  d <- with_pdf(plot(x, "bounds", scale = "hr", x_axis = "time", main = "Bounds",
                     ylim = c(0.3, 1.5)))
  expect_identical(d$upper, x$looks$hr_upper)
  expect_identical(d$x, x$looks$time)
  d <- with_pdf(plot(x, type = "power"))
  expect_identical(d$cum_efficacy, x$looks$cum_power)
  d <- with_pdf(plot(x, type = "survival", max_time = 36))
  expect_equal(max(d$time), 36)
})

test_that("plot.ms_power draws only the survival functions of a design with two endpoints", {
  x <- plot_design_b()
  y <- x
  y$settings$endpoint <- c("pfs", "os")
  expect_error(plot(y), "one endpoint only")
  expect_error(plot(y, type = "power"), "one endpoint only")
  d <- with_pdf(plot(y, type = "survival"))
  d0 <- with_pdf(plot(x, type = "survival"))
  attr(d, "pdf_size") <- NULL
  attr(d0, "pdf_size") <- NULL
  expect_identical(d, d0)
})

test_that("plot.ms_power checks its input", {
  x <- plot_design_b()
  expect_error(plot(x, type = "other"), "should be one of")
  expect_error(plot(x, scale = "log"), "should be one of")
  expect_error(plot(x, x_axis = "look"), "should be one of")
  expect_error(plot(x, max_time = 30), "only with type = \"survival\"")
  expect_error(plot(x, type = "power", max_time = 30), "only with type = \"survival\"")
  expect_error(plot(x, "bounds", "z", "events", NULL, "main"), "must be named")
})
