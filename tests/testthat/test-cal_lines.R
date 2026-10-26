test_that(".cal_lines lists the lines of scenario C and none without switching", {
  ln <- .cal_lines(models_c("forward")$control)
  expect_identical(colnames(ln), c("c0", "c1"))
  # forward clock: o = 0, 4, 3 (breaks of time and time_post) and o = s
  expect_equal(unname(ln), cbind(c(0, 4, 3, 0), c(0, 0, 0, 1)))
  expect_equal(nrow(.cal_lines(models_c("forward")$treatment)), 0L)
  # reset clock: o = 0, 3, 4, 7 and o = s - 3, s, s + 3
  lr <- unname(.cal_lines(models_c("reset")$control))
  expect_equal(lr[order(lr[, 2], lr[, 1]), ], cbind(c(0, 3, 4, 7, -3, 0, 3), c(0, 0, 0, 0, 1, 1, 1)))
})

test_that("the OS density given the availability has a kink on a line and not off it", {
  # difference of the one-sided differences in o (step h): the jump of the derivative on a line,
  # of order h off the lines
  slope_jump <- function(m, s, o, h = 1e-3) {
    d <- .os_avail(m, rep(s, 3L), o + c(-h, 0, h))$dens
    (d[3L] - d[2L]) / h - (d[2L] - d[1L]) / h
  }
  cases <- list(list("forward_aft1.5", 2.5, 1.5, 1.2), list("forward_switch", 12, 9, 8.5),
                list("reset_same", 12, 9, 8.5), list("reset_aft0.6", 12, 7, 6.5))
  for (cs in cases) {
    m <- model_post_cal(cs[[1L]])
    ln <- .cal_lines(m)
    on <- abs(ln[, 1L] + ln[, 2L] * cs[[2L]] - cs[[3L]]) < 1e-12
    off <- abs(ln[, 1L] + ln[, 2L] * cs[[2L]] - cs[[4L]]) < 1e-12
    expect_true(any(on))
    expect_false(any(off))
    expect_gt(abs(slope_jump(m, cs[[2L]], cs[[3L]])), 1e-5)
    expect_lt(abs(slope_jump(m, cs[[2L]], cs[[4L]])), 1e-6)
  }
})
