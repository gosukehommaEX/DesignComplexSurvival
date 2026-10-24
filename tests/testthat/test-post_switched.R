test_that(".post_switched gives the survival with the hazard death_switch", {
  m <- models_c("forward")$control
  fr <- c(0, 1, 2.5)
  to <- c(2, 6, 25)
  p <- .post_switched(m, fr, to)
  tr <- .post_trans(m, fr, to)
  expect_equal(p$surv, tr$r33, tolerance = 1e-15)
  expect_equal(p$dens, tr$r33 * tr$d_to, tolerance = 1e-15)
  # with a switch clock or an AFT effect the hazard is still death_switch on the clock
  w <- models_s3("W", "forward")
  expect_equal(.post_switched(w, fr, to), p)
  a <- models_s3("A06")
  expect_equal(.post_switched(a, 0, 4)$surv, exp(-0.6 * (0.05 * 3 + 0.06)), tolerance = 1e-15)
})
