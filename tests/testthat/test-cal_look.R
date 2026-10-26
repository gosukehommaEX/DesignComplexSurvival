# Expected values: inst/validation/python/results/calendar_looks.csv (ref_calendar.py: group
# functions from matrix exponentials over the entry time, then Gauss-Legendre rules over the
# follow-up time)

test_that(".cal_look matches the reference (scenario C, switching from month 15)", {
  ex <- expected_csv("calendar_looks.csv")
  for (d in c("alt", "null")) {
    scn <- scn_cal(null = d == "null")
    acc <- .acc_setup(scn$accrual, 500)
    for (t in c(24, 36)) {
      e <- ex[ex$design == d & ex$t == t, ]
      lk <- .cal_look(scn, acc, t)
      expect_lt(max(abs(c(lk$ev, lk$do) - c(e$ev0, e$ev1, e$do0, e$do1))), 1e-10)
      expect_lt(abs(lk$mu - e$mu), 1e-10)
      expect_lt(abs(lk$v - e$v), 1e-10)
    }
  }
})

test_that(".cal_look covers the follow-up time and gives the risk-set functions", {
  scn <- scn_cal(null = TRUE)
  acc <- .acc_setup(scn$accrual, 500)
  lk <- .cal_look(scn, acc, 30)
  expect_equal(sum(lk$g$w), 30)
  expect_true(all(c(15, 24) %in% lk$g$lo))
  # identical groups: half of the risk set in each, no score mean
  expect_equal(lk$q1, rep(0.5, length(lk$g$x)))
  expect_identical(lk$Y[[1L]], lk$Y[[2L]])
  expect_equal(lk$mu, 0)
  expect_equal(lk$lbar, lk$D[[1L]] / lk$Y[[1L]], tolerance = 1e-14)
})

test_that(".cal_look agrees with a finer rule over the follow-up time", {
  # control model "forward_switch" of calendar_os.csv, treatment group of scenario C; without
  # the points of .cal_cross() the rule of .cal_look() would be off by about 1.6e-8 (the reset
  # clock, slower here, is checked with the Python port in chk_tests_calendar.py)
  m <- model_post_cal("forward_switch")
  scn <- ms_scenario(m, models_c("forward")$treatment, accrual = scn_cal()$accrual,
                     switch_start = 13.7)
  acc <- .acc_setup(scn$accrual, 500)
  lk <- .cal_look(scn, acc, 26.3)
  g <- .gl_pieces(0, 26.3, c(.scn_breaks(scn), 26.3 - acc$alpha, 26.3 - 13.7,
                             13.7 - acc$alpha, 13.7), hmax = 0.5, n = 12L)
  gr <- .cal_group(m, acc, 26.3, g$x, 13.7, lk$hmax)
  ev0 <- sum(g$w * scn$prob[[1L]] * exp(-.pw_cumhaz(g$x, m$time, m$dropout)) * gr$D)
  expect_lt(abs(lk$ev[1L] - ev0), 1e-10)
})
