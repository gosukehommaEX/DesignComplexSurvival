# Expected values: inst/validation/python/results/calendar_group.csv (forward clock: matrix
# exponentials over the whole entry range; reset clock: adaptive quadrature for the OS)

test_that(".cal_group matches the reference (scenario C, switching from month 15)", {
  ex <- expected_csv("calendar_group.csv")
  acc <- .acc_setup(scn_cal()$accrual, 500)
  for (clock in c("forward", "reset")) {
    m <- models_c(clock)
    e <- ex[ex$clock == clock, ]
    for (j in unique(e$group)) {
      ej <- e[e$group == j, ]
      got <- .cal_group(if (j == 0) m$control else m$treatment, acc, 36, ej$s, 15, 2)
      expect_lt(max(abs(got$Y - ej$Y)), 1e-9)
      expect_lt(max(abs(got$D - ej$D)), 1e-10)
    }
  }
})

test_that(".cal_group agrees with the integral over the whole entry range", {
  m <- models_c()$control
  acc <- .acc_setup(scn_cal()$accrual, 500)
  s <- c(2, 9, 14, 21, 33)
  got <- .cal_group(m, acc, 36, s, 15, 2)
  for (k in seq_along(s)) {
    g <- .gl_pieces(0, min(36 - s[k], 18), c(6, 15 - c(0, 3, 4), 15 - s[k] + c(0, 3)),
                    hmax = 0.5, n = 12L)
    av <- .os_avail(m, rep(s[k], length(g$x)), pmax(15 - g$x, 0))
    w <- g$w * .pw_value(g$x, c(0, 6), acc$rho)
    expect_equal(got$Y[k], sum(w * av$surv), tolerance = 1e-10)
    expect_equal(got$D[k], sum(w * av$dens), tolerance = 1e-10)
  }
})

test_that(".cal_group of a group without switching is S(s) A(t - s)", {
  m <- models_c()$treatment
  acc <- .acc_setup(scn_cal()$accrual, 500)
  s <- c(2, 20, 30)
  ep <- .ep_eval(m, s, "os")
  got <- .cal_group(m, acc, 36, s, 15, 2)
  expect_equal(got$Y, ep$surv * .acc_cum(acc, 36 - s))
  expect_equal(got$D, ep$dens * .acc_cum(acc, 36 - s))
})

test_that(".cal_group gives the same values with a much finer rule (reset clock)", {
  # the forward clock does not use the rule over the entry time (test-cal_mid_fwd.R)
  acc <- .acc_setup(scn_cal()$accrual, 500)
  s <- c(2.5, 7.3, 12.9, 19.1)
  for (nm in c("reset_same", "reset_switch", "reset_aft0.6", "reset_aft1.5")) {
    m <- model_post_cal(nm)
    a <- .cal_group(m, acc, 33.1, s, 13.7, 20)
    b <- .cal_group(m, acc, 33.1, s, 13.7, 0.25)
    expect_lt(max(abs(c(a$Y - b$Y, a$D - b$D))), 1e-12)
  }
})
