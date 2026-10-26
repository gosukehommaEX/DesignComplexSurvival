# Expected values: inst/validation/python/results/calendar_os.csv (adaptive quadrature over the
# progression time with the post-progression quantities of ref_model5.post(), itself adaptive
# quadrature over the switch time; ref_calendar.os_quad())

test_that(".os_avail matches the reference for eight post-progression models", {
  ex <- expected_csv("calendar_os.csv")
  for (nm in unique(ex$model)) {
    e <- ex[ex$model == nm, ]
    got <- .os_avail(model_post_cal(nm), e$s, e$o)
    expect_lt(max(abs(got$surv - e$surv)), 1e-10)
    expect_lt(max(abs(got$dens - e$dens)), 1e-11)
  }
})

test_that(".os_avail reduces to the model and to the model without switching", {
  s <- c(0.5, 4, 12.3, 25, 40)
  for (nm in c("reset_same", "forward_switch", "reset_aft0.6")) {
    m <- model_post_cal(nm)
    ep <- .ep_eval(m, s, "os")
    av <- .os_avail(m, s, numeric(length(s)))
    expect_lt(max(abs(av$surv - ep$surv)), 1e-12)
    expect_lt(max(abs(av$dens - ep$dens)), 1e-13)
    m0 <- m
    m0$switch_prob[] <- 0
    m0$switch_rate[] <- 0
    e0 <- .ep_eval(m0, s, "os")
    av <- .os_avail(m, s, s + 1)
    expect_lt(max(abs(av$surv - e0$surv)), 1e-12)
    expect_lt(max(abs(av$dens - e0$dens)), 1e-13)
  }
  expect_equal(.os_avail(m, numeric(0), numeric(0)), list(surv = numeric(0), dens = numeric(0)))
})

test_that(".os_avail handles the pairs in chunks", {
  m <- model_post_cal("reset_switch")
  s <- seq(1, 30, length.out = 9)
  o <- seq(12, 0, length.out = 9)
  expect_equal(.os_avail(m, s, o, chunk = 2L), .os_avail(m, s, o), tolerance = 1e-14)
})
