# Expected values: inst/validation/python/results/model5_post.csv (ref_model5.post:
# adaptive quadrature over the time of the switch)

test_that(".post_general matches the reference values", {
  ex <- expected_csv("model5_post.csv")
  for (name in unique(ex$model)) {
    e <- ex[ex$model == name, ]
    g <- .post_general(model_post_s3(name), e$from, e$to)
    expect_lt(max(abs(g$r22 - e$r22)), 1e-12)
    expect_lt(max(abs(g$r23 - e$r23)), 1e-12)
    expect_lt(max(abs(g$r33 - e$r33)), 1e-12)
    expect_lt(max(abs(g$f23 - e$f23)), 1e-12)
    expect_lt(max(abs(g$r33 * g$d_to - e$f33)), 1e-12)
  }
})

test_that(".post_general reduces to the closed forms of the current model", {
  fr <- c(0, 0, 0, 1.5, 2)
  to <- c(1, 4, 20, 4, 7.5)
  args <- list(death = 0.01, time_post = c(0, 3), death_prog = c(0.05, 0.08),
               switch_rate = c(0.04, 0.06))
  # constant death_switch: the clock of the switch does not matter
  a <- .post_general(do.call(ms_model, c(args, death_switch = 0.03, switch_clock = "switch")),
                     fr, to)
  b <- .post_trans(do.call(ms_model, c(args, death_switch = 0.03)), fr, to)
  expect_lt(max(abs(unlist(a) - unlist(b[names(a)]))), 1e-15)
  # constant death_prog: the AFT effect is a proportional hazards effect
  args$death_prog <- 0.05
  a <- .post_general(do.call(ms_model, c(args, switch_af = 0.6)), fr, to)
  b <- .post_trans(do.call(ms_model, c(args, death_switch = 0.6 * 0.05)), fr, to)
  expect_lt(max(abs(unlist(a) - unlist(b[names(a)]))), 1e-15)
  # acceleration factor 1: no effect of the switch (default closed form)
  m1 <- do.call(ms_model, c(args, switch_af = 1))
  expect_equal(.post_trans(m1, fr, to),
               .post_trans(do.call(ms_model, args), fr, to))
})

test_that(".post_general stays finite where the default closed form overflows", {
  m0 <- ms_model(prog = 0.1, death = 0.01, death_prog = 1, switch_rate = 2, death_switch = 0.5)
  expect_error(.post_trans(m0, 400, 410), "larger time unit")
  m <- ms_model(prog = 0.1, death = 0.01, death_prog = 1, switch_rate = 2, death_switch = 0.5,
                switch_clock = "switch")
  g <- .post_general(m, c(400, 0), c(410, 10))
  expect_equal(g$r23, rep(2 / 2.5 * (exp(-5) - exp(-30)), 2), tolerance = 1e-14)
  expect_equal(g$f23, 0.5 * g$r23, tolerance = 1e-14)
})
