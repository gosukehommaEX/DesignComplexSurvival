# Expected values: inst/validation/python/results/model5_states.csv

test_that(".pre_states matches the reference values", {
  ex <- expected_csv("model5_states.csv")
  ex <- ex[ex$pfs_switch == "event" & ex$clock == "reset", ]
  for (name in c("P", "P2", "Z")) {
    e <- ex[ex$model == name, ]
    ps <- .pre_states(models_s3(name), e$time)
    expect_lt(max(abs(ps$p1 - e$p1)), 1e-12)
    expect_lt(max(abs(ps$p5 - e$p5)), 1e-12)
  }
})

test_that(".pre_states is exponential without a switch and stays finite for large hazards", {
  m <- models_c()$control
  ps <- .pre_states(m, c(0, 3, 10))
  expect_identical(ps$p1, exp(-.pw_cumhaz(c(0, 3, 10), m$time, m$prog + m$death)))
  expect_identical(ps$p5, c(0, 0, 0))
  # hazards per day over years: no overflow, and the constant-hazard closed form
  big <- ms_model(prog = 5, death = 0.01, death_prog = 1, switch_pre = 40, prog_switch = 0.001,
                  death_switch_pre = 0.0001)
  s <- c(0.5, 10, 900)
  ps <- .pre_states(big, s)
  h1 <- 45.01
  h5 <- 0.0011
  expect_true(all(is.finite(ps$p5)))
  expect_equal(ps$p5, 40 / (h1 - h5) * (exp(-h5 * s) - exp(-h1 * s)), tolerance = 1e-13)
})
