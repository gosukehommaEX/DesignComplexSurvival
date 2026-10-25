test_that("ms_model expands the intensities and sets defaults", {
  m <- ms_model(time = c(0, 6), prog = c(0.1, 0.08), death = 0.01,
                death_prog = 0.05, switch_prob = 0.3, time_post = c(0, 3))
  expect_s3_class(m, "ms_model")
  expect_equal(m$death, c(0.01, 0.01))
  expect_equal(m$switch_prob, c(0.3, 0.3))
  expect_equal(m$death_prog, c(0.05, 0.05))
  expect_equal(m$death_switch, m$death_prog)
  expect_equal(m$switch_rate, c(0, 0))
  expect_equal(m$dropout, c(0, 0))
  expect_identical(m$clock, "reset")
  expect_identical(ms_model(death = 0.1, clock = "forward")$clock, "forward")
})

test_that("ms_model without progression needs no post-progression hazard", {
  m <- ms_model(death = 0.05)
  expect_equal(m$prog, 0)
  expect_equal(m$death_prog, 0)
})

test_that("ms_model rejects invalid input", {
  expect_error(ms_model(prog = 0.1, death = 0.01), "death_prog")
  expect_error(ms_model(), "death")
  expect_error(ms_model(time = c(1, 2), death = 0.1), "start at 0")
  expect_error(ms_model(time = c(0, 2, 1), death = 0.1), "increasing")
  expect_error(ms_model(death = -0.1), "between 0")
  expect_error(ms_model(prog = 0.1, death = 0.1, death_prog = 0.1,
                        switch_prob = 1.5), "switch_prob")
  expect_error(ms_model(time = c(0, 6), death = c(0.1, 0.2, 0.3)), "length 1 or 2")
})

test_that("ms_model sets the options of version 0.3.0", {
  m <- ms_model(prog = 0.1, death = 0.01, death_prog = 0.05)
  expect_identical(m$switch_clock, "same")
  expect_null(m$switch_af)
  expect_equal(c(m$switch_pre, m$prog_switch, m$death_switch_pre, m$dropout_pfs),
               c(0, 0.1, 0.01, 0))
  expect_identical(m$pfs_switch, "event")
  a <- ms_model(prog = 0.1, death = 0.01, death_prog = c(0.05, 0.08), time_post = c(0, 3),
                switch_rate = 0.04, switch_af = 0.7)
  expect_equal(a$switch_af, 0.7)
  expect_equal(a$death_switch, 0.7 * c(0.05, 0.08))
  p <- ms_model(time = c(0, 4), prog = 0.1, death = 0.01, death_prog = 0.05,
                switch_pre = c(0.03, 0.02), prog_switch = 0.07, death_switch_pre = 0.007,
                dropout_pfs = 0.01)
  expect_equal(p$switch_pre, c(0.03, 0.02))
  expect_equal(p$prog_switch, c(0.07, 0.07))
  expect_equal(p$death_switch_pre, c(0.007, 0.007))
  expect_equal(p$dropout_pfs, c(0.01, 0.01))
  expect_identical(ms_model(death = 0.1, switch_clock = "switch")$switch_clock, "switch")
  z <- models_s3("Z")
  expect_equal(z$death_prog, c(0, 0))
  expect_equal(z$death_switch, c(0.03, 0.035))
})

test_that("ms_model rejects invalid options of version 0.3.0", {
  expect_error(ms_model(prog = 0.1, death = 0.01, death_prog = 0.05, switch_af = 0),
               "switch_af")
  expect_error(ms_model(prog = 0.1, death = 0.01, death_prog = 0.05, switch_af = c(1, 2)),
               "switch_af")
  expect_error(ms_model(prog = 0.1, death = 0.01, death_prog = 0.05, switch_af = 0.7,
                        death_switch = 0.03), "death_switch")
  expect_error(ms_model(prog = 0.1, death = 0.01, death_prog = 0.05, switch_af = 0.7,
                        switch_clock = "switch"), "switch_clock")
  expect_error(ms_model(death = 0.01, switch_pre = 0.02, prog_switch = 0.05), "death_switch")
  expect_error(ms_model(death = 0.01, switch_pre = -1), "switch_pre")
  expect_error(ms_model(time = c(0, 4), death = 0.01, dropout_pfs = c(1, 2, 3)), "dropout_pfs")
  expect_error(ms_model(death = 0.01, prog_switch = NA), "prog_switch")
  expect_error(ms_model(death = 0.01, switch_clock = "other"), "should be one of")
})
