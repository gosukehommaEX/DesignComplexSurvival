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
