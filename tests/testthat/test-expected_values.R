test_that(".expected_values sums the values over the looks", {
  # stopping probabilities 0.15, 0.3 and 0.55 at the three looks
  ex <- .expected_values(c(0.1, 0.2, 0.3), c(0.05, 0.1, 0.25), time = c(10, 20, 30),
                         events = c(50, 100, 150), enrolled = c(200, 300, 300))
  expect_equal(ex$prob_efficacy_early, 0.3)
  expect_equal(ex$prob_futility_early, 0.15)
  expect_equal(ex$study_duration, 0.15 * 10 + 0.3 * 20 + 0.55 * 30)
  expect_equal(ex$events, 0.15 * 50 + 0.3 * 100 + 0.55 * 150)
  expect_equal(ex$enrolled, 0.15 * 200 + 0.3 * 300 + 0.55 * 300)
})

test_that(".expected_values gives the values of the look of a design with one look", {
  ex <- .expected_values(0.8, 0.2, time = 36, events = 300, enrolled = 400)
  expect_equal(unlist(ex), c(prob_efficacy_early = 0, prob_futility_early = 0,
                             study_duration = 36, events = 300, enrolled = 400))
})

test_that(".expected_values takes the stopping probability of the last look by difference", {
  a <- .expected_values(c(0.2, 0.5), c(0.1, 0.19), c(1, 2), c(10, 20), c(5, 5))
  b <- .expected_values(c(0.2, 0.6), c(0.1, 0.1), c(1, 2), c(10, 20), c(5, 5))
  expect_equal(a, b)
  expect_equal(a$study_duration, 0.3 * 1 + 0.7 * 2)
  expect_equal(a$enrolled, 5)
})
