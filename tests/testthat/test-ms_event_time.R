# Expected values: inst/validation/python/results/event_time.csv

test_that("ms_event_time matches the reference values", {
  expect_equal(ms_event_time(scn_a(), n = 400, events = 150), 20.835972482521797,
               tolerance = 1e-8)
  tb <- ms_event_time(scn_b(), n = 300, events = c(60, 120))
  expect_equal(tb, c(12.068961461916697, 21.652644979986587), tolerance = 1e-8)
  e <- ms_events(scn_b(), n = 300, calendar_time = tb)
  expect_equal(e$events, c(60, 120), tolerance = 1e-8)
})

test_that("ms_event_time stops when the target cannot be reached", {
  expect_error(ms_event_time(scn_a(), n = 400, events = 400), "cannot reach")
  expect_error(ms_event_time(scn_a(), n = 400, events = 0), "positive")
})

test_that("ms_event_time with switching from a calendar time inverts the reference events", {
  # Expected values: calendar_looks.csv (ref_calendar.py, per patient)
  ex <- expected_csv("calendar_looks.csv")
  e <- ex[ex$design == "alt" & ex$t == 24, ]
  expect_equal(ms_event_time(scn_cal(), n = 500, events = 500 * (e$ev0 + e$ev1)), 24,
               tolerance = 1e-9)
})
