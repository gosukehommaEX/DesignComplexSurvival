# The expected events change sign around the time found (steps of 1e-8, about 8e-8 events)
# and are equal to the target at that time
check_event_time <- function(scn, targets) {
  acc <- .acc_setup(scn$accrual, 500)
  got <- .cal_event_times(scn, acc, targets)
  for (k in seq_along(got)) {
    ev <- .cal_core(scn, acc, got[k] + c(-1e-8, 0, 1e-8))$events
    expect_lt(ev[1L], targets[k])
    expect_gt(ev[3L], targets[k])
    expect_lt(abs(ev[2L] - targets[k]), 1e-9)
  }
}

test_that(".cal_event_times finds the times of the events (forward clock)", {
  check_event_time(scn_cal(), c(60, 230))
})

test_that(".cal_event_times finds the times of the events (reset clock)", {
  check_event_time(scn_cal(clock = "reset", tc = 13.7), 230)
})

test_that(".cal_event_times gives NA without a starting time", {
  acc <- .acc_setup(scn_cal()$accrual, 500)
  # more OS events than patients: no starting time from the scenario with switch_start = 0
  got <- .cal_event_times(scn_cal(), acc, c(100, 600))
  expect_true(got[1L] > 0)
  expect_identical(got[2L], NA_real_)
})
