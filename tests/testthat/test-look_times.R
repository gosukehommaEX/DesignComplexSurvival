test_that(".look_times reproduces the event-driven looks of one trigger", {
  scn <- scn_b()
  acc <- .acc_setup(scn$accrual, 300)
  lt <- .look_times(scn, acc, .looks_table(ms_looks("os", fraction = c(0.5, 0.75, 1))),
                    events = 120)
  expect_identical(lt$time, .event_times(scn, acc, c(60, 90, 120), "os"))
  expect_equal(lt$target, c(60, 90, 120))
  expect_identical(lt$trigger, c("os", "os", "os"))
  expect_equal(lt$final, c(os = 120))
  lf <- .look_times(scn, acc, .looks_table(ms_looks("os", fraction = c(0.5, 1))),
                    followup = 18)
  expect_equal(lf$time[2], acc$duration + 18)
  expect_identical(lf$trigger, c("os", NA))
  ev <- .events_core(scn, acc, acc$duration + 18, "os")
  expect_equal(lf$final, c(os = ev$ev0 + ev$ev1))
  lc <- .look_times(scn, acc, .looks_table(ms_looks("calendar", time = c(12, 24))))
  expect_equal(lc$time, c(12, 24))
  expect_identical(lc$trigger, c(NA_character_, NA_character_))
  expect_equal(lc$target, c(NA_real_, NA_real_))
})

test_that(".look_times handles looks driven by the events of PFS and OS", {
  # Expected values: design2.csv (design D1, root finding on
  # ref_logrank.events(); time0 under the null scenario)
  scn <- scn_f()
  acc <- .acc_setup(scn$accrual, 400)
  lk <- .looks_table(design_def("D1")$looks)
  lt <- .look_times(scn, acc, lk, events = c(os = 220, pfs = 250))
  expect_equal(lt$time, design_expected("D1", "time")$value, tolerance = 1e-8)
  expect_equal(lt$target, c(150, 250, 220))
  expect_equal(lt$final, c(pfs = 250, os = 220))
  l0 <- .look_times(ms_null(scn), acc, lk, targets = lt$target)
  expect_equal(l0$time, design_expected("D1", "time0")$value, tolerance = 1e-8)
  # calendar-driven final look without events: final numbers at the final analysis
  lm <- .looks_table(ms_looks(c("pfs", "calendar"), fraction = c(0.5, NA), time = c(NA, 30)))
  lt2 <- .look_times(scn, acc, lm)
  ev <- .events_core(scn, acc, 30, "pfs")
  expect_equal(lt2$final, c(pfs = ev$ev0 + ev$ev1))
  expect_equal(.events_core(scn, acc, lt2$time[1], "pfs")$ev0 +
                 .events_core(scn, acc, lt2$time[1], "pfs")$ev1, 0.5 * lt2$final[[1]],
               tolerance = 1e-8)
})

test_that(".look_times checks the events and the order of the looks", {
  scn <- scn_f()
  acc <- .acc_setup(scn$accrual, 400)
  lk <- .looks_table(ms_looks(c("pfs", "os"), fraction = c(1, 1)))
  expect_error(.look_times(scn, acc, lk), "exactly one")
  expect_error(.look_times(scn, acc, lk, events = c(pfs = 380, os = 100)),
               "analysis times under the scenario are not increasing")
  expect_error(.look_times(scn, acc, .looks_table(ms_looks("calendar", time = 12)),
                           followup = 6), "not used")
  lm <- .looks_table(ms_looks(c("pfs", "calendar"), fraction = c(0.5, NA), time = c(NA, 30)))
  expect_error(.look_times(scn, acc, lm, followup = 6), "calendar-driven")
})
