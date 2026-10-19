# Expected values: inst/validation/python/results/events.csv (adaptive
# quadrature; scenario A was also checked against the closed form for
# exponential events with uniform accrual)

test_that("ms_events matches the reference values (scenarios A and B)", {
  ea <- ms_events(scn_a(), n = 400, calendar_time = c(12, 24, 36))
  expect_equal(ea$events_control, c(32.851265591565216, 102.49878342546697, 145.65789883128582),
               tolerance = 1e-10)
  expect_equal(ea$events_treatment, c(24.470615021098048, 80.49240506214556, 121.46823450253052),
               tolerance = 1e-10)
  expect_equal(ea$dropouts_control, c(2.8436614773526436, 8.872469192704385, 12.608395698612146),
               tolerance = 1e-10)
  expect_equal(ea$enrolled, c(240, 400, 400), tolerance = 1e-12)
  eb <- ms_events(scn_b(), n = 300, calendar_time = c(15, 25))
  expect_equal(eb$events_control, c(29.82352042599617, 50.17360295400336), tolerance = 1e-10)
  expect_equal(eb$events_treatment, c(54.823684903182446, 83.3439018447882), tolerance = 1e-10)
  expect_equal(eb$dropouts_treatment, c(15.083262221574788, 27.980075821919847),
               tolerance = 1e-10)
  expect_equal(eb$events, eb$events_control + eb$events_treatment)
})

test_that("ms_events matches the reference values with progression (scenario C)", {
  ep <- ms_events(scn_c(), n = 500, calendar_time = 20, endpoint = "pfs")
  expect_equal(c(ep$events_control, ep$events_treatment),
               c(148.52838106972982, 111.9857983296598), tolerance = 1e-10)
  expect_equal(c(ep$dropouts_control, ep$dropouts_treatment),
               c(7.383078745578053, 8.926059814554732), tolerance = 1e-10)
  eo <- ms_events(scn_c(), n = 500, calendar_time = 20, endpoint = "os")
  expect_equal(c(eo$events_control, eo$events_treatment),
               c(51.30630709556631, 41.00086383313494), tolerance = 1e-9)
  expect_equal(c(eo$dropouts_control, eo$dropouts_treatment),
               c(11.383757558367165, 11.719348332575649), tolerance = 1e-9)
  ef <- ms_events(scn_c("forward"), n = 500, calendar_time = 20, endpoint = "os")
  expect_equal(c(ef$events_control, ef$events_treatment),
               c(52.79892124389747, 41.697963259701105), tolerance = 1e-9)
})

test_that("ms_events agrees with lrstat::lrstat()", {
  skip_if_not_installed("lrstat")
  lr <- lrstat::lrstat(time = c(15, 25), allocationRatioPlanned = 2,
                       accrualTime = c(0, 3), accrualIntensity = c(10, 30),
                       piecewiseSurvivalTime = c(0, 6),
                       lambda1 = c(0.05, 0.025), lambda2 = c(0.05, 0.04),
                       gamma1 = 0.012, gamma2 = 0.01, accrualDuration = 12,
                       followupTime = 100, fixedFollowup = FALSE,
                       predictTarget = 1)
  eb <- ms_events(scn_b(), n = 300, calendar_time = c(15, 25))
  expect_lt(rel_diff(eb$events_treatment, lr$nevents1), 1e-6)
  expect_lt(rel_diff(eb$events_control, lr$nevents2), 1e-6)
  expect_lt(rel_diff(eb$dropouts_treatment, lr$ndropouts1), 1e-6)
})
