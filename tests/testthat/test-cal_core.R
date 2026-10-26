# Expected values: inst/validation/python/results/calendar_looks.csv (per patient, n = 500)

test_that(".cal_core matches the reference and gives zeros before the start", {
  ex <- expected_csv("calendar_looks.csv")
  e <- ex[ex$design == "alt", ]
  scn <- scn_cal()
  acc <- .acc_setup(scn$accrual, 500)
  core <- .cal_core(scn, acc, c(0, 24, 36))
  expect_named(core, c("time", "events_control", "events_treatment", "events",
                       "dropouts_control", "dropouts_treatment", "score_mean", "info", "z_mean"))
  expect_equal(core$events_control, c(0, 500 * e$ev0), tolerance = 1e-10)
  expect_equal(core$events_treatment, c(0, 500 * e$ev1), tolerance = 1e-10)
  expect_equal(core$dropouts_control, c(0, 500 * e$do0), tolerance = 1e-10)
  expect_equal(core$events, core$events_control + core$events_treatment)
  expect_equal(core$info[-1L], 500 * e$v, tolerance = 1e-10)
  expect_equal(core$z_mean[-1L], -sqrt(500) * e$mu / sqrt(e$v), tolerance = 1e-9)
})

test_that(".cal_core agrees with .lr_core when nobody or everybody can switch", {
  acc <- .acc_setup(scn_cal()$accrual, 500)
  m <- models_c("forward")
  m0 <- m$control
  m0$switch_prob[] <- 0
  m0$switch_rate[] <- 0
  no_sw <- ms_scenario(m0, m$treatment, accrual = scn_cal()$accrual)
  # switching from month 40: nobody switches by month 36
  expect_equal(.cal_core(scn_cal(tc = 40), acc, c(20, 36)),
               .lr_core(no_sw, acc, c(20, 36), "os", "direct"), tolerance = 1e-9)
  # switching from month 1e-8: everybody can switch from randomization
  expect_equal(.cal_core(scn_cal(tc = 1e-8), acc, c(20, 36)),
               .lr_core(scn_cal(tc = 0), acc, c(20, 36), "os", "direct"), tolerance = 1e-9)
})
