# Expected values: inst/validation/python/results/model_states.csv

test_that("ms_surv matches the reference OS and PFS functions", {
  m <- models_c("reset")
  tt <- c(2, 5, 10, 20)
  os <- ms_surv(m$control, tt, "os")
  expect_equal(os$surv, c(0.9741726102902674, 0.9186583951846874, 0.8043840382493628,
                          0.5794171047106101), tolerance = 1e-10)
  expect_equal(os$dens, c(0.015467401290874, 0.021043711964349145, 0.023702415460585033,
                          0.02034214486498131), tolerance = 1e-9)
  pfs <- ms_surv(m$control, tt, "pfs")
  expect_equal(pfs$surv, c(0.8025187979624785, 0.5886049696783552, 0.3753110988513995,
                           0.1525901057568839), tolerance = 1e-12)
  expect_equal(pfs$dens, c(0.08827706777587263, 0.052974447271051966, 0.03377799889662596,
                           0.01373310951811955), tolerance = 1e-12)
  os_t <- ms_surv(m$treatment, tt, "os")
  expect_equal(os_t$surv, c(0.9800301614865685, 0.9377268277195632, 0.8477190810847506,
                            0.6482651577261405), tolerance = 1e-10)
  expect_equal(os_t$dens, c(0.01180617986481297, 0.016063434997460248, 0.019380858627952503,
                            0.01953274253444884), tolerance = 1e-9)
  mf <- models_c("forward")
  os_f <- ms_surv(mf$control, tt, "os")
  expect_equal(os_f$surv, c(0.9741726102902674, 0.9153595511221271, 0.7972824742773675,
                            0.5710672505738932), tolerance = 1e-10)
  expect_equal(os_f$dens, c(0.015467401290874, 0.022352089626070954, 0.024135313034450857,
                            0.02027050502473384), tolerance = 1e-9)
})

test_that("ms_surv returns consistent hazards and handles a single endpoint", {
  m <- models_c("reset")$control
  s <- ms_surv(m, c(3, 12), "os")
  expect_equal(s$haz, s$dens / s$surv)
  e <- ms_surv(ms_model(time = c(0, 6), death = c(0.05, 0.02)), c(3, 10), "os")
  expect_equal(e$surv, exp(-c(0.15, 0.3 + 0.08)), tolerance = 1e-14)
  expect_equal(e$haz, c(0.05, 0.02), tolerance = 1e-14)
  expect_error(ms_surv(m, -1), "nonnegative")
})

# Expected values of version 0.3.0: inst/validation/python/results/model5_states.csv

test_that("ms_surv matches the reference values with the options of version 0.3.0", {
  ex <- expected_csv("model5_states.csv")
  ex <- ex[ex$pfs_switch == "event", ]
  for (name in unique(ex$model)) {
    for (clock in c("reset", "forward")) {
      e <- ex[ex$model == name & ex$clock == clock, ]
      m <- models_s3(name, clock)
      os <- ms_surv(m, e$time, "os")
      expect_lt(max(abs(os$surv - e$os_surv)), 1e-10)
      expect_lt(max(abs(os$dens - e$os_dens)), 1e-10)
      pfs <- ms_surv(m, e$time, "pfs")
      expect_lt(max(abs(pfs$surv - e$pfs_surv)), 1e-10)
      expect_lt(max(abs(pfs$dens - e$pfs_dens)), 1e-10)
    }
  }
})

test_that("ms_surv gives PFS to progression or death whatever pfs_switch", {
  s <- scn_s3("PRC")
  tt <- c(2, 10)
  expect_equal(ms_surv(s$control, tt, "pfs"), ms_surv(models_s3("P"), tt, "pfs"))
  h <- 1e-5
  a <- ms_surv(models_s3("P2", "forward"), c(2.5, 12.5) + h, "os")$surv
  b <- ms_surv(models_s3("P2", "forward"), c(2.5, 12.5) - h, "os")$surv
  expect_equal(ms_surv(models_s3("P2", "forward"), c(2.5, 12.5), "os")$dens, (b - a) / (2 * h),
               tolerance = 1e-7)
})
