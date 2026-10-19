# Expected values: inst/validation/python/results/moments.csv and
# alt_cov.csv (variance under the alternative from the formulas of Luo et al.
# 2019 as printed, by nested adaptive quadrature)

test_that("lr_moments matches the reference values without progression", {
  ma <- lr_moments(scn_a(), n = 400, calendar_time = c(24, 36))
  expect_equal(ma$score_mean, c(-16.22651665423436, -23.655807223764107), tolerance = 1e-10)
  expect_equal(ma$info, c(45.511550039021955, 66.00320276127533), tolerance = 1e-10)
  expect_equal(ma$z_mean, c(2.40527360996092, 2.911757933368407), tolerance = 1e-10)
  sa <- lr_moments(scn_a(), n = 400, calendar_time = c(24, 36), mean = "schoenfeld")
  expect_equal(sa$z_mean, c(2.41244652203202, 2.9147465361716014), tolerance = 1e-10)
  # constant hazard ratio: the Schoenfeld mean is -log(HR) sqrt(p0 p1 D)
  expect_equal(sa$z_mean, -log(0.7) * sqrt(0.25 * sa$events), tolerance = 1e-12)
  mb <- lr_moments(scn_b(), n = 300, calendar_time = c(15, 25))
  expect_equal(mb$score_mean, c(-1.5842658069726598, -6.37016594222324), tolerance = 1e-10)
  expect_equal(mb$info, c(18.817550458247293, 29.418218596536327), tolerance = 1e-10)
  expect_equal(mb$z_mean, c(0.3652132322828265, 1.174471710657412), tolerance = 1e-10)
  sb <- lr_moments(scn_b(), n = 300, calendar_time = c(15, 25), mean = "schoenfeld")
  expect_equal(sb$z_mean, c(0.344886671925964, 1.1163384801879728), tolerance = 1e-10)
  expect_equal(mb$z_var, c(1, 1))
})

test_that("lr_moments matches the reference values with progression", {
  mp <- lr_moments(scn_c(), n = 500, calendar_time = c(10, 20), endpoint = "pfs")
  expect_equal(mp$score_mean, c(-9.523426495379917, -29.79313778086134), tolerance = 1e-10)
  expect_equal(mp$info, c(20.337952807172893, 64.31999243326673), tolerance = 1e-10)
  sp <- lr_moments(scn_c(), n = 500, calendar_time = c(10, 20), endpoint = "pfs",
                   mean = "schoenfeld")
  expect_equal(sp$z_mean, c(2.1295812836868877, 3.7306089255900674), tolerance = 1e-10)
  mo <- lr_moments(scn_c(), n = 500, calendar_time = 20, endpoint = "os")
  expect_equal(mo$score_mean, -5.997014230007823, tolerance = 1e-9)
  expect_equal(mo$info, 23.064811036964382, tolerance = 1e-9)
  so <- lr_moments(scn_c(), n = 500, calendar_time = 20, endpoint = "os",
                   mean = "schoenfeld")
  expect_equal(so$z_mean, 1.2536530740621616, tolerance = 1e-9)
  mf <- lr_moments(scn_c("forward"), n = 500, calendar_time = 20, endpoint = "os")
  expect_equal(mf$score_mean, -6.486544839804843, tolerance = 1e-9)
  expect_equal(mf$info, 23.60987009337309, tolerance = 1e-9)
})

test_that("lr_moments with the alternative variance matches Luo et al. (2019)", {
  vb <- lr_moments(scn_b(), n = 300, calendar_time = c(15, 25), variance = "alternative")
  expect_equal(vb$z_var, c(0.06317798711200127 / 0.06272516819415765,
                            0.09823341210532265 / 0.09806072865512112), tolerance = 1e-8)
  vc <- lr_moments(scn_c(), n = 500, calendar_time = 20, endpoint = "pfs",
                   variance = "alternative")
  expect_equal(vc$z_var, 0.12554210853294953 / 0.12863998486653347, tolerance = 1e-8)
})

test_that("lr_moments agrees with lrstat::lrstat()", {
  skip_if_not_installed("lrstat")
  lr <- lrstat::lrstat(time = c(15, 25), allocationRatioPlanned = 2,
                       accrualTime = c(0, 3), accrualIntensity = c(10, 30),
                       piecewiseSurvivalTime = c(0, 6),
                       lambda1 = c(0.05, 0.025), lambda2 = c(0.05, 0.04),
                       gamma1 = 0.012, gamma2 = 0.01, accrualDuration = 12,
                       followupTime = 100, fixedFollowup = FALSE,
                       predictTarget = 2)
  mb <- lr_moments(scn_b(), n = 300, calendar_time = c(15, 25))
  expect_lt(rel_diff(mb$score_mean, lr$uscore), 1e-5)
  expect_lt(rel_diff(mb$info, lr$vscore), 1e-5)
})
