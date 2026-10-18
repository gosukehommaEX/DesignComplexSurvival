# Expected values: inst/validation/python/results/power_B.csv (event times by
# root finding on quadrature-based events, bounds and probabilities by exact
# recursive integration of the multivariate normal density)

test_that("ms_power matches the reference design for scenario B", {
  te <- ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2, beta = 0.1)
  res <- ms_power(scn_b(), n = 300, looks = ms_looks("os", fraction = c(0.5, 0.75, 1)),
                  test = te, events = 120)
  lk <- res$looks
  expect_equal(lk$time, c(12.068961461916697, 15.778673164508445, 21.652644979986587),
               tolerance = 1e-8)
  expect_equal(lk$events, c(60, 90, 120), tolerance = 1e-8)
  expect_equal(lk$z_mean, c(0.16368475726758816, 0.42971426938773155, 0.9535298330508851),
               tolerance = 1e-8)
  expect_lt(max(abs(lk$z_upper - c(2.9625880427275697, 2.3590177072405387,
                                   2.0140836682405143))), 1e-5)
  expect_lt(max(abs(lk$z_lower - c(-1.7648529056501607, -1.2779474607837409,
                                   2.0140836682405143))), 1e-5)
  expect_lt(max(abs(lk$prob_upper - c(0.0025638246619758, 0.024706789964818957,
                                      0.1184960368015021))), 1e-6)
  expect_lt(max(abs(lk$prob_lower - c(0.026894142136999435, 0.027600434470659064,
                                      0.799738771964045))), 1e-6)
  expect_lt(abs(res$power - 0.14576665142829687), 1e-6)
  expect_lt(abs(res$type1[["nonbinding"]] - 0.025), 1e-6)
  expect_lt(abs(res$type1[["binding"]] - 0.024999975684393556), 1e-6)
  expect_equal(lk$hr_upper, exp(-lk$z_upper / sqrt(2 / 9 * lk$events)))
})

test_that("ms_power with follow-up agrees with lrstat::lrpower()", {
  skip_if_not_installed("lrstat")
  message("lrstat is available: comparing ms_power() with lrstat::lrpower()")
  common <- list(kMax = 3, informationRates = c(0.5, 0.75, 1), alpha = 0.025,
                 typeAlphaSpending = "sfOF", allocationRatioPlanned = 2,
                 accrualTime = c(0, 3), accrualIntensity = c(10, 30),
                 piecewiseSurvivalTime = c(0, 6), lambda1 = c(0.05, 0.025),
                 lambda2 = c(0.05, 0.04), gamma1 = 0.012, gamma2 = 0.01,
                 accrualDuration = 12, followupTime = 18, fixedFollowup = FALSE,
                 typeOfComputation = "direct")
  lr <- do.call(lrstat::lrpower, common)
  res <- ms_power(scn_b(), n = 300, looks = ms_looks("os", fraction = c(0.5, 0.75, 1)),
                  test = ms_test(), followup = 18)
  expect_lt(abs(res$power - lr$overallResults$overallReject), 1e-5)
  expect_lt(max(abs(res$looks$z_upper - lr$byStageResults$efficacyBounds)), 1e-4)
  expect_lt(max(abs(res$looks$time - lr$byStageResults$analysisTime)), 1e-4)
  expect_lt(max(abs(res$looks$events - lr$byStageResults$numberOfEvents)), 1e-4)
  # beta spending: lrpower uses the type II error of the design itself, so the
  # same total beta is passed to ms_test()
  lr2 <- do.call(lrstat::lrpower, c(common, typeBetaSpending = "sfHSD",
                                    parameterBetaSpending = -2))
  pw <- lr2$overallResults$overallReject
  te <- ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2, beta = 1 - pw)
  res2 <- ms_power(scn_b(), n = 300, looks = ms_looks("os", fraction = c(0.5, 0.75, 1)),
                   test = te, followup = 18)
  expect_lt(abs(res2$power - pw), 1e-4)
  expect_lt(max(abs(res2$looks$z_lower[1:2] - lr2$byStageResults$futilityBounds[1:2])), 1e-3)
})

test_that("ms_power handles calendar-driven looks and fixed futility on the HR scale", {
  te <- ms_test(futility = "fixed", fut_bound = 1, fut_scale = "hr", binding = TRUE)
  res <- ms_power(scn_a(), n = 400, looks = ms_looks("calendar", time = c(24, 36)),
                  test = te)
  expect_equal(res$looks$time, c(24, 36))
  expect_equal(res$looks$z_lower[1], 0)
  expect_lt(res$type1[["binding"]], res$type1[["nonbinding"]])
  expect_error(ms_power(scn_a(), n = 400, looks = ms_looks("calendar", time = 24),
                        test = ms_test(), events = 100), "not used")
  expect_error(ms_power(scn_a(), n = 400, looks = ms_looks("os"), test = ms_test()),
               "exactly one")
})

test_that("ms_power with the alternative variance runs and differs little", {
  lo <- ms_looks("os", fraction = c(0.5, 1))
  r0 <- ms_power(scn_b(), n = 300, looks = lo, test = ms_test(), events = 120)
  r1 <- ms_power(scn_b(), n = 300, looks = lo, test = ms_test(variance = "alternative"),
                 events = 120)
  expect_true(all(r1$looks$z_sd > 0.9 & r1$looks$z_sd < 1.1))
  expect_lt(abs(r1$power - r0$power), 0.02)
  expect_equal(r1$looks$z_upper, r0$looks$z_upper)
})
