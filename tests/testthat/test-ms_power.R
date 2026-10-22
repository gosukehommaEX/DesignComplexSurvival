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

test_that("ms_power with the alternative variance and five looks matches the exact design", {
  # design V2_K5_T2_alternative of tools/miwa_accuracy.R (helper-pmvn.R);
  # version 0.2.0 (Miwa algorithm) gave a power of 0.2762 against 0.2724, a
  # beta spent at look 4 off by 7.4e-5 and a futility bound at look 4 off by
  # 1.0e-3. The futility bound at look 4 (four coordinates, Genz-Bretz
  # algorithm with abseps 1e-5) moves by about 1.4e-4 per 1e-5 of beta
  res <- ms_power(scn_b(), n = 300, looks = ms_looks("os", fraction = (1:5) / 5),
                  test = ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2,
                                 beta = 0.1, variance = "alternative"),
                  events = 150)
  x <- pmvn_v2k5
  expect_lt(max(abs(res$looks$z_upper - x$bounds_upper)), 1e-5)
  expect_lt(max(abs(res$looks$z_lower - x$bounds_lower)), 5e-4)
  expect_lt(max(abs(res$looks$prob_upper - x$bounds_prob_upper)), 5e-5)
  expect_lt(max(abs(res$looks$prob_lower - x$bounds_prob_lower)), 5e-5)
  expect_lt(abs(res$power - x$bounds_power), 5e-5)
})

test_that("ms_power with the studentized variance uses the event-driven covariance", {
  te <- ms_test(variance = "studentized")
  lk <- ms_looks("os", fraction = c(0.6, 1))
  res <- ms_power(scn_f(), n = 400, looks = lk, test = te, events = 200)
  lc <- lr_cov(scn_f(), n = 400, calendar_time = res$looks$time, endpoint = "os",
               variance = "studentized", trigger = "os")
  expect_equal(res$looks$z_sd, sqrt(diag(lc$sigma)), tolerance = 1e-10)
  r0 <- ms_power(ms_null(scn_f()), n = 400, looks = lk, test = te, events = 200)
  expect_equal(r0$looks$z_sd, c(1, 1), tolerance = 1e-8)
})

# Expected values: inst/validation/python/results/design2.csv (designs D1 to
# D4 of helper-designs.R; analysis times, moments, covariance, bounds and
# rejection probabilities computed independently in expected_design.py).
# Designs with two endpoints are run with the error tolerance 1e-6 of the
# Genz-Bretz algorithm (with_abseps()) and their probabilities are compared
# at a tolerance of 1e-5; design D4 (one endpoint) at 1e-6.

compare_design <- function(name, res, tol_prob = 1e-5) {
  lk <- res$looks
  e <- design_expected(name, "time")
  expect_equal(lk$time[match(e$look, lk$look)], e$value, tolerance = 1e-8)
  for (q in c("events", "z_mean", "z_sd")) {
    e <- design_expected(name, q)
    expect_equal(looks_value(res, e, q), e$value, tolerance = 1e-6, label = q)
  }
  for (q in c("z_upper", "z_lower")) {
    e <- design_expected(name, q)
    v <- looks_value(res, e, q)
    expect_identical(is.finite(v), is.finite(e$value), label = q)
    expect_lt(max(abs(v - e$value)[is.finite(v)]), 1e-5, label = q)
  }
  for (q in c("prob_upper", "prob_lower")) {
    e <- design_expected(name, q)
    expect_lt(max(abs(looks_value(res, e, q) - e$value)), tol_prob, label = q)
  }
}

test_that("ms_power with two endpoints matches the reference design D1", {
  res <- with_abseps(run_design("D1"))
  compare_design("D1", res)
  e <- design_expected("D1", "sigma")
  expect_equal(res$sigma[cbind(e$row, e$col)], e$value, tolerance = 1e-6)
  e <- design_expected("D1", "power")
  expect_lt(max(abs(res$power[e$endpoint] - e$value)), 1e-5)
  e <- design_expected("D1", "cum_alpha")
  expect_lt(max(abs(looks_value(res, e, "cum_alpha") - e$value)), 1e-5)
  e <- design_expected("D1", "final_events")
  expect_equal(res$final_events[e$endpoint], stats::setNames(e$value, e$endpoint))
  for (b in c("nonbinding", "binding")) {
    e <- design_expected("D1", paste0("fwer_", b))
    expect_lt(max(abs(res$fwer[[b]] - e$value)), 1e-5)
  }
  expect_identical(res$fwer$null_endpoints, c("PFS, OS", "OS"))
  expect_identical(res$looks$endpoint, c("pfs", "pfs", "os", "os"))
  expect_identical(res$looks$look, c(1L, 2L, 2L, 3L))
  expect_identical(res$looks$trigger, c("pfs", "pfs", "pfs", "os"))
  expect_equal(res$looks$events_trigger, c(150, 250, 250, 220))
  expect_equal(res$type1[["nonbinding"]], res$fwer$nonbinding[1])
})

test_that("ms_power with two endpoints matches the reference design D2", {
  res <- with_abseps(run_design("D2"))
  compare_design("D2", res)
  e <- design_expected("D2", "power")
  expect_lt(max(abs(res$power[e$endpoint] - e$value)), 1e-5)
  e <- design_expected("D2", "final_events")
  expect_equal(res$final_events[e$endpoint], stats::setNames(e$value, e$endpoint),
               tolerance = 1e-8)
  for (b in c("nonbinding", "binding")) {
    e <- design_expected("D2", paste0("fwer_", b))
    expect_lt(max(abs(res$fwer[[b]] - e$value)), 1e-5)
  }
  expect_lt(res$type1[["binding"]], res$type1[["nonbinding"]])
  # the final analysis takes place 12 months after the end of accrual
  expect_equal(res$study_duration, 18.4 + 12)
  expect_true(all(is.na(res$looks$events_trigger[res$looks$look == 1])))
})

test_that("ms_power with two endpoints matches the reference design D3", {
  res <- with_abseps(run_design("D3"))
  compare_design("D3", res)
  e <- design_expected("D3", "sigma")
  expect_equal(res$sigma[cbind(e$row, e$col)], e$value, tolerance = 1e-6)
  e <- design_expected("D3", "power")
  expect_lt(max(abs(res$power[e$endpoint] - e$value)), 1e-5)
  e <- design_expected("D3", "fwer_nonbinding")
  expect_lt(abs(res$type1[["nonbinding"]] - e$value), 1e-5)
})

test_that("ms_power with looks driven by the events of PFS and OS matches design D4", {
  res <- run_design("D4")
  compare_design("D4", res, tol_prob = 1e-6)
  e <- design_expected("D4", "power")
  expect_lt(abs(res$power - e$value), 1e-6)
  for (b in c("nonbinding", "binding")) {
    e <- design_expected("D4", paste0("type1_", b))
    expect_lt(abs(res$type1[[b]] - e$value), 1e-6)
  }
  expect_identical(res$looks$trigger, c("pfs", "os"))
  expect_equal(res$looks$events_trigger, c(250, 220))
})

test_that("ms_power with two endpoints checks its input", {
  scn <- scn_f()
  te <- ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence")
  lk <- ms_looks(c("pfs", "os"), fraction = c(1, 1))
  expect_error(ms_power(scn, 400, lk, te, events = 250), "named")
  expect_error(ms_power(scn, 400, lk, te, events = c(os = 220)), "PFS and OS")
  expect_error(ms_power(scn, 400, lk, te, events = c(pfs = 380, os = 100)), "not increasing")
  expect_error(ms_power(scn, 400, lk, te, followup = 6), "not increasing")
  expect_error(ms_power(scn, 400, lk, te, events = c(pfs = -1, os = 220)), "events")
  expect_error(ms_power(scn, 400, lk, ms_test(), events = c(pfs = 250, os = 220),
                        null_scenarios = list(os = scn_f_null("os"))), "two endpoints")
  expect_error(ms_power(scn, 400, lk, te, events = c(pfs = 250, os = 220),
                        null_scenarios = list(dfs = scn_f_null("os"))), "named by")
  s_acc <- scn_f_null("os")
  s_acc$accrual <- ms_accrual(rate = 20)
  expect_error(ms_power(scn, 400, lk, te, events = c(pfs = 250, os = 220),
                        null_scenarios = list(os = s_acc)), "allocation and accrual")
  lp <- ms_looks(c("pfs", "os"), fraction = c(1, 1), test = list("pfs", "pfs"))
  expect_error(ms_power(scn, 400, lp, te, events = c(pfs = 250, os = 220)),
               "OS is not tested")
  lo <- ms_looks("os", fraction = c(0.5, 1), test = list("os", c("pfs", "os")))
  expect_error(ms_power(scn, 400, lo, ms_test(), events = 220), "PFS")
})
