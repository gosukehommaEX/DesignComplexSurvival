# Expected values: inst/validation/python/results/cond_power_stage4.csv (expected_stage4.py):
# conditional moments from the independent increments or from the precision matrix, probabilities
# of the later looks by recursive conditioning with adaptive quadrature. The bounds of the package
# agree with those of the reference within 1e-5 (test-ms_power.R), which moves the probabilities by
# less than 1e-5; with the bounds of the reference (stage4_cp(, "z_upper")) the probabilities are
# compared at 1e-6.

cp_design_b <- function() {
  te <- ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2, beta = 0.1)
  ms_power(scn_b(), n = 300, looks = ms_looks("os", fraction = c(0.5, 0.75, 1)), test = te,
           events = 120)
}

# design G of cond_power_stage4.csv: proportional hazards, mean = "schoenfeld", futility bounds by
# beta spending that are high at look 2
cp_design_g <- function() {
  te <- ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2, beta = 0.1, mean = "schoenfeld")
  ms_power(scn_a(), n = 400, looks = ms_looks("os", fraction = c(1, 2, 3) / 3), test = te,
           events = 300)
}

# scenario B with the death hazards 0.05 and 0.03 in the treatment group
scn_b_trt3 <- function() {
  ms_scenario(ms_model(time = c(0, 6), death = c(0.05, 0.04), dropout = 0.01),
              ms_model(time = c(0, 6), death = c(0.05, 0.03), dropout = 0.012),
              allocation = 2, accrual = ms_accrual(time = c(0, 3), rate = c(10, 30)))
}

compare_cp <- function(res, case, tol_prob = 1e-5) {
  lk <- res$looks
  expect_lt(abs(res$cond_power - stage4_cp(case, "cond_power")), tol_prob, label = case)
  expect_lt(abs(res$cond_power_nonbinding - stage4_cp(case, "cond_power_nonbinding")), tol_prob,
            label = case)
  expect_equal(lk$cond_mean, stage4_cp(case, "cond_mean"), tolerance = 1e-7, label = case)
  expect_equal(lk$cond_sd, stage4_cp(case, "cond_sd"), tolerance = 1e-7, label = case)
  expect_lt(max(abs(lk$prob_upper - stage4_cp(case, "prob_upper"))), tol_prob, label = case)
  expect_lt(max(abs(lk$prob_lower - stage4_cp(case, "prob_lower"))), tol_prob, label = case)
  expect_equal(lk$time, stage4_cp(case, "time"), tolerance = 1e-8, label = case)
  expect_equal(lk$events, stage4_cp(case, "events"), tolerance = 1e-8, label = case)
  expect_equal(lk$cum_upper, cumsum(lk$prob_upper), label = case)
}

# the probabilities with the bounds of the reference design (case_design)
compare_cp_ref_bounds <- function(x, res, case, case_design, tol_prob = 1e-6) {
  mo <- .cp_moments(x, res$look, res$z[length(res$z)], res$under,
                    hr = if (res$under == "hr") res$hr)
  up <- stage4_cp(case_design, "z_upper")
  lo <- if (x$settings$futility == "none") {
    c(rep(-Inf, length(up) - 1L), up[length(up)])
  } else {
    stage4_cp(case_design, "z_lower")
  }
  r <- .cond_probs(mo$mean, mo$sigma, res$look, res$z, lo, up)
  expect_lt(abs(r$cond_power - stage4_cp(case, "cond_power")), tol_prob, label = case)
  expect_lt(abs(r$cond_power_nonbinding - stage4_cp(case, "cond_power_nonbinding")), tol_prob,
            label = case)
  expect_lt(max(abs(r$looks$prob_lower - stage4_cp(case, "prob_lower"))), tol_prob, label = case)
}

test_that("ms_cond_power matches the reference for design B", {
  x <- cp_design_b()
  for (look in 1:2) {
    z <- c(0.5, 1)[look]
    for (u in c("design", "null", "trend", "hr")) {
      case <- sprintf("B_l%d_%s", look, u)
      res <- ms_cond_power(x, look, z, under = u, hr = if (u == "hr") 0.75)
      compare_cp(res, case)
      compare_cp_ref_bounds(x, res, case, "B_design")
      expect_identical(res$looks$z_upper, x$looks$z_upper[(look + 1):3])
      expect_identical(res$looks$z_lower, x$looks$z_lower[(look + 1):3])
      if (u %in% c("trend", "hr")) {
        expect_equal(res$hr, stage4_cp(case, "hr"), tolerance = 1e-10, label = case)
      } else {
        expect_identical(res$hr, NA_real_)
      }
    }
    res <- ms_cond_power(x, look, z, under = "scenario", scenario = scn_b_trt3())
    compare_cp(res, sprintf("B_l%d_scenario", look))
  }
  # below the non-binding futility bound of look 2
  res <- ms_cond_power(x, 2, -1.5)
  expect_lt(-1.5, x$looks$z_lower[2])
  compare_cp(res, "B_l2low_design")
  expect_identical(res[c("look", "z", "under")], list(look = 2L, z = -1.5, under = "design"))
})

test_that("ms_cond_power drops the futility bounds of the later interim looks", {
  x <- cp_design_g()
  res <- ms_cond_power(x, 1, 0.5)
  lk <- x$looks
  r <- .cond_probs(x$mean, x$sigma, 1L, 0.5, c(lk$z_lower[1], -Inf, lk$z_lower[3]), lk$z_upper)
  expect_equal(res$cond_power_nonbinding, r$cond_power)
  expect_gt(res$cond_power_nonbinding - res$cond_power, 0.01)
  expect_equal(res$cond_power, sum(res$looks$prob_upper))
  compare_cp(res, "G_l1_design")
  # with one later look the two coincide
  r2 <- ms_cond_power(x, 2, 1.5)
  expect_identical(r2$cond_power_nonbinding, r2$cond_power)
  expect_equal(r2$looks$prob_upper + r2$looks$prob_lower, 1)
})

test_that("ms_cond_power given the earlier statistics with independent increments", {
  x <- cp_design_b()
  for (u in c("design", "null", "trend", "hr")) {
    h <- ms_cond_power(x, 2, c(0.5, 1), under = u, hr = if (u == "hr") 0.75)
    cu <- ms_cond_power(x, 2, 1, under = u, hr = if (u == "hr") 0.75)
    expect_equal(h$cond_power, cu$cond_power, tolerance = 1e-10, label = u)
    expect_equal(h$looks$cond_mean, cu$looks$cond_mean, tolerance = 1e-10, label = u)
    expect_equal(h$looks$cond_sd, cu$looks$cond_sd, tolerance = 1e-10, label = u)
    expect_identical(h$z, c(0.5, 1))
  }
})

test_that("ms_cond_power under a scenario gives the design and the null hypothesis", {
  x <- cp_design_b()
  for (look in 1:2) {
    z <- c(0.5, 1)[look]
    a <- ms_cond_power(x, look, z, under = "scenario", scenario = x$scenario)
    b <- ms_cond_power(x, look, z)
    expect_equal(a$cond_power, b$cond_power, tolerance = 1e-10)
    expect_equal(a$looks, b$looks, tolerance = 1e-10)
    a <- ms_cond_power(x, look, z, under = "scenario", scenario = ms_null(x$scenario))
    b <- ms_cond_power(x, look, z, under = "null")
    expect_equal(a$cond_power, b$cond_power, tolerance = 1e-10)
    expect_equal(a$looks, b$looks, tolerance = 1e-10)
  }
})

test_that("ms_cond_power with the studentized variance uses the earlier statistics", {
  x <- ms_power(scn_b(), n = 300, looks = ms_looks("os", fraction = c(1, 2, 3) / 3),
                test = ms_test(variance = "studentized"), events = 120)
  lk <- x$looks
  expect_equal(lk$time, stage4_cp("S_design", "time"), tolerance = 1e-8)
  expect_equal(lk$z_mean, stage4_cp("S_design", "z_mean"), tolerance = 1e-7)
  expect_equal(x$sigma, stage4_cp("S_design", "sigma"), tolerance = 1e-6)
  expect_lt(max(abs(lk$z_upper - stage4_cp("S_design", "z_upper"))), 1e-5)
  expect_false(.is_markov(x$sigma))
  cases <- list(list("S_l1", 1, 0.8), list("S_l2_hist", 2, c(0.8, 1.2)), list("S_l2_cur", 2, 1.2))
  out <- list()
  for (cs in cases) {
    res <- ms_cond_power(x, cs[[2]], cs[[3]])
    compare_cp(res, cs[[1]])
    compare_cp_ref_bounds(x, res, cs[[1]], "S_design")
    out[[cs[[1]]]] <- res$cond_power
  }
  expect_gt(abs(out$S_l2_hist - out$S_l2_cur), 1e-4)
})

test_that("ms_cond_power of a design with two looks has a closed form", {
  x <- ms_power(scn_a(), n = 400, looks = ms_looks("os", fraction = c(0.5, 1)), test = ms_test(),
                events = 300)
  u2 <- x$looks$z_upper[2]
  z <- 1.2
  s <- x$sigma
  m2 <- x$mean[2] + s[1, 2] / s[1, 1] * (z - x$mean[1])
  s2 <- sqrt(s[2, 2] - s[1, 2]^2 / s[1, 1])
  res <- ms_cond_power(x, 1, z)
  expect_equal(res$cond_power, pnorm((m2 - u2) / s2), tolerance = 1e-12)
  expect_equal(res$looks$prob_lower, pnorm((u2 - m2) / s2), tolerance = 1e-12)
  expect_identical(res$cond_power_nonbinding, res$cond_power)
  d <- x$looks$events
  res <- ms_cond_power(x, 1, z, under = "trend")
  expect_equal(res$cond_power, pnorm((z * sqrt(d[2] / d[1]) - u2) / sqrt(1 - d[1] / d[2])),
               tolerance = 1e-12)
  res <- ms_cond_power(x, 1, z, under = "hr", hr = 0.8)
  m2 <- -log(0.8) * sqrt(d[2] / 4) + sqrt(d[1] / d[2]) * (z + log(0.8) * sqrt(d[1] / 4))
  expect_equal(res$cond_power, pnorm((m2 - u2) / sqrt(1 - d[1] / d[2])), tolerance = 1e-12)
  res <- ms_cond_power(x, 1, z, under = "null")
  s <- x$sigma_null
  expect_equal(res$cond_power, pnorm((s[1, 2] * z - u2) / sqrt(1 - s[1, 2]^2)), tolerance = 1e-12)
})

test_that("ms_cond_power agrees with gsDesign::gsCP()", {
  skip_if_not_installed("gsDesign")
  message("Conditional power compared with gsDesign ", utils::packageVersion("gsDesign"))
  # gsCP() takes the correlations from the events, and "trend", "hr" and "null" use the same
  # correlations (under the null hypothesis with the same dropout in both groups the information
  # is proportional to the events). Under "design" the package takes them from the information
  # under the design scenario, which changes the conditional power by about 1e-3: both are compared
  # with the reference (G_*_design and G_*_events)
  x <- cp_design_g()
  lk <- x$looks
  expect_equal(lk$time, stage4_cp("G_design", "time"), tolerance = 1e-8)
  expect_equal(lk$z_mean, stage4_cp("G_design", "z_mean"), tolerance = 1e-8)
  expect_lt(max(abs(lk$z_upper - stage4_cp("G_design", "z_upper"))), 1e-5)
  expect_lt(max(abs(lk$z_lower - stage4_cp("G_design", "z_lower"))), 1e-5)
  g <- gsDesign::gsProbability(k = 3, theta = 0, n.I = lk$events, a = pmax(lk$z_lower, -20),
                               b = lk$z_upper)
  delta <- -log(0.7) / 2
  for (look in 1:2) {
    z <- c(0.5, 1.5)[look]
    gcp <- function(theta) {
      sum(gsDesign::gsCP(g, theta = theta, i = look, zi = z)$upper$prob)
    }
    lab <- paste("look", look)
    expect_lt(abs(ms_cond_power(x, look, z, under = "null")$cond_power - gcp(0)), 1e-5,
              label = lab)
    expect_lt(abs(ms_cond_power(x, look, z, under = "trend")$cond_power -
                    gcp(z / sqrt(lk$events[look]))), 1e-5, label = lab)
    expect_lt(abs(ms_cond_power(x, look, z, under = "hr", hr = 0.7)$cond_power - gcp(delta)), 1e-5,
              label = lab)
    expect_lt(abs(gcp(delta) - stage4_cp(sprintf("G_l%d_events", look), "cond_power")), 1e-5,
              label = lab)
    compare_cp(ms_cond_power(x, look, z), sprintf("G_l%d_design", look))
  }
})

test_that("ms_cond_power checks its input", {
  x <- cp_design_b()
  expect_error(ms_cond_power(list(), 1, 0), "'ms_power' object")
  y <- x
  y$settings$endpoint <- c("pfs", "os")
  expect_error(ms_cond_power(y, 1, 0), "one endpoint only")
  y <- x
  y$sigma <- NULL
  expect_error(ms_cond_power(y, 1, 0), "compute it again")
  x1 <- ms_power(scn_a(), n = 400, looks = ms_looks("os", fraction = 1), test = ms_test(),
                 events = 300)
  expect_error(ms_cond_power(x1, 1, 0), "no interim look")
  for (lk in list(0, 3, 1.5, NA, c(1, 2), "1")) {
    expect_error(ms_cond_power(x, lk, 0), "an integer from 1 to 2")
  }
  for (z in list(c(0.1, 0.2, 0.3), NA, Inf, "0", numeric(0))) {
    expect_error(ms_cond_power(x, 2, z), "statistics at looks 1 to 2")
  }
  expect_error(ms_cond_power(x, 1, x$looks$z_upper[1]), "crosses the efficacy bound")
  expect_error(ms_cond_power(x, 2, c(3, 1)), "crosses the efficacy bound")
  expect_error(ms_cond_power(x, 1, 0, under = "hr"), "positive number")
  expect_error(ms_cond_power(x, 1, 0, under = "hr", hr = -1), "positive number")
  expect_error(ms_cond_power(x, 1, 0, hr = 0.7), "only with under = \"hr\"")
  expect_error(ms_cond_power(x, 1, 0, under = "scenario"), "'ms_scenario' object")
  expect_error(ms_cond_power(x, 1, 0, scenario = scn_b()), "only with under = \"scenario\"")
  expect_error(ms_cond_power(x, 1, 0, under = "scenario", scenario = scn_a()),
               "allocation and accrual")
  expect_error(ms_cond_power(x, 1, 0, under = "other"), "should be one of")
  # binding futility bound at look 1 (Z = 0)
  te <- ms_test(futility = "fixed", fut_bound = 1, fut_scale = "hr", binding = TRUE)
  xb <- ms_power(scn_a(), n = 400, looks = ms_looks("calendar", time = c(24, 36)), test = te)
  expect_error(ms_cond_power(xb, 1, -0.1), "binding futility bound")
  expect_gt(ms_cond_power(xb, 1, 0.1)$cond_power, 0)
})
