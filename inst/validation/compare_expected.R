# Compare the expected values when the trial stops at the bounds
# (ms_power()$expected and the probabilities of stopping at each look under
# the design scenario, looks, and under the null scenario, looks_null;
# stage 4, Planning/stage4_plan.md, section 2.5) with rpact, lrstat,
# gsDesign, the values printed in the rpact vignette and the simulation of
# python/simulate_expected.py.
#
# 1. rpact::getPowerSurvival() (source "rpact"): the design of the rpact
#    vignette "Designing Group Sequential Trials with Two Groups and a
#    Survival Endpoint with rpact" (piecewise exponential death hazards
#    0.025, 0.04, 0.015, 0.01 and 0.007 from months 0, 6, 9, 15 and 21 in the
#    control group, hazard ratio 0.7, exponential dropout of 5 per cent in 12
#    months, 42 patients per month, n = 1000, looks at 50, 75 and 100 per
#    cent of 387 OS events, O'Brien-Fleming type efficacy bounds, mean =
#    "schoenfeld"), without futility bounds (R0) and with the fixed futility
#    bounds Z = 0 and 0.5, binding (R1) and non-binding (R2), under the
#    hazard ratio 0.7 (the design scenario) and 1 (the null scenario). rpact
#    finds its analysis times to about 2e-4 and takes the correlations of the
#    statistics from the event fractions, whereas variance = "null" takes
#    them from the information under the design scenario (stopping
#    probabilities differ by about 1.2e-5 under the hazard ratio 0.7; see
#    tests/testthat/test-ms_power.R). Variant "package": the values of the
#    package (analysis times abs:5e-4, stopping probabilities abs:5e-5,
#    expected events and patients abs:2e-3, expected study duration
#    abs:5e-4, efficacy bounds abs:1e-5). Variant "events_corr": the
#    stopping probabilities with the bounds and means of the package and the
#    correlations of the event fractions (gs_prob()), and the expected values
#    from them, which check the definitions (probabilities abs:1e-6, events
#    abs:1e-4, patients abs:2e-3 and study duration abs:5e-4, the last two
#    because of the precision of the analysis times of rpact).
# 2. The values printed in the same vignette (source "rpact_vignette";
#    rpact 4.5.0.9315, page dated 2026-08-24, retrieved 2026-10-31), judged
#    within half a unit of their last digit: getPowerSurvival() for design R0
#    (11 values) and getSampleSizeSurvival() with the hazard ratio 0.75 and
#    power 0.8 (16 values). As rpact, they take the correlations of the
#    statistics from the event fractions: the analysis times, bounds, means
#    and enrolled patients come from ms_power(), the probabilities from
#    gs_prob() with these correlations, and the number of events of
#    getSampleSizeSurvival() is solved for the power 0.8 with them
#    (uniroot()). With the correlations of the information (variance =
#    "null") the overall power of getPowerSurvival() is 0.93545 (printed
#    0.9355); the values of the package itself are compared with rpact in
#    part 1, and the solution of ms_size(solve_for = "events") is given as an
#    INFO row.
# 3. lrstat::lrpower() with typeOfComputation = "direct" (source
#    "lrstat_expected"), under the design scenario only (its null hypothesis
#    takes the final analysis at the follow-up, not at the same numbers of
#    events): scenario V2 (helpers.R; 2:1, delayed effect, piecewise accrual),
#    looks at 50, 75 and 100 per cent of the events at 18 months after the
#    end of accrual, O'Brien-Fleming type efficacy bounds, without futility
#    bounds (B0) and with the non-binding futility bounds Z = -0.5 and 0
#    (B1); and a delayed effect in piecewise exponential form (D0, D1: death
#    hazard log(2) / 12 in the control group and in the treatment group up
#    to month 4, 0.6 times it afterwards, dropout 0.005 per month, 25
#    patients per month, n = 450, looks at 40, 70 and 100 per cent of the
#    events at 12 months after the end of accrual, without futility bounds
#    and with the non-binding futility bounds Z = 0 and 0.3). Analysis times
#    and events abs:1e-4, stopping probabilities abs:1e-5, efficacy bounds
#    abs:1e-4, expected events, patients and study duration rel:1e-5 (the
#    tolerances of tests/testthat/test-ms_power.R).
# 4. gsDesign::gsProbability() (source "gsdesign_expected"): proportional
#    hazards designs of scenario V1 with looks at 150, 225 and 300 OS events
#    and mean = "schoenfeld" (P1: O'Brien-Fleming type efficacy bounds and
#    non-binding futility bounds by beta spending, Hwang-Shih-DeCani with
#    gamma -2 and beta 0.1; P2: fixed binding futility bounds Z = 0 at both
#    interim looks), with the bounds and events of the package and theta 0
#    and -log(0.7) / 2: the expected numbers of events ($en) and the stopping
#    probabilities. Under the null scenario both use the correlations of the
#    events (abs:1e-6 for the probabilities, abs:1e-4 for the events); under
#    the design scenario the package uses those of the information
#    (abs:2e-4 and abs:5e-2, above the differences of about 7e-5 and 0.01
#    that the reference implementation gives for P1), and the variant
#    "events_corr" those of the
#    events (abs:1e-6 and abs:1e-4).
# 5. Simulation (source "python_expected"; python/simulate_expected.py,
#    results/sim_python_expected.csv): designs E1-E3 defined here and,
#    independently, in Python. Judged (4se_p): the probabilities of stopping
#    early for efficacy and for futility and of stopping at each look for
#    each reason; (4se): the expected study duration, events and enrolled
#    patients; (abs:1e-6): the efficacy bounds. INFO: the mean look times,
#    events and enrolled patients at the looks, and the expected values
#    without the correlation between the look time and the decision to stop
#    (variant "indep"). The expected values of the package ignore the
#    variability of the analysis times and its correlation with the decision
#    to stop (see ms_power()); the differences that this causes are
#    registered in explained.csv.
#      E1  V1, n = 400, looks at 150, 225 and 300 OS events
#      E2  V6, n = 500, looks at 150 and 250 OS events
#      E3  V1, n = 400, looks at calendar months 24, 30 and 36
#    each with O'Brien-Fleming type efficacy bounds and the non-binding
#    futility bound Z = 0 at the interim looks.

if (!exists("val_scenarios", mode = "function")) {
  source(file.path("inst", "validation", "helpers.R"))
}

val_expected_rows <- function() {
  out <- list()
  emit <- function(source, scn, hyp, look_type, look, time, quantity, variant, criterion,
                   simulated, package, se = 0, n = NA_real_, note = "") {
    out[[length(out) + 1L]] <<- data.frame(
      id = paste(source, scn, hyp, look_type, look, ifelse(is.na(time), "", format(time)),
                 quantity, "os", "", variant, sep = "|"),
      source = source, scenario = scn, hypothesis = hyp, look_type = look_type,
      look = as.character(look), time = time, quantity = quantity, endpoint = "os", group = "",
      variant = variant, criterion = criterion, simulated = simulated, se = se, n = n,
      package = package, note = note, stringsAsFactors = FALSE)
  }
  list(emit = emit, get = function() out)
}

# expected values from the stopping probabilities (the definition of
# ms_power(): one minus the probability of stopping earlier at the last look)
val_expected_from <- function(pu, pl, time, events, enrolled) {
  k <- length(pu)
  s <- c(pu[-k] + pl[-k], 1 - sum(pu[-k] + pl[-k]))
  c(study_duration = sum(s * time), events = sum(s * events), enrolled = sum(s * enrolled))
}

val_corr_events <- function(d) {
  sqrt(outer(d, d, pmin) / outer(d, d, pmax))
}

val_rpact_scenario <- function(hr) {
  haz <- c(0.025, 0.04, 0.015, 0.01, 0.007)
  drop <- -log(0.95) / 12
  ms_scenario(ms_model(time = c(0, 6, 9, 15, 21), death = haz, dropout = drop),
              ms_model(time = c(0, 6, 9, 15, 21), death = hr * haz, dropout = drop),
              accrual = ms_accrual(rate = 42))
}

val_rpact_test <- function(futility) {
  if (futility == "R0") {
    return(ms_test(mean = "schoenfeld"))
  }
  ms_test(futility = "fixed", fut_bound = c(0, 0.5), binding = futility == "R1", mean = "schoenfeld")
}

compare_expected_rpact <- function(rw) {
  if (!requireNamespace("rpact", quietly = TRUE)) {
    message("rpact is not installed: the comparison with rpact::getPowerSurvival() is not made.")
    return(invisible(FALSE))
  }
  message("Expected values compared with rpact ", utils::packageVersion("rpact"))
  pst <- list("0 - <6" = 0.025, "6 - <9" = 0.04, "9 - <15" = 0.015, "15 - <21" = 0.01,
              ">= 21" = 0.007)
  for (des_id in c("R0", "R1", "R2")) {
    res <- ms_power(val_rpact_scenario(0.7), n = 1000,
                    looks = ms_looks("os", fraction = c(0.5, 0.75, 1)),
                    test = val_rpact_test(des_id), events = 387)
    des <- if (des_id == "R0") {
      rpact::getDesignGroupSequential(kMax = 3, alpha = 0.025, sided = 1, typeOfDesign = "asOF",
                                      informationRates = c(0.5, 0.75, 1))
    } else {
      rpact::getDesignGroupSequential(kMax = 3, alpha = 0.025, sided = 1, typeOfDesign = "asOF",
                                      informationRates = c(0.5, 0.75, 1),
                                      futilityBounds = c(0, 0.5),
                                      bindingFutility = des_id == "R1")
    }
    for (k in 1:3) {
      rw$emit("rpact", des_id, "design", "event", k, NA_real_, "bound_upper", "package",
              "abs:1e-5", as.numeric(des$criticalValues)[k], res$looks$z_upper[k])
    }
    for (hr in c(0.7, 1)) {
      pw <- rpact::getPowerSurvival(des, piecewiseSurvivalTime = pst, hazardRatio = hr,
                                    dropoutRate1 = 0.05, dropoutRate2 = 0.05, dropoutTime = 12,
                                    accrualTime = 0, accrualIntensity = 42,
                                    maxNumberOfSubjects = 1000, maxNumberOfEvents = 387,
                                    directionUpper = FALSE)
      h <- if (hr == 1) "null" else "design"
      lk <- if (hr == 1) res$looks_null else res$looks
      ex <- res$expected[res$expected$hypothesis == h, ]
      r_time <- as.numeric(pw$analysisTime)
      r_up <- as.numeric(pw$rejectPerStage)
      r_fut <- if (des_id == "R0") c(0, 0) else as.numeric(pw$futilityPerStage)
      d <- lk$events
      pr <- gs_prob(res$looks$z_lower, res$looks$z_upper, if (hr == 1) 0 else res$mean,
                    val_corr_events(d))
      ev <- val_expected_from(pr$prob_upper, pr$prob_lower, lk$time, d, lk$enrolled)
      for (k in 1:3) {
        rw$emit("rpact", des_id, h, "event", k, NA_real_, "time", "package", "abs:5e-4",
                r_time[k], lk$time[k])
        rw$emit("rpact", des_id, h, "event", k, NA_real_, "stop_efficacy", "package", "abs:5e-5",
                r_up[k], lk$prob_upper[k])
        rw$emit("rpact", des_id, h, "event", k, NA_real_, "stop_efficacy", "events_corr",
                "abs:1e-6", r_up[k], pr$prob_upper[k])
        if (k < 3 && des_id != "R0") {
          rw$emit("rpact", des_id, h, "event", k, NA_real_, "stop_futility", "package",
                  "abs:5e-5", r_fut[k], lk$prob_lower[k])
          rw$emit("rpact", des_id, h, "event", k, NA_real_, "stop_futility", "events_corr",
                  "abs:1e-6", r_fut[k], pr$prob_lower[k])
        }
      }
      r_ex <- c(events = as.numeric(pw$expectedNumberOfEvents),
                enrolled = as.numeric(pw$expectedNumberOfSubjects),
                study_duration = as.numeric(pw$studyDuration))
      tol <- c(events = "abs:2e-3", enrolled = "abs:2e-3", study_duration = "abs:5e-4")
      tol_c <- c(events = "abs:1e-4", enrolled = "abs:2e-3", study_duration = "abs:5e-4")
      for (q in names(r_ex)) {
        rw$emit("rpact", des_id, h, "overall", "", NA_real_, q, "package", tol[[q]], r_ex[[q]],
                ex[[q]])
        rw$emit("rpact", des_id, h, "overall", "", NA_real_, q, "events_corr", tol_c[[q]],
                r_ex[[q]], ev[[q]])
      }
    }
  }
  invisible(TRUE)
}

compare_expected_rpact_vignette <- function(rw) {
  src <- "rpact_vignette"
  half <- function(digits) sprintf("abs:%s", format(0.5 * 10^(-digits), scientific = FALSE))
  looks <- ms_looks("os", fraction = c(0.5, 0.75, 1))
  te <- val_rpact_test("R0")
  # probabilities and expected values with the correlations of the event
  # fractions, as in rpact
  ev_corr <- function(res) {
    d <- res$looks$events
    pu <- gs_prob(rep(-Inf, 3), res$looks$z_upper, res$mean, val_corr_events(d))$prob_upper
    list(pu = pu, ex = val_expected_from(pu, rep(0, 3), res$looks$time, d, res$looks$enrolled))
  }
  # getPowerSurvival(): hazard ratio 0.7, 387 events
  res <- ms_power(val_rpact_scenario(0.7), n = 1000, looks = looks, test = te, events = 387)
  e <- ev_corr(res)
  pr <- list(
    list("time", 1, 23.58, 2, res$looks$time[1]), list("time", 2, 34.72, 2, res$looks$time[2]),
    list("time", 3, 63.37, 2, res$looks$time[3]),
    list("stop_efficacy", 1, 0.3150, 4, e$pu[1]), list("stop_efficacy", 2, 0.4392, 4, e$pu[2]),
    list("stop_efficacy", 3, 0.1813, 4, e$pu[3]), list("power", "", 0.9355, 4, sum(e$pu)),
    list("prob_efficacy_early", "", 0.7542, 4, sum(e$pu[1:2])),
    list("events", "", 283.6, 1, e$ex[["events"]]), list("enrolled", "", 997, 0, e$ex[["enrolled"]]),
    list("study_duration", "", 38.26, 2, e$ex[["study_duration"]]))
  for (p in pr) {
    rw$emit(src, "getPowerSurvival", "design", if (nzchar(p[[2]])) "event" else "overall",
            p[[2]], NA_real_, p[[1]], "printed", half(p[[4]]), p[[3]], p[[5]])
  }
  # getSampleSizeSurvival(): hazard ratio 0.75, power 0.8, n = 1000
  scn75 <- val_rpact_scenario(0.75)
  des_d <- function(d) ms_power(scn75, n = 1000, looks = looks, test = te, events = d)
  d_sol <- stats::uniroot(function(d) sum(ev_corr(des_d(d))$pu) - 0.8, c(300, 500),
                          tol = 1e-8)$root
  sz <- des_d(d_sol)
  e <- ev_corr(sz)
  pr <- list(
    list("design", "events_at_look", 1, 193.4, 1, sz$looks$events[1]),
    list("design", "events_at_look", 2, 290.1, 1, sz$looks$events[2]),
    list("design", "events_at_look", 3, 386.8, 1, sz$looks$events[3]),
    list("design", "time", 1, 23.17, 2, sz$looks$time[1]),
    list("design", "time", 2, 33.28, 2, sz$looks$time[2]),
    list("design", "time", 3, 60.00, 2, sz$looks$time[3]),
    list("design", "accrual_duration", "", 23.81, 2, sz$accrual_duration),
    list("design", "followup", "", 36.19, 2, sz$study_duration - sz$accrual_duration),
    list("null", "stop_efficacy", 1, 0.0015, 4, sz$looks_null$prob_upper[1]),
    list("null", "stop_efficacy", 2, 0.0081, 4, sz$looks_null$prob_upper[2]),
    list("design", "stop_efficacy", 1, 0.1680, 4, e$pu[1]),
    list("design", "stop_efficacy", 2, 0.3720, 4, e$pu[2]),
    list("design", "prob_efficacy_early", "", 0.5400, 4, sum(e$pu[1:2])),
    list("design", "events", "", 318.3, 1, e$ex[["events"]]),
    list("design", "enrolled", "", 995.5, 1, e$ex[["enrolled"]]),
    list("design", "study_duration", "", 43.87, 2, e$ex[["study_duration"]]))
  for (p in pr) {
    rw$emit(src, "getSampleSizeSurvival", p[[1]], if (nzchar(p[[3]])) "event" else "overall",
            p[[3]], NA_real_, p[[2]], "printed", half(p[[5]]), p[[4]], p[[6]])
  }
  sol <- ms_size(scn75, looks = looks, test = te, power = 0.8, solve_for = "events", n = 1000,
                 round_up = FALSE)$solution
  rw$emit(src, "getSampleSizeSurvival", "design", "event", 3, NA_real_, "events_at_look",
          "ms_size", "info", 386.8, sol,
          note = "ms_size() with the correlations of the information (variance = \"null\")")
}

compare_expected_lrstat <- function(rw) {
  if (!requireNamespace("lrstat", quietly = TRUE)) {
    message("lrstat is not installed: the comparison with lrstat::lrpower() is not made.")
    return(invisible(FALSE))
  }
  message("Expected values compared with lrstat ", utils::packageVersion("lrstat"))
  defs <- val_scenarios()
  lam <- log(2) / 12
  scn_d <- ms_scenario(ms_model(death = lam, dropout = 0.005),
                       ms_model(time = c(0, 4), death = c(lam, 0.6 * lam), dropout = 0.005),
                       accrual = ms_accrual(rate = 25))
  common_b <- list(kMax = 3, informationRates = c(0.5, 0.75, 1), alpha = 0.025,
                   typeAlphaSpending = "sfOF", allocationRatioPlanned = 2,
                   accrualTime = c(0, 3), accrualIntensity = c(10, 30),
                   piecewiseSurvivalTime = c(0, 6), lambda1 = c(0.05, 0.025),
                   lambda2 = c(0.05, 0.04), gamma1 = 0.012, gamma2 = 0.01,
                   accrualDuration = 12, followupTime = 18, fixedFollowup = FALSE,
                   typeOfComputation = "direct")
  common_d <- list(kMax = 3, informationRates = c(0.4, 0.7, 1), alpha = 0.025,
                   typeAlphaSpending = "sfOF", allocationRatioPlanned = 1,
                   accrualTime = 0, accrualIntensity = 25, piecewiseSurvivalTime = c(0, 4),
                   lambda1 = c(lam, 0.6 * lam), lambda2 = c(lam, lam), gamma1 = 0.005,
                   gamma2 = 0.005, accrualDuration = 18, followupTime = 12,
                   fixedFollowup = FALSE, typeOfComputation = "direct")
  des <- list(
    B0 = list(scn = defs$V2$scenario, n = 300, frac = c(0.5, 0.75, 1), fu = 18, fut = NULL,
              args = common_b),
    B1 = list(scn = defs$V2$scenario, n = 300, frac = c(0.5, 0.75, 1), fu = 18, fut = c(-0.5, 0),
              args = common_b),
    D0 = list(scn = scn_d, n = 450, frac = c(0.4, 0.7, 1), fu = 12, fut = NULL, args = common_d),
    D1 = list(scn = scn_d, n = 450, frac = c(0.4, 0.7, 1), fu = 12, fut = c(0, 0.3),
              args = common_d))
  for (nm in names(des)) {
    d <- des[[nm]]
    te <- if (is.null(d$fut)) ms_test() else ms_test(futility = "fixed", fut_bound = d$fut)
    res <- ms_power(d$scn, n = d$n, looks = ms_looks("os", fraction = d$frac), test = te,
                    followup = d$fu)
    args <- d$args
    if (!is.null(d$fut)) {
      args$futilityBounds <- d$fut
    }
    lr <- do.call(lrstat::lrpower, args)
    bs <- lr$byStageResults
    o <- lr$overallResults
    lk <- res$looks
    ex <- res$expected[res$expected$hypothesis == "design", ]
    for (k in 1:3) {
      rw$emit("lrstat_expected", nm, "design", "event", k, NA_real_, "time", "package",
              "abs:1e-4", bs$analysisTime[k], lk$time[k])
      rw$emit("lrstat_expected", nm, "design", "event", k, NA_real_, "events_at_look", "package",
              "abs:1e-4", bs$numberOfEvents[k], lk$events[k])
      rw$emit("lrstat_expected", nm, "design", "event", k, NA_real_, "bound_upper", "package",
              "abs:1e-4", bs$efficacyBounds[k], lk$z_upper[k])
      rw$emit("lrstat_expected", nm, "design", "event", k, NA_real_, "stop_efficacy", "package",
              "abs:1e-5", bs$rejectPerStage[k], lk$prob_upper[k])
      if (k < 3 && !is.null(d$fut)) {
        rw$emit("lrstat_expected", nm, "design", "event", k, NA_real_, "stop_futility",
                "package", "abs:1e-5", bs$futilityPerStage[k], lk$prob_lower[k])
      }
    }
    rw$emit("lrstat_expected", nm, "design", "overall", "", NA_real_, "events", "package",
            "rel:1e-5", o$expectedNumberOfEvents, ex$events)
    rw$emit("lrstat_expected", nm, "design", "overall", "", NA_real_, "enrolled", "package",
            "rel:1e-5", o$expectedNumberOfSubjects, ex$enrolled)
    rw$emit("lrstat_expected", nm, "design", "overall", "", NA_real_, "study_duration",
            "package", "rel:1e-5", o$expectedStudyDuration, ex$study_duration)
  }
  invisible(TRUE)
}

compare_expected_gsdesign <- function(rw) {
  if (!requireNamespace("gsDesign", quietly = TRUE)) {
    message("gsDesign is not installed: the comparison with gsDesign::gsProbability() is not made.")
    return(invisible(FALSE))
  }
  message("Expected events compared with gsDesign ", utils::packageVersion("gsDesign"))
  v1 <- val_scenarios()$V1
  looks <- ms_looks("os", fraction = c(150, 225, 300) / 300)
  tests <- list(P1 = ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2, beta = 0.1,
                             binding = FALSE, mean = "schoenfeld"),
                P2 = ms_test(futility = "fixed", fut_bound = c(0, 0), binding = TRUE,
                             mean = "schoenfeld"))
  delta <- -log(0.7) / 2
  for (nm in names(tests)) {
    res <- ms_power(v1$scenario, n = v1$n, looks = looks, test = tests[[nm]], events = 300)
    lk <- res$looks
    g <- gsDesign::gsProbability(k = 3, theta = c(0, delta), n.I = lk$events,
                                 a = pmax(lk$z_lower, -20), b = lk$z_upper)
    for (j in 1:2) {
      h <- if (j == 1) "null" else "design"
      pk_lk <- if (j == 1) res$looks_null else lk
      ex <- res$expected[res$expected$hypothesis == h, ]
      pr <- gs_prob(lk$z_lower, lk$z_upper, if (j == 1) 0 else res$mean,
                    val_corr_events(lk$events))
      ev <- val_expected_from(pr$prob_upper, pr$prob_lower, lk$time, lk$events, lk$enrolled)
      tol_p <- if (j == 1) "abs:1e-6" else "abs:2e-4"
      tol_e <- if (j == 1) "abs:1e-4" else "abs:5e-2"
      for (k in 1:3) {
        rw$emit("gsdesign_expected", nm, h, "event", k, NA_real_, "stop_efficacy", "package",
                tol_p, g$upper$prob[k, j], pk_lk$prob_upper[k])
        if (k < 3) {
          rw$emit("gsdesign_expected", nm, h, "event", k, NA_real_, "stop_futility", "package",
                  tol_p, g$lower$prob[k, j], pk_lk$prob_lower[k])
        }
        if (j == 2) {
          rw$emit("gsdesign_expected", nm, h, "event", k, NA_real_, "stop_efficacy",
                  "events_corr", "abs:1e-6", g$upper$prob[k, j], pr$prob_upper[k])
        }
      }
      rw$emit("gsdesign_expected", nm, h, "overall", "", NA_real_, "events", "package", tol_e,
              g$en[j], ex$events)
      if (j == 2) {
        rw$emit("gsdesign_expected", nm, h, "overall", "", NA_real_, "events", "events_corr",
                "abs:1e-4", g$en[j], ev[["events"]])
      }
    }
  }
  invisible(TRUE)
}

val_expected_designs <- function() {
  defs <- val_scenarios()
  te <- function(k) ms_test(futility = "fixed", fut_bound = rep(0, k - 1L), binding = FALSE)
  list(E1 = list(def = defs$V1, looks = ms_looks("os", fraction = c(150, 225, 300) / 300),
                 events = 300, test = te(3)),
       E2 = list(def = defs$V6, looks = ms_looks("os", fraction = c(150, 250) / 250),
                 events = 250, test = te(2)),
       E3 = list(def = defs$V1, looks = ms_looks("calendar", time = c(24, 30, 36)),
                 events = NULL, test = te(3)))
}

compare_expected_sim <- function(rw) {
  p <- file.path(val_results_dir(), "sim_python_expected.csv")
  if (!file.exists(p)) {
    return(invisible(FALSE))
  }
  sims <- utils::read.csv(p, stringsAsFactors = FALSE,
                          colClasses = c(look = "character", group = "character",
                                         endpoint = "character", target_events = "character"))
  sims$look[is.na(sims$look)] <- ""
  des <- val_expected_designs()
  pw <- lapply(des, function(d) {
    ms_power(d$def$scenario, n = d$def$n, looks = d$looks, test = d$test, events = d$events)
  })
  for (i in seq_len(nrow(sims))) {
    r <- sims[i, ]
    x <- pw[[r$scenario]]
    h <- if (r$hypothesis == "alternative") "design" else "null"
    lk <- if (h == "design") x$looks else x$looks_null
    ex <- x$expected[x$expected$hypothesis == h, ]
    k <- suppressWarnings(as.integer(r$look))
    q <- r$quantity
    em <- function(variant, criterion, package, se = r$se) {
      rw$emit("python_expected", r$scenario, r$hypothesis, r$look_type, r$look, r$time, q,
              variant, criterion, r$mean, package, se = se, n = r$n)
    }
    if (q %in% c("prob_efficacy_early", "prob_futility_early")) {
      em("", "4se_p", ex[[q]])
    } else if (q %in% c("study_duration", "events", "enrolled")) {
      em("", "4se", ex[[q]])
    } else if (q %in% c("study_duration_indep", "events_indep", "enrolled_indep")) {
      em("indep", "info", ex[[sub("_indep$", "", q)]], se = NA_real_)
    } else if (q == "stop_efficacy") {
      em("", "4se_p", lk$prob_upper[k])
    } else if (q == "stop_futility") {
      em("", "4se_p", lk$prob_lower[k])
    } else if (q == "look_time") {
      em("", "info", lk$time[k])
    } else if (q == "events_at_look") {
      em("", "info", lk$events[k])
    } else if (q == "enrolled_at_look") {
      em("", "info", lk$enrolled[k])
    } else if (q == "bound") {
      em("", "abs:1e-6", x$looks$z_upper[k], se = 0)
    } else if (q == "not_reached") {
      em("", "info", 0, se = NA_real_)
    }
  }
  invisible(TRUE)
}

compare_expected <- function() {
  rw <- val_expected_rows()
  compare_expected_rpact(rw)
  compare_expected_rpact_vignette(rw)
  compare_expected_lrstat(rw)
  compare_expected_gsdesign(rw)
  compare_expected_sim(rw)
  res <- do.call(rbind, rw$get())
  rownames(res) <- NULL
  res
}
