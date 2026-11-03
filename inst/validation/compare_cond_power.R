# Compare the conditional power of ms_cond_power() (stage 4,
# Planning/stage4_plan.md, section 3.4) with gsDesign::gsCP().
#
# gsCP() gives, for a design of gsDesign and the statistic zi at analysis i,
# the probabilities of crossing the efficacy bound at the later analyses
# before the futility bound, with the drift theta of the canonical form
# (mean theta sqrt(n.I) and correlations sqrt(n.I_j / n.I_k) of the events).
# ms_cond_power() has the same definition: under = "trend" corresponds to
# theta = zi / sqrt(n.I[i]), "null" to theta = 0 and "hr" with the hazard
# ratio hr to theta = -log(hr) sqrt(p0 p1). Under "design" the package takes
# the correlations from the information under the design scenario
# (variance = "null"), which differ from the event fractions when the hazard
# ratio is not 1, so that these rows are INFO (differences of 2e-4 to 6e-3).
#
# 1. The design of the gsDesign vignette "Conditional power and conditional
#    error" (source "gsdesign_vignette"): gsSurv(hr = 0.7, lambdaC =
#    log(2) / 12, minfup = 24, T = 36) with toInteger(), updated by
#    gsDesign() to the first interim analysis at 117 events (117, 235 and 353
#    events, Hwang-Shih-DeCani alpha spending with gamma -4 and beta spending
#    with gamma -2, beta 0.1, non-binding), and the conditional power at
#    analysis 1 with the nominal one-sided p-value 0.04. The package design
#    has the same numbers of events and spending functions, mean =
#    "schoenfeld", 466 patients accrued uniformly over 12 months and the
#    hazard ratio exp(-2 delta) of the design of gsDesign (0.6995487 in the
#    vignette). Judged: the efficacy bounds against gsDesign (abs:1e-5;
#    alpha spending at the event fractions, identical by construction), the
#    futility bounds of the package's beta spending against gsDesign (the
#    package takes the correlations of the design scenario from the
#    information, gsDesign from the events; abs:5e-4), and, for the design
#    with the futility bounds of gsDesign (futility = "fixed"), the
#    conditional powers at analyses 2 and 3 and their sum under the three
#    assumptions of the vignette against gsCP() (abs:1e-5). The values
#    printed in the vignette (variant "printed"; retrieved 2026-10-31 from
#    https://cran.r-universe.dev/gsDesign/doc/ConditionalPowerPlot.html) are
#    judged within half a unit of their last digit: the bounds of the
#    package's design, the hazard ratio of the interim estimate and the
#    conditional powers (sums over analyses 2 and 3) of the design with the
#    futility bounds of gsDesign. The conditional powers at single analyses,
#    printed with seven or eight digits from the numerical integration of
#    gsCP(), are judged within 1e-7 rather than half a unit of the last digit
#    (the independent recomputation of
#    Planning/stage4_feasibility/feas_gsdesign_cp.py in the project folder
#    reproduced them to the printed digits, apart from two values that depend
#    on the rounding of the printed hazard ratio of the design).
# 2. A grid (source "gsdesign_cp"): three designs, both interim analyses, z
#    in -1, 0, 1 and 2 (below the efficacy bound), and under = "null",
#    "trend" and "hr" (judged, abs:1e-5), against gsCP() applied to
#    gsProbability() with the bounds and events of the package design (as in
#    tests/testthat/test-ms_cond_power.R); for the proportional hazards
#    designs G1 and G3 with mean = "schoenfeld" also "design" (INFO, against
#    gsCP() with the drift of the hazard ratio of the scenario).
#      G1  scenario V1 (helpers.R), n = 400, looks at 100, 200 and 300 OS
#          events, O'Brien-Fleming type efficacy bounds, no futility bounds,
#          mean = "schoenfeld"; "hr" with the hazard ratio 0.7
#      G2  scenario V2 (2:1, delayed effect), n = 300, looks at 60, 90 and
#          120 OS events, Pocock type efficacy bounds, fixed non-binding
#          futility bounds Z = 0 and 0.5; "hr" with the hazard ratio 0.75
#      G3  scenario V1, looks as G1, Hwang-Shih-DeCani efficacy bounds (gamma
#          -4) and non-binding futility bounds by beta spending
#          (Hwang-Shih-DeCani, gamma -2, beta 0.1), mean = "schoenfeld"; "hr"
#          with the hazard ratio 0.7

if (!exists("val_scenarios", mode = "function")) {
  source(file.path("inst", "validation", "helpers.R"))
}

compare_cond_power <- function() {
  if (!requireNamespace("gsDesign", quietly = TRUE)) {
    message("gsDesign is not installed: the comparison with gsDesign::gsCP() is not made.")
    return(NULL)
  }
  message("Conditional power compared with gsDesign ", utils::packageVersion("gsDesign"))
  out <- list()
  emit <- function(source, scn, hyp, look, quantity, variant, criterion, simulated, package,
                   note = "") {
    out[[length(out) + 1L]] <<- data.frame(
      id = paste(source, scn, hyp, "interim", look, "", quantity, "os", "", variant, sep = "|"),
      source = source, scenario = scn, hypothesis = hyp, look_type = "interim",
      look = as.character(look), time = NA_real_, quantity = quantity, endpoint = "os",
      group = "", variant = variant, criterion = criterion, simulated = simulated, se = 0,
      n = NA_real_, package = package, note = note, stringsAsFactors = FALSE)
  }
  half <- function(digits) sprintf("abs:%s", format(0.5 * 10^(-digits), scientific = FALSE))
  gcp <- function(g, i, zi, theta) {
    gsDesign::gsCP(g, theta = theta, i = i, zi = zi)$upper$prob
  }

  # 1. the design of the vignette
  if (!exists("toInteger", envir = asNamespace("gsDesign"))) {
    message("gsDesign ", utils::packageVersion("gsDesign"), " has no toInteger(): ",
            "the design of the vignette is not compared.")
  } else {
    design <- gsDesign::toInteger(gsDesign::gsSurv(hr = 0.7, lambdaC = log(2) / 12, minfup = 24,
                                                   T = 36))
    update <- gsDesign::gsDesign(k = design$k, test.type = design$test.type,
                                 alpha = design$alpha, beta = design$beta,
                                 sfu = design$upper$sf, sfupar = design$upper$param,
                                 sfl = design$lower$sf, sflpar = design$lower$param,
                                 n.I = c(117, design$n.I[2:3]),
                                 maxn.IPlan = design$n.I[design$k], delta = design$delta,
                                 delta1 = design$delta1, delta0 = design$delta0)
    hr_d <- exp(-design$delta * sqrt((1 + design$ratio)^2 / design$ratio))
    lam <- log(2) / 12
    scn <- ms_scenario(ms_model(death = lam), ms_model(death = hr_d * lam),
                       accrual = ms_accrual(time = 0, prop = 1, duration = 12))
    ev <- update$n.I
    looks <- ms_looks("os", fraction = ev / ev[3])
    beta_des <- ms_power(scn, n = 466, looks = looks,
                         test = ms_test(sf_eff = "HSD", param_eff = -4, futility = "beta",
                                        sf_fut = "HSD", param_fut = -2, beta = 0.1,
                                        binding = FALSE, mean = "schoenfeld"),
                         events = ev[3])
    fixed_des <- ms_power(scn, n = 466, looks = looks,
                          test = ms_test(sf_eff = "HSD", param_eff = -4, futility = "fixed",
                                         fut_bound = update$lower$bound[1:2], binding = FALSE,
                                         mean = "schoenfeld"),
                          events = ev[3])
    src <- "gsdesign_vignette"
    for (k in 1:3) {
      emit(src, "V117", "design", k, "bound_upper", "gsDesign", "abs:1e-5", update$upper$bound[k],
           beta_des$looks$z_upper[k])
    }
    for (k in 1:2) {
      emit(src, "V117", "design", k, "bound_lower", "gsDesign", "abs:5e-4", update$lower$bound[k],
           beta_des$looks$z_lower[k],
           note = "beta spending with the correlations of the information (package) and of the events (gsDesign)")
    }
    z1 <- stats::qnorm(1 - 0.04)
    cp <- gcp(update, 1, z1, NULL)
    th <- c(trend = z1 / sqrt(ev[1]), null = 0, hr = design$delta)
    res <- list(trend = ms_cond_power(fixed_des, 1, z1, under = "trend"),
                null = ms_cond_power(fixed_des, 1, z1, under = "null"),
                hr = ms_cond_power(fixed_des, 1, z1, under = "hr", hr = hr_d))
    for (j in seq_along(res)) {
      u <- names(res)[j]
      pk <- res[[u]]$looks$prob_upper
      for (k in 2:3) {
        emit(src, "V117", u, k, "cond_power_look", "gsCP", "abs:1e-5", cp[k - 1L, j], pk[k - 1L])
      }
      emit(src, "V117", u, "", "cond_power", "gsCP", "abs:1e-5", sum(cp[, j]),
           res[[u]]$cond_power)
    }
    emit(src, "V117", "trend", 1, "hr_interim", "gsCP", "abs:1e-10",
         exp(-th[["trend"]] * sqrt((1 + design$ratio)^2 / design$ratio)), res$trend$hr)
    # values printed in the vignette
    pr <- list(
      list("design", 1, "bound_upper", 3.0139, 4, beta_des$looks$z_upper[1]),
      list("design", 2, "bound_upper", 2.5478, 4, beta_des$looks$z_upper[2]),
      list("design", 3, "bound_upper", 1.9991, 4, beta_des$looks$z_upper[3]),
      list("design", 1, "bound_lower", -0.2451, 4, beta_des$looks$z_lower[1]),
      list("design", 2, "bound_lower", 0.9413, 4, beta_des$looks$z_lower[2]),
      list("trend", 1, "hr_interim", 0.7234658, 7, res$trend$hr),
      list("trend", 2, "cond_power_look", 0.4625266, 7, res$trend$looks$prob_upper[1]),
      list("trend", 3, "cond_power_look", 0.4364291, 7, res$trend$looks$prob_upper[2]),
      list("trend", "", "cond_power", 0.899, 3, res$trend$cond_power),
      list("null", 2, "cond_power_look", 0.03199697, 8, res$null$looks$prob_upper[1]),
      list("null", 3, "cond_power_look", 0.08886339, 8, res$null$looks$prob_upper[2]),
      list("null", "", "cond_power", 0.121, 3, res$null$cond_power),
      list("hr", 2, "cond_power_look", 0.5352688, 7, res$hr$looks$prob_upper[1]),
      list("hr", 3, "cond_power_look", 0.4017165, 7, res$hr$looks$prob_upper[2]),
      list("hr", "", "cond_power", 0.937, 3, res$hr$cond_power))
    for (p in pr) {
      crit <- if (p[[3]] == "cond_power_look") "abs:1e-7" else half(p[[5]])
      emit(src, "V117", p[[1]], p[[2]], p[[3]], "printed", crit, p[[4]], p[[6]])
    }
  }

  # 2. grid
  defs <- val_scenarios()
  fr <- c(100, 200, 300) / 300
  des <- list(
    G1 = list(def = defs$V1, looks = ms_looks("os", fraction = fr), events = 300, hr = 0.7,
              ph = TRUE, test = ms_test(mean = "schoenfeld")),
    G2 = list(def = defs$V2, looks = ms_looks("os", fraction = c(60, 90, 120) / 120), events = 120,
              hr = 0.75, ph = FALSE,
              test = ms_test(sf_eff = "Pocock", futility = "fixed", fut_bound = c(0, 0.5),
                             binding = FALSE)),
    G3 = list(def = defs$V1, looks = ms_looks("os", fraction = fr), events = 300, hr = 0.7,
              ph = TRUE,
              test = ms_test(sf_eff = "HSD", param_eff = -4, futility = "beta", sf_fut = "HSD",
                             param_fut = -2, beta = 0.1, binding = FALSE, mean = "schoenfeld"))
  )
  for (nm in names(des)) {
    d <- des[[nm]]
    x <- ms_power(d$def$scenario, n = d$def$n, looks = d$looks, test = d$test, events = d$events)
    lk <- x$looks
    p <- d$def$scenario$prob
    g <- gsDesign::gsProbability(k = 3, theta = 0, n.I = lk$events, a = pmax(lk$z_lower, -20),
                                 b = lk$z_upper)
    for (look in 1:2) {
      for (z in c(-1, 0, 1, 2)) {
        if (z >= lk$z_upper[look]) {
          next
        }
        vz <- sprintf("z=%g", z)
        theta <- c(null = 0, trend = z / sqrt(lk$events[look]),
                   hr = -log(d$hr) * sqrt(p[[1L]] * p[[2L]]))
        for (u in c("null", "trend", "hr", if (d$ph) "design")) {
          r <- if (u == "hr") {
            ms_cond_power(x, look, z, under = "hr", hr = d$hr)
          } else {
            ms_cond_power(x, look, z, under = u)
          }
          ref <- sum(gcp(g, look, z, theta[[if (u == "design") "hr" else u]]))
          emit("gsdesign_cp", nm, u, look, "cond_power", vz,
               if (u == "design") "info" else "abs:1e-5", ref, r$cond_power,
               note = if (u == "design") "correlations of the information under the design scenario" else "")
        }
      }
    }
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
