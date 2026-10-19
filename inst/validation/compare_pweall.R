# Comparison with PWEALL (Luo et al. 2019) (Phase 3 plan, sections 4.3 and 12).
#
# Cases: Luo et al. (2019) Examples 1-3 and Barthel et al. (2006) Table IV,
# rows 5 and 12. Quantities at fixed calendar times: group-wise expected
# events (pwecxpwu()), the mean of the log-rank score, its information and
# its variance under the alternative (overallvar() with beta = 0, per patient,
# times n), and the power with the variance under the alternative and under
# the null hypothesis (pwepower(), columns 2 and 3). PWEALL integrates with a
# tolerance eps; eps = veps = 0.001 is used.
#
# Criteria: power within 0.002 (abs:0.002); the other quantities within a
# relative difference of 0.2% (rel:0.002).
#
# PWEALL parameterization: all hazards share the breakpoints tchange;
# type = 1 is the Markov crossover (clock = "forward"), type = 3 the hybrid
# crossover whose hazard after crossover is rp2 rate2(t - u) + (1 - rp2)
# rate4(t). The accrual is given by densities u on intervals ending at ut.

if (!exists("val_scenarios", mode = "function")) {
  source(file.path("inst", "validation", "helpers.R"))
}

compare_pweall <- function(eps = 0.001) {
  if (!requireNamespace("PWEALL", quietly = TRUE)) {
    return(NULL)
  }
  cases <- list()
  zero <- function(k) rep(0, k)
  # Luo et al. (2019), Example 1 (months)
  tm <- c(0, 6, 12, 24)
  lam0 <- c(0.09, 0.09, 0.07, 0.05) / 12
  drop <- c(0, 0, 0.015, 0.015) / 12
  disc <- c(0.04, 0.04, 0.03, 0.03) / 12
  hr_d <- c(1, 0.75, 0.75, 0.75)
  acc1 <- list(taur = 14, ut = 1:14, u = c(100, 200, 300, 300, rep(400, 8), 300, 200) / 4600)
  ctl1 <- ms_model(time = tm, death = lam0, dropout = drop)
  pctl1 <- list(rate10 = lam0, rate20 = lam0, rate30 = zero(4), ratec0 = drop, type0 = 1)
  a1 <- ms_accrual(time = 0:13, rate = c(100, 200, 300, 300, rep(400, 8), 300, 200))
  for (v in c("constHR", "delay-w1", "delay-w0")) {
    trt <- switch(v,
      "constHR" = ms_model(time = tm, death = 0.8 * lam0, dropout = drop),
      "delay-w1" = ms_model(time = tm, death = hr_d * lam0, dropout = drop),
      "delay-w0" = ms_model(time = tm, death = hr_d * lam0, dropout = drop, prog = disc,
                            death_prog = lam0, time_post = tm, clock = "forward"))
    ptrt <- switch(v,
      "constHR" = list(rate11 = 0.8 * lam0, rate21 = 0.8 * lam0, rate31 = zero(4)),
      "delay-w1" = list(rate11 = hr_d * lam0, rate21 = hr_d * lam0, rate31 = zero(4)),
      "delay-w0" = list(rate11 = hr_d * lam0, rate21 = lam0, rate31 = disc))
    cases[[paste0("Luo2019-Ex1-", v)]] <- list(
      scenario = ms_scenario(ctl1, trt, accrual = a1), n = 4600, pi1 = 0.5, alpha = 0.025,
      times = c(30, 38.9, 48), tchange = tm, acc = acc1,
      pw = c(pctl1, ptrt, list(ratec1 = drop, type1 = 1)))
  }
  # Example 2 (2:1, one-sided 0.5%)
  tm2 <- c(0, 1, 2, 3, 4, 5, 6, 12)
  lam2 <- c(0.9, 0.75, 0.6, 0.55, 0.5, 0.4, 0.3, 0.1) / 12
  d2 <- rep(0.03 / 12, 8)
  acc2 <- list(taur = 15, ut = c(3, 6, 9, 12, 15), u = c(10, 30, 60, 100, 100) / 900)
  a2 <- ms_accrual(time = c(0, 3, 6, 9, 12), rate = c(10, 30, 60, 100, 100))
  ctl2 <- ms_model(time = tm2, death = lam2, dropout = 0.03 / 12)
  for (v in c("nodil", "dil")) {
    trt <- if (v == "nodil") {
      ms_model(time = tm2, death = 0.6 * lam2, dropout = 0.03 / 12)
    } else {
      ms_model(time = tm2, death = 0.6 * lam2, dropout = 0.03 / 12, prog = 0.2 / 12,
               death_prog = lam2, time_post = tm2, clock = "forward")
    }
    ptrt <- if (v == "nodil") {
      list(rate11 = 0.6 * lam2, rate21 = 0.6 * lam2, rate31 = zero(8))
    } else {
      list(rate11 = 0.6 * lam2, rate21 = lam2, rate31 = rep(0.2 / 12, 8))
    }
    cases[[paste0("Luo2019-Ex2-", v)]] <- list(
      scenario = ms_scenario(ctl2, trt, allocation = 2, accrual = a2), n = 900, pi1 = 2 / 3,
      alpha = 0.005, times = c(18, 24, 34.5), tchange = tm2, acc = acc2,
      pw = c(list(rate10 = lam2, rate20 = lam2, rate30 = zero(8), ratec0 = d2, type0 = 1),
             ptrt, list(ratec1 = d2, type1 = 1)))
  }
  # Example 3: hybrid crossover, hazard after crossover 0.7 times the
  # experimental hazard (from crossover) plus 0.3 times the control hazard
  l0 <- 1.5 / 12
  dp <- c(0.7 * 0.8 * l0 + 0.3 * l0, 0.7 * 0.5 * l0 + 0.3 * l0)
  acc3 <- list(taur = 15, ut = c(3, 6, 9, 12, 15), u = c(10, 20, 30, 30, 30) / 360)
  cases[["Luo2019-Ex3"]] <- list(
    scenario = ms_scenario(
      ms_model(death = l0, dropout = 0.3 / 12, prog = 0.8 / 12, death_prog = dp,
               time_post = c(0, 3), clock = "reset"),
      ms_model(time = c(0, 3), death = c(0.8, 0.5) * l0, dropout = 0.3 / 12),
      accrual = ms_accrual(time = c(0, 3, 6, 9, 12), rate = c(10, 20, 30, 30, 30))),
    n = 360, pi1 = 0.5, alpha = 0.025, times = c(19.7, 22, 28.8, 60), tchange = c(0, 3),
    acc = acc3,
    pw = list(rate10 = c(l0, l0), rate20 = c(0.8, 0.5) * l0, rate30 = rep(0.8 / 12, 2),
              rate40 = c(l0, l0), ratec0 = rep(0.3 / 12, 2), type0 = 3, rp20 = 0.7,
              rate11 = c(0.8, 0.5) * l0, rate21 = c(0.8, 0.5) * l0, rate31 = zero(2),
              ratec1 = rep(0.3 / 12, 2), type1 = 1))
  # Barthel et al. (2006), Table IV, rows 5 and 12 (years)
  lc <- log(2)
  dh <- function(x) -log(1 - x) / 4
  for (rw in list(list(id = "r05", c1 = 0, n = 489), list(id = "r12", c1 = 0.3, n = 606))) {
    ctl <- if (rw$c1 > 0) {
      ms_model(death = lc, prog = dh(rw$c1), death_prog = 0.7 * lc, clock = "forward")
    } else {
      ms_model(death = lc)
    }
    trt <- ms_model(death = 0.7 * lc, prog = dh(0.3), death_prog = lc, clock = "forward")
    cases[[paste0("Barthel2006-T4-", rw$id)]] <- list(
      scenario = ms_scenario(ctl, trt, accrual = ms_accrual(time = 0, prop = 1, duration = 2)),
      n = rw$n, pi1 = 0.5, alpha = 0.025, times = c(3, 4), tchange = c(0, 2),
      acc = list(taur = 2, ut = 2, u = 0.5),
      pw = list(rate10 = c(lc, lc), rate20 = c(0.7, 0.7) * lc, rate30 = rep(dh(rw$c1), 2),
                ratec0 = zero(2), type0 = 1,
                rate11 = c(0.7, 0.7) * lc, rate21 = c(lc, lc), rate31 = rep(dh(0.3), 2),
                ratec1 = zero(2), type1 = 1))
  }

  rows <- list()
  for (cn in names(cases)) {
    cs <- cases[[cn]]
    pw <- cs$pw
    get <- function(nm, default) if (is.null(pw[[nm]])) default else pw[[nm]]
    common <- list(taur = cs$acc$taur, u = cs$acc$u, ut = cs$acc$ut, pi1 = cs$pi1,
                   rate11 = pw$rate11, rate21 = pw$rate21, rate31 = pw$rate31,
                   rate41 = get("rate41", pw$rate21), rate51 = get("rate51", pw$rate21),
                   ratec1 = pw$ratec1,
                   rate10 = pw$rate10, rate20 = pw$rate20, rate30 = pw$rate30,
                   rate40 = get("rate40", pw$rate20), rate50 = get("rate50", pw$rate20),
                   ratec0 = pw$ratec0, tchange = cs$tchange, type1 = pw$type1, type0 = pw$type0,
                   rp21 = get("rp21", 0.5), rp20 = get("rp20", 0.5), eps = eps)
    grp <- function(g) {
      if (g == "control") {
        list(rate1 = common$rate10, rate2 = common$rate20, rate3 = common$rate30,
             rate4 = common$rate40, rate5 = common$rate50, ratec = common$ratec0,
             type = common$type0, rp2 = common$rp20)
      } else {
        list(rate1 = common$rate11, rate2 = common$rate21, rate3 = common$rate31,
             rate4 = common$rate41, rate5 = common$rate51, ratec = common$ratec1,
             type = common$type1, rp2 = common$rp21)
      }
    }
    scn <- cs$scenario
    m <- lr_moments(scn, n = cs$n, calendar_time = cs$times, variance = "alternative")
    ev <- ms_events(scn, n = cs$n, calendar_time = cs$times)
    za <- stats::qnorm(1 - cs$alpha)
    pkg <- list(events_control = ev$events_control, events_treatment = ev$events_treatment,
                score_mean = m$score_mean, info = m$info, score_var = m$z_var * m$info,
                power_alternative = stats::pnorm((m$z_mean - za) / sqrt(m$z_var)),
                power_null = stats::pnorm(m$z_mean - za))
    # each PWEALL function is called separately, so that an error in one
    # leaves the values of the others
    msgs <- character(0)
    try_pw <- function(expr) {
      tryCatch(expr, error = function(e) {
        msgs <<- c(msgs, conditionMessage(e))
        NULL
      })
    }
    pwv <- list()
    e0 <- try_pw(do.call(PWEALL::pwecxpwu, c(list(t = cs$times, taur = common$taur, u = common$u,
                                                  ut = common$ut, tchange = common$tchange,
                                                  eps = eps), grp("control")))$du)
    e1 <- try_pw(do.call(PWEALL::pwecxpwu, c(list(t = cs$times, taur = common$taur, u = common$u,
                                                  ut = common$ut, tchange = common$tchange,
                                                  eps = eps), grp("treatment")))$du)
    if (!is.null(e0)) pwv$events_control <- cs$n * (1 - cs$pi1) * e0
    if (!is.null(e1)) pwv$events_treatment <- cs$n * cs$pi1 * e1
    ov <- try_pw(lapply(cs$times, function(tt) {
      do.call(PWEALL::overallvar, c(list(tfix = tt, veps = eps, beta = 0), common))
    }))
    if (!is.null(ov)) {
      pwv$score_mean <- cs$n * vapply(ov, function(o) o$EA, numeric(1))
      pwv$info <- cs$n * vapply(ov, function(o) o$xdenom, numeric(1))
      pwv$score_var <- cs$n * vapply(ov, function(o) o$vs, numeric(1))
    }
    pp <- try_pw(do.call(PWEALL::pwepower, c(list(t = cs$times, alpha = cs$alpha, twosided = 0,
                                                  veps = eps, n = cs$n), common))$power)
    if (!is.null(pp)) {
      pwv$power_alternative <- pp[, 2]
      pwv$power_null <- pp[, 3]
    }
    pwv$msg <- paste(unique(msgs), collapse = "; ")
    for (q in names(pkg)) {
      for (i in seq_along(cs$times)) {
        val <- if (is.null(pwv[[q]])) NA_real_ else pwv[[q]][i]
        crit <- if (grepl("^power", q)) "abs:0.002" else "rel:0.002"
        if (!is.finite(val)) {
          crit <- "info"
        }
        rows[[length(rows) + 1L]] <- data.frame(
          id = paste("pweall", cn, format(cs$times[i]), q, sep = "|"), source = "pweall",
          scenario = cn, hypothesis = "alternative", look_type = "calendar", look = as.character(i),
          time = cs$times[i], quantity = q, endpoint = "os", group = "", variant = "",
          criterion = crit, simulated = val, se = 0, n = NA_real_, package = pkg[[q]][i],
          note = if (nzchar(pwv$msg)) paste("PWEALL error:", pwv$msg) else if (!is.finite(val)) "PWEALL value not finite" else "",
          stringsAsFactors = FALSE)
      }
    }
  }
  res <- do.call(rbind, rows)
  rownames(res) <- NULL
  res
}
