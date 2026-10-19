# Reproduction of Lu (2021), Sample size calculation for logrank test and
# prediction of number of events over time, Pharmaceutical Statistics
# (doi:10.1002/pst.2069), computed with lrstat.
#
# Table 1: exponential survival, accrual uniform over 5 years, analysis at
# 10 years, no dropout, one-sided 2.5%, power 90%. Column L (Lakatos) with
# mean = "direct", column S (Schoenfeld) with mean = "schoenfeld"; the
# published values are rounded up to even numbers.
# Table 2: Weibull survival, approximated by piecewise exponential hazards on
# half-year intervals (INFO).
# Table 3: expected numbers of events over time (300 patients accrued over
# 1 year, hazards 1 and 2/3, dropout hazard 0.2).

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "Lu2021"
acc5 <- repro_unif(5)
lab <- function(x) if (abs(x - 0.667) < 1e-9) "0.667" else sprintf("%.2f", x)
n_lu <- function(ctl, trt, mean) {
  repro_n_fixed(ms_scenario(ctl, trt, accrual = acc5), t_end = 10, alpha = 0.025,
                power = 0.9, mean = mean)
}
rows <- list()

# Table 1
for (s10 in c(0.8, 0.2)) {
  lc <- -log(s10) / 10
  for (th in c(0.667, 0.5, 0.25)) {
    for (col in c("L", "S")) {
      mean <- if (col == "L") "direct" else "schoenfeld"
      id <- sprintf("Lu2021-T1-%s-%.1f-%s", col, s10, lab(th))
      n <- n_lu(ms_model(death = lc), ms_model(death = th * lc), mean)
      is667 <- abs(th - 0.667) < 1e-9
      n23 <- if (is667) n_lu(ms_model(death = lc), ms_model(death = 2 / 3 * lc), mean) else NA
      note <- if (is667) sprintf("theta = 2/3: %.2f", n23) else ""
      if (id == "Lu2021-T1-S-0.8-0.667") {
        rows[[id]] <- repro_row(id, n, note = note, evidence = n23)
      } else {
        rows[[id]] <- repro_row(id, n, note = note)
      }
    }
  }
}

# Table 2 (INFO): Weibull with shape kappa, S(10) = s10, half-year pieces
tm <- seq(0, 9.5, by = 0.5)
for (kp in c(3, 1 / 3)) {
  klab <- if (kp == 3) "k3" else "k1_3"
  for (s10 in c(0.8, 0.2)) {
    lc <- -log(s10) / 10^kp
    hz <- diff(lc * seq(0, 10, by = 0.5)^kp) / 0.5
    for (th in c(0.667, 0.5, 0.25)) {
      for (col in c("L", "S")) {
        mean <- if (col == "L") "direct" else "schoenfeld"
        id <- sprintf("Lu2021-T2-%s-%s-%.1f-%s", col, klab, s10, lab(th))
        n <- n_lu(ms_model(time = tm, death = hz), ms_model(time = tm, death = th * hz), mean)
        rows[[id]] <- repro_row(id, n, note = "half-year piecewise exponential approximation")
      }
    }
  }
}

# Table 3
scn3 <- ms_scenario(ms_model(death = 1, dropout = 0.2), ms_model(death = 2 / 3, dropout = 0.2),
                    accrual = repro_unif(1))
yr <- seq(0.5, 4, by = 0.5)
ev <- ms_events(scn3, n = 300, calendar_time = yr)
for (i in seq_along(yr)) {
  id <- sprintf("Lu2021-T3-y%.1f", yr[i])
  rows[[id]] <- repro_row(id, ev$events[i])
}

result <- repro_finish(src, rows)
