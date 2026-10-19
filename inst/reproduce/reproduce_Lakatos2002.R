# Reproduction of Lakatos (2002), Designing complex group sequential survival
# trials, Statistics in Medicine 21:1969-1989, Tables VIII and IX.
#
# Times are in months. Yearly rates are yearly probabilities, as in the
# article (Section 4: 1 - (1 - 0.04035)^12 = 0.39), so a yearly rate r gives
# the monthly hazard -log(1 - r) / 12. One-sided 2.5%, power 90%, Lan-DeMets
# spending of O'Brien-Fleming and Pocock type with the event fractions as
# information fractions. The article computes with a monthly discrete-time
# Markov chain, so all items are INFO.

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "Lakatos2002"
pv <- repro_published(src)
pub <- function(id) pv$value[pv$id == id]
f <- function(r) -log(1 - r) / 12
a <- 0.025
# The early looks of the 6-monthly designs have very small information
# fractions, so that almost no alpha is spent there (Table VI prints the
# bound 5.0000 at months 6-18). An error of the bound computation is reported
# in the note instead of stopping the run.
try_num <- function(f) {
  tryCatch(list(value = f(), msg = ""),
           error = function(e) list(value = NA_real_, msg = paste("error:", conditionMessage(e))))
}
n_gs <- function(scn, looks, sf) {
  test <- if (is.null(sf)) ms_test(alpha = a) else ms_test(alpha = a, sf_eff = sf)
  try_num(function() {
    ms_size(scn, looks = ms_looks("calendar", time = looks), test = test, power = 0.9,
            solve_for = "n")$solution
  })
}
pow_gs <- function(scn, n, looks, sf) {
  try_num(function() {
    ms_power(scn, n = n, looks = ms_looks("calendar", time = looks),
             test = ms_test(alpha = a, sf_eff = sf))$power
  })
}
row_try <- function(id, r, note = "") {
  repro_row(id, r$value, note = trimws(paste(note, r$msg)))
}
rows <- list()

# Table IX: accrual uniform over 12 months, yearly analyses over 60 months
acc9 <- repro_unif(12)
ctl9 <- ms_model(death = f(0.09))
trt9 <- list(trial1 = ms_model(time = c(0, 12, 24, 36, 48),
                               death = f(c(0.03, 0.03, 0.08, 0.08, 0.08))),
             trial2 = ms_model(time = c(0, 12, 24, 36, 48),
                               death = f(c(0.08, 0.08, 0.056, 0.03, 0.03))))
yearly <- c(12, 24, 36, 48, 60)
for (tr in names(trt9)) {
  scn <- ms_scenario(ctl9, trt9[[tr]], accrual = acc9)
  id <- sprintf("Lakatos2002-T9-%s-fixed-n", tr)
  rows[[id]] <- row_try(id, n_gs(scn, 60, NULL))
  for (sf in c("OBF", "Pocock")) {
    sfp <- if (sf == "OBF") "OF" else "Pocock"
    id <- sprintf("Lakatos2002-T9-%s-%s-n", tr, sf)
    rows[[id]] <- row_try(id, n_gs(scn, yearly, sfp))
    id_p <- sprintf("Lakatos2002-T9-%s-%s-power", tr, sf)
    rows[[id_p]] <- row_try(id_p, pow_gs(scn, pub(id), yearly, sfp),
                            note = "power at the published sample size")
  }
}

# Table VIII, exponential column: 27% per year, (1 - 0.225) times the
# control rate, accrual uniform over 24 months. The analysis times of the
# group sequential designs are not given; 6-monthly analyses (as in Table VI)
# are used, and the sample sizes with yearly analyses are given in the note.
scn8 <- ms_scenario(ms_model(death = f(0.27)), ms_model(death = f(0.775 * 0.27)),
                    accrual = repro_unif(24))
six <- seq(6, 60, by = 6)
rows[["Lakatos2002-T8-exp-fixed-n"]] <- row_try("Lakatos2002-T8-exp-fixed-n",
                                                 n_gs(scn8, 60, NULL))
for (sf in c("OBF", "Pocock")) {
  sfp <- if (sf == "OBF") "OF" else "Pocock"
  id <- sprintf("Lakatos2002-T8-exp-%s-n", sf)
  ny <- n_gs(scn8, yearly, sfp)
  rows[[id]] <- row_try(id, n_gs(scn8, six, sfp),
                        note = sprintf("6-monthly analyses; yearly analyses: %.1f %s",
                                       ny$value, ny$msg))
}

# Table VIII, Markov column (RALES): control rates from the last column of
# Table I (39, 26, 25, 23, 20 and 20% per year in months 0-3, 3-6, 6-12,
# 12-24, 24-36 and 36-60), experimental rates (1 - 0.225) times the control
# rates, non-compliance 10% in the first year and 5% afterwards, drop-in 5%
# per year, no loss to follow-up, recruitment of Table V (relative monthly
# rates 10, 20, 40, 60, 80 and 100% in months 0-3, 3-6, 6-9, 9-12, 12-15 and
# 15-24), 6-monthly analyses (Tables VI and VII). In the article the patients
# may switch back and forth any number of times; here at most two switches
# (non-compliance or drop-in, then back) are allowed.
tmr <- c(0, 3, 6, 12, 24, 36)
yr_rate <- c(0.39, 0.26, 0.25, 0.23, 0.20, 0.20)
hc <- f(yr_rate)
he <- f(0.775 * yr_rate)
ncomp <- ifelse(tmr < 12, f(0.10), f(0.05))
din <- rep(f(0.05), length(tmr))
ctl_r <- ms_model(time = tmr, prog = din, death = hc, death_prog = he,
                  switch_rate = ncomp, death_switch = hc, time_post = tmr,
                  clock = "forward")
trt_r <- ms_model(time = tmr, prog = ncomp, death = he, death_prog = hc,
                  switch_rate = din, death_switch = he, time_post = tmr,
                  clock = "forward")
rel <- c(10, 20, 40, 60, 80, 100)
len <- c(3, 3, 3, 3, 3, 9)
acc_r <- ms_accrual(time = c(0, 3, 6, 9, 12, 15), prop = rel * len / sum(rel * len),
                    duration = 24)
scn_r <- ms_scenario(ctl_r, trt_r, accrual = acc_r)
nt <- "approximation with at most two switches"
rows[["Lakatos2002-T8-Markov-fixed-n"]] <- row_try("Lakatos2002-T8-Markov-fixed-n",
                                                    n_gs(scn_r, 60, NULL), note = nt)
for (sf in c("OBF", "Pocock")) {
  sfp <- if (sf == "OBF") "OF" else "Pocock"
  id <- sprintf("Lakatos2002-T8-Markov-%s-n", sf)
  rows[[id]] <- row_try(id, n_gs(scn_r, six, sfp), note = nt)
}

result <- repro_finish(src, rows)
