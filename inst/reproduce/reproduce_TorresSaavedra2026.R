# Reproduction of Torres-Saavedra et al. (2026), A note on rank-preserving
# structural failure time models to account for crossover, Clinical Trials
# 23:219-224: the probability of progression before death on the standard
# arm (Simulations section, p. 221) and the average proportions of censored
# observations under the null and alternative scenarios (Results section,
# p. 223).
#
# Model of the Simulations section: on the standard arm, the time to
# progression and the tentative overall survival T_S are independent and
# exponential with medians 20 and 40 months. A patient who progresses before
# death crosses over to the experimental treatment at the time of
# progression U and dies at U + (T_S - U) / c (model (1), the rank
# preserving structural failure time model). On the experimental arm,
# overall survival is exponential with median 40 / c, with log(c) = 0 (null
# scenario) or -0.40 (alternative scenario). In ms_model() the crossover is
# switch_prob = 1 with switch_af = c; since T_S is exponential, the hazard
# of death after progression without crossover is that of T_S. Accrual is
# uniform over 60 months with 60 months of additional follow-up, so that the
# data are analysed at month 120.
#
# The proportions of censored observations are printed as averages over
# 10,000 simulated trials with 100 patients per group. The expected
# proportions computed here do not depend on the number of patients; the
# Monte Carlo standard error of the printed averages is about 0.03
# percentage points.

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "TorresSaavedra2026"
rows <- list()
h_prog <- log(2) / 20
h_os <- log(2) / 40

models <- function(log_c) {
  cf <- exp(log_c)
  list(control = ms_model(prog = h_prog, death = h_os, death_prog = h_os,
                          switch_prob = 1, switch_af = cf),
       treatment = ms_model(death = cf * h_os))
}

# Simulations section (p. 221): probability of progression before death on
# the standard arm, the integral of the progression hazard times the
# probability of being progression-free
ctl <- models(0)$control
p_prog <- stats::integrate(function(t) h_prog * ms_surv(ctl, t, endpoint = "pfs")$surv,
                           0, Inf, rel.tol = 1e-10)$value
rows[["TorresSaavedra2026-Sim-prog"]] <- repro_row("TorresSaavedra2026-Sim-prog", 100 * p_prog)

# Results section (p. 223): proportions of censored observations (per cent)
# at month 120, both groups together and by group
censored <- function(log_c) {
  m <- models(log_c)
  scn <- ms_scenario(m$control, m$treatment, allocation = 1, accrual = repro_unif(60))
  ev <- ms_events(scn, n = 200, calendar_time = 120)
  100 * c(all = 1 - ev$events / ev$enrolled,
          control = 1 - ev$events_control / ev$enrolled_control,
          treatment = 1 - ev$events_treatment / ev$enrolled_treatment)
}
cn <- censored(0)
ca <- censored(-0.40)
c067 <- censored(log(0.67))
rows[["TorresSaavedra2026-Res-cens-null"]] <- repro_row(
  "TorresSaavedra2026-Res-cens-null", cn[["all"]],
  note = sprintf("standard arm %.2f, experimental arm %.2f", cn[["control"]],
                 cn[["treatment"]]))
rows[["TorresSaavedra2026-Res-cens-alt"]] <- repro_row(
  "TorresSaavedra2026-Res-cens-alt", ca[["all"]],
  note = sprintf("standard arm %.2f, experimental arm %.2f; with c = 0.67 (caption of Table 1) %.2f",
                 ca[["control"]], ca[["treatment"]], c067[["all"]]))

result <- repro_finish(src, rows)
