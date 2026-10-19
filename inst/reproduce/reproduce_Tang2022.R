# Reproduction of Tang (2022), Complex survival trial design by the product
# integration method, Statistics in Medicine (doi:10.1002/sim.9256),
# Example 4, Tables 1 and 2.
#
# Times are in years. Control: progression hazard 0.429, death before
# progression 0.089, death after progression 0.089 * 11.64 (measured from
# progression), with a proportion of 0 or 0.4 switching at progression and a
# hazard ratio of 0.7 from 0.4 years after the switch. Experimental: hazard
# ratio 0.6 for progression after a lag, 0.7 for death (before and after
# progression) from 0.4 years after randomization. Dropout hazard 0.02,
# accrual uniform over 2 years, 1:1, one-sided 2.5%, power 80%.
#
# The trial duration tau is the time at which the proportion of patients
# with an event reaches the data maturity (0.65 or 0.80). The sample size is
# that of the log-rank test at tau, and the nominal power is computed with
# the published sample size at tau.
#
# The text prints a progression lag of 0.025 years (D0.025(0.6)). The rows
# are judged with 0.025 years; a lag of 0.25 years is the evidence of
# Tang2022-Ex4-lag (see explained.csv).

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "Tang2022"
pv <- repro_published(src)
pub <- function(id) pv$value[pv$id == id]
hpost <- 0.089 * 11.64
acc <- repro_unif(2)
a <- 0.025

trt_model <- function(lag) {
  tm <- sort(unique(c(0, lag, 0.4)))
  th_p <- ifelse(tm < lag - 1e-12, 1, 0.6)
  th_d <- ifelse(tm < 0.4 - 1e-12, 1, 0.7)
  ms_model(time = tm, death = 0.089 * th_d, prog = 0.429 * th_p,
           death_prog = hpost * th_d, time_post = tm, dropout = 0.02, clock = "forward")
}
ctl_model <- function(ps) {
  ms_model(death = 0.089, prog = 0.429, death_prog = hpost, switch_prob = ps,
           death_switch = c(hpost, 0.7 * hpost), time_post = c(0, 0.4),
           dropout = 0.02, clock = "reset")
}
design <- function(ps, de, lag, n_pub) {
  scn <- ms_scenario(ctl_model(ps), trt_model(lag), accrual = acc)
  tau <- ms_event_time(scn, n = 1, events = de)
  n <- ms_size(scn, looks = ms_looks("calendar", time = tau), test = ms_test(alpha = a),
               power = 0.8, solve_for = "n")$solution
  pw <- 100 * repro_power_fixed(scn, n_pub, tau, a)
  c(tau = tau, n = n, power = pw)
}

rows <- list()
for (sw in c("NoSwi", "P40Swi")) {
  ps <- if (sw == "NoSwi") 0 else 0.4
  for (de in c(0.65, 0.80)) {
    p <- sprintf("Tang2022-Ex4-%s-%.2f-", sw, de)
    n_pub <- pub(paste0(p, "n"))
    d1 <- design(ps, de, 0.025, n_pub)
    d2 <- design(ps, de, 0.25, n_pub)
    rows[[paste0(p, "tau")]] <- repro_row(paste0(p, "tau"), d1[["tau"]], evidence = d2[["tau"]])
    rows[[paste0(p, "n")]] <- repro_row(paste0(p, "n"), d1[["n"]], evidence = d2[["n"]])
    for (q in c("nomPI", "nomRI")) {
      rows[[paste0(p, q)]] <- repro_row(paste0(p, q), d1[["power"]], evidence = d2[["power"]])
    }
  }
}

result <- repro_finish(src, rows)
