# Reproduction of Erdmann et al. (2025), Biometrical Journal 67:e70017,
# doi:10.1002/bimj.70017: Tables 2 and 3 (PFS and OS as coprimary endpoints,
# illness-death model with constant transition hazards).
#
# The settings follow the authors' R code (GlobalSigLevel.R for Table 2 and
# GroupSeq_Sec1_2IA.R to GroupSeq_Sec4_2IA.R for Table 3; simulations with
# simIDM, 10,000 trials): 800 patients per group, uniform accrual of 25
# (Table 2) or 100 (Table 3) patients per time unit in each group, that is
# over 32 or 8 time units (simIDM applies the accrual intensity to each group;
# the article states 25 patients per time unit), exponential dropout with
# probability 0.1 within 12 time units, and the transition hazards of Table 1
# (h01 progression, h02 death before progression, h12 death after
# progression). The log-rank tests are two-sided (survdiff) at the levels 1%
# for PFS and 4% for OS (Bonferroni); the package uses the one-sided levels
# 0.005 and 0.02, and the probability of rejecting in the wrong direction is
# ignored for the power.
#
# Table 2: PFS is analysed at d_PFS PFS events and OS at d_OS OS events.
#   Column 3, PFS (PASS/FAIL): Schoenfeld's number of events for the hazard
#   ratio of PFS (ratio of h01 + h02), ms_size() with mean = "schoenfeld".
#   Column 3, OS (INFO): the article uses the average hazard ratio of
#   Kalbfleisch and Prentice estimated from 100 simulated trials; the package
#   value is the number of events with mean = "schoenfeld", which uses the
#   event-weighted average of the log hazard ratio.
#   Column 4 (INFO): numbers of events for 80% power from the simulation;
#   ms_size() for each endpoint.
#   Column 6 (INFO): two-sided global type I error, simulated under the
#   global null hypothesis at the numbers of events of column 3 (before the
#   search of column 4); the package gives the normal probability with the
#   correlation of lr_cov() under ms_null() at the times of these numbers of
#   expected events.
#   Text (INFO): power of 88.5% with 770 OS events in Scenario 1 (Figure 6 of
#   the article shows 87.22% for the same design).
# Table 3: the final PFS analysis and an interim OS analysis at d_PFS PFS
#   events (Schoenfeld's number of rpact rounded up, with the hazard ratios of
#   PFS of the code: 0.72, 0.725, 0.764 and 0.8) and the final OS analysis at
#   d_OS OS events. The OS bounds are those of rpact for a Lan-DeMets design of
#   O'Brien-Fleming type with two-sided level 4% (equal to the one-sided design
#   at 2%) at the planned information fractions 0.4, 0.3, 0.4 and 0.25; the
#   simulation keeps them fixed.
#   Column 4 (PASS/FAIL): final number of OS events of rpact for 80% power
#   with the average hazard ratio of Table 1, from the canonical form with
#   gs_bounds() and gs_prob(), rounded up as in the code.
#   Column 6 (INFO): mean number of OS events at the PFS analysis; package:
#   expected OS events at the time of d_PFS expected PFS events.
#   Column 7 (INFO): final number of OS events for 80% power of OS (rejection
#   at either analysis) with the fixed bounds; package: lr_cov() with the
#   looks driven by PFS and OS events.
#   Columns 8 and 9 (INFO): two-sided type I error of OS and global type I
#   error under the global null hypothesis, at d_PFS and the final number of
#   OS events of rpact.
#
# The package values use variance = "null" (the default); the notes give the
# values with "studentized".
#
# Not reproduced: the empirical type I errors of each endpoint in Table 2
# (column 5; the nominal levels in the asymptotic approximation), the columns
# 2, 3 and 5 of Table 3 (numbers of events of rpact that the simulated design
# does not use; the printed values 644 for PFS in Scenario 3 and 279 for the
# interim OS analysis in Scenario 4 do not agree with the code, which gives
# 644.7 and 0.25 times 1118.3 rounded up, 280), and Figures 6 and 7, whose
# numbers are in images.

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "Erdmann2025"
n <- 1600
pw <- 0.8
a_pfs <- 0.005
a_os <- 0.02
c_pfs <- stats::qnorm(1 - a_pfs)
c_os <- stats::qnorm(1 - a_os)
dr <- -log(0.9) / 12
pv <- repro_published(src)
pub <- function(id) pv$value[pv$id == id]
rows <- list()

# Transition hazards (h01, h02, h12) of Table 1
erd_hazards <- list(
  list(control = c(0.10, 0.40, 0.30), treatment = c(0.06, 0.30, 0.30)),
  list(control = c(0.50, 0.30, 0.60), treatment = c(0.30, 0.28, 0.50)),
  list(control = c(0.180, 0.150, 0.255), treatment = c(0.140, 0.112, 0.250)),
  list(control = c(0.23, 0.07, 0.19), treatment = c(0.18, 0.06, 0.17)))
# Inputs of the code for Table 3: hazard ratios of PFS and OS (average) and
# planned information fraction of the interim OS analysis
erd_hr_pfs <- c(0.72, 0.725, 0.764, 0.8)
erd_hr_os <- c(0.812, 0.804, 0.821, 0.841)
erd_t_ia <- c(0.4, 0.3, 0.4, 0.25)

erd_scenario <- function(k, duration) {
  h <- erd_hazards[[k]]
  mk <- function(x) ms_model(prog = x[1], death = x[2], death_prog = x[3], dropout = dr)
  ms_scenario(mk(h$control), mk(h$treatment), accrual = repro_unif(duration))
}

# Expected number of events of one endpoint for 80% power with one analysis
erd_events <- function(scn, endpoint, alpha, ...) {
  ms_size(scn, looks = ms_looks(endpoint),
          test = ms_test(alpha = alpha, endpoints = endpoint, ...),
          power = pw, solve_for = "events", n = n)$solution
}

# Correlation matrix, under the global null hypothesis, of the statistics of
# the endpoints ep at the times of the expected numbers of events d of the
# endpoints trig
erd_null_corr <- function(scn, ep, d, trig) {
  s0 <- ms_null(scn)
  tm <- vapply(seq_along(d), function(i) {
    ms_event_time(s0, n, events = d[i], endpoint = trig[i])
  }, numeric(1))
  tt <- sort(unique(tm))
  jd <- lr_cov(s0, n, calendar_time = tt, endpoint = c("pfs", "os"), variance = "null")
  idx <- vapply(seq_along(ep), function(i) which(jd$endpoint == ep[i] & jd$time == tm[i]),
                integer(1))
  jd$sigma[idx, idx, drop = FALSE]
}

# Two-sided rejection probability under the null hypothesis: one minus the
# probability that every statistic stays within its bounds
erd_reject_two_sided <- function(bounds, r) {
  1 - as.numeric(mvtnorm::pmvnorm(lower = -bounds, upper = bounds, sigma = r,
                                  algorithm = mvtnorm::Miwa(steps = 512)))
}

# Table 2 (accrual over 32 time units)
for (k in 1:4) {
  id <- function(q) sprintf("%s-T2-S%d-%s", src, k, q)
  scn <- erd_scenario(k, 32)
  h <- erd_hazards[[k]]
  hr_pfs <- sum(h$treatment[1:2]) / sum(h$control[1:2])
  rows[[id("PFS-Schoenfeld")]] <- repro_row(
    id("PFS-Schoenfeld"), erd_events(scn, "pfs", a_pfs, mean = "schoenfeld"),
    note = sprintf("hazard ratio of PFS %.5f", hr_pfs))
  rows[[id("OS-Schoenfeld")]] <- repro_row(
    id("OS-Schoenfeld"), erd_events(scn, "os", a_os, mean = "schoenfeld"),
    note = paste("event-weighted average of the log hazard ratio",
                 "(the article: average hazard ratio of 100 simulated trials)"))
  for (ep in c("pfs", "os")) {
    a <- if (ep == "pfs") a_pfs else a_os
    q <- paste0(toupper(ep), "-IDM")
    rows[[id(q)]] <- repro_row(
      id(q), erd_events(scn, ep, a),
      note = sprintf("simulation; variance \"studentized\": %.1f",
                     erd_events(scn, ep, a, variance = "studentized")))
  }
  r0 <- erd_null_corr(scn, c("pfs", "os"),
                      c(pub(id("PFS-Schoenfeld")), pub(id("OS-Schoenfeld"))), c("pfs", "os"))
  rows[[id("alpha-global")]] <- repro_row(
    id("alpha-global"), 100 * erd_reject_two_sided(c(c_pfs, c_os), r0),
    note = sprintf("simulation; correlation of the PFS and OS statistics under the null hypothesis %.4f",
                   r0[1, 2]))
}

# Text of Section 5.3: power of OS with 770 OS events in Scenario 1
scn1 <- erd_scenario(1, 32)
p770 <- vapply(c("null", "studentized"), function(v) {
  ms_power(scn1, n = n, looks = ms_looks("os"), test = ms_test(alpha = a_os, variance = v),
           events = 770)$power
}, numeric(1))
rows[["Erdmann2025-T2-S1-OS-power770"]] <- repro_row(
  "Erdmann2025-T2-S1-OS-power770", 100 * p770[["null"]],
  note = sprintf("simulation (Figure 6: 87.22); variance \"studentized\": %.2f",
                 100 * p770[["studentized"]]))

# Table 3 (accrual over 8 time units)
for (k in 1:4) {
  id <- function(q) sprintf("%s-T3-S%d-%s", src, k, q)
  scn <- erd_scenario(k, 8)
  # Schoenfeld's number of PFS events of rpact, rounded up as in the code
  d_pfs <- ceiling((c_pfs + stats::qnorm(pw))^2 / (0.25 * log(erd_hr_pfs[k])^2))
  # rpact: OS design of O'Brien-Fleming type in the canonical form
  tf <- c(erd_t_ia[k], 1)
  b <- gs_bounds(tf, alpha = a_os, sf_eff = "OF")
  corr <- sqrt(outer(tf, tf, pmin) / outer(tf, tf, pmax))
  rp_power <- function(d) {
    sum(gs_prob(b$lower, b$upper, mean = -log(erd_hr_os[k]) * sqrt(d * tf / 4),
                sigma = corr)$prob_upper)
  }
  d_rp <- stats::uniroot(function(d) rp_power(d) - pw, c(100, 5000), tol = 1e-8)$root
  rows[[id("OS-FA-rpact")]] <- repro_row(
    id("OS-FA-rpact"), d_rp,
    note = sprintf("bounds %.4f and %.4f", b$upper[1], b$upper[2]))
  # Expected OS events at the PFS analysis
  t1 <- ms_event_time(scn, n, events = d_pfs, endpoint = "pfs")
  ev_os_t1 <- ms_events(scn, n, calendar_time = t1, endpoint = "os")$events
  rows[[id("OS-IA-IDM")]] <- repro_row(
    id("OS-IA-IDM"), ev_os_t1,
    note = sprintf("simulation; expected OS events at %d expected PFS events", d_pfs))
  # Final number of OS events for 80% power of OS with the fixed bounds
  os_power <- function(d, variance) {
    t2 <- ms_event_time(scn, n, events = d, endpoint = "os")
    jd <- lr_cov(scn, n, calendar_time = c(t1, t2), endpoint = "os", variance = variance,
                 trigger = c("pfs", "os"))
    sum(gs_prob(c(-Inf, -Inf), b$upper, mean = jd$mean, sigma = jd$sigma)$prob_upper)
  }
  d_hi <- 0.99 * ms_events(scn, n, calendar_time = 200, endpoint = "os")$events
  d_fa <- vapply(c("null", "studentized"), function(v) {
    stats::uniroot(function(d) os_power(d, v) - pw, c(1.05 * ev_os_t1, d_hi),
                   tol = 1e-6)$root
  }, numeric(1))
  rows[[id("OS-FA-IDM")]] <- repro_row(
    id("OS-FA-IDM"), d_fa[["null"]],
    note = sprintf("simulation; variance \"studentized\": %.1f", d_fa[["studentized"]]))
  # Type I errors under the global null hypothesis at d_PFS and the final
  # number of OS events of rpact (rounded up as in the code)
  r0 <- erd_null_corr(scn, c("pfs", "os", "os"), c(d_pfs, d_pfs, ceiling(d_rp)),
                      c("pfs", "pfs", "os"))
  rows[[id("alpha-OS")]] <- repro_row(
    id("alpha-OS"), 100 * erd_reject_two_sided(b$upper, r0[2:3, 2:3]),
    note = "simulation")
  rows[[id("alpha-global")]] <- repro_row(
    id("alpha-global"), 100 * erd_reject_two_sided(c(c_pfs, b$upper), r0),
    note = sprintf("simulation; correlation of PFS and OS at the PFS analysis %.4f", r0[1, 2]))
}

result <- repro_finish(src, rows)
