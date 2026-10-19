# Reproduction of the design of innovaTV 301 (Vergote et al. 2024, New
# England Journal of Medicine, doi:10.1056/NEJMoa2313811): statistical
# analysis plan version 2.0, Sections 6.2 and 6.7 (computed with EAST).
#
# OS: hazard ratio 0.7, two-sided 5% level, power 90%, interim analysis at
# 75% of the events, Lan-DeMets spending of O'Brien-Fleming type. The number
# of events, the power and the bounds follow from the canonical form of the
# group sequential design (drift -log(0.7) sqrt(D t / 4) at information
# fraction t), computed with gs_bounds() and gs_prob().
# PFS: Lan-DeMets spending of Pocock type at 413 and 457 events.
# The approximate sample size, timing and PFS events (INFO) are computed with
# the exponential model of Fleischer et al. (2009) for PFS and OS (medians
# 3 and 9 months for control, 4.3 and 12.9 months for tisotumab vedotin),
# dropout of 5% per year and uniform accrual over 23 months (SAP Section
# 6.2: "accrual duration of approximately 23 months"); the accrual pattern is
# not given.

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "Vergote2024"
a <- 0.025
rows <- list()

# OS, canonical form
tf <- c(0.75, 1)
corr <- sqrt(outer(tf, tf, pmin) / outer(tf, tf, pmax))
b_os <- gs_bounds(tf, alpha = a, sf_eff = "OF")
gs_power <- function(d) {
  m <- -log(0.7) * sqrt(d * tf / 4)
  sum(gs_prob(b_os$lower, b_os$upper, mean = m, sigma = corr)$prob_upper)
}
d_req <- stats::uniroot(function(d) gs_power(d) - 0.9, c(100, 1000), tol = 1e-8)$root
rows[["Vergote2024-OS-events"]] <- repro_row("Vergote2024-OS-events", d_req)
rows[["Vergote2024-OS-power"]] <- repro_row("Vergote2024-OS-power", 100 * gs_power(336))
d_os <- c(252, 336)
nom_os <- 2 * b_os$nominal_p
hr_os <- exp(-b_os$upper / sqrt(d_os / 4))
rows[["Vergote2024-OS-nominal-IA"]] <- repro_row("Vergote2024-OS-nominal-IA", nom_os[1])
rows[["Vergote2024-OS-nominal-FA"]] <- repro_row("Vergote2024-OS-nominal-FA", nom_os[2])
# evidence: rounding to 3 decimals and then to 2 decimals (half up)
rows[["Vergote2024-OS-HR-IA"]] <- repro_row(
  "Vergote2024-OS-HR-IA", hr_os[1],
  evidence = repro_round_half_up(repro_round_half_up(hr_os[1], 3), 2),
  evidence_perturbed = repro_round_half_up(repro_round_half_up(hr_os[1], 3), 2))
rows[["Vergote2024-OS-HR-FA"]] <- repro_row("Vergote2024-OS-HR-FA", hr_os[2])

# PFS bounds at the OS analyses
d_pfs <- c(413, 457)
b_pfs <- gs_bounds(d_pfs / d_pfs[2], alpha = a, sf_eff = "Pocock")
nom_pfs <- 2 * b_pfs$nominal_p
hr_pfs <- exp(-b_pfs$upper / sqrt(d_pfs / 4))
rows[["Vergote2024-PFS-nominal-IA"]] <- repro_row("Vergote2024-PFS-nominal-IA", nom_pfs[1])
rows[["Vergote2024-PFS-nominal-FA"]] <- repro_row("Vergote2024-PFS-nominal-FA", nom_pfs[2])
rows[["Vergote2024-PFS-HR-IA"]] <- repro_row("Vergote2024-PFS-HR-IA", hr_pfs[1])
rows[["Vergote2024-PFS-HR-FA"]] <- repro_row("Vergote2024-PFS-HR-FA", hr_pfs[2])

# INFO: Fleischer model with uniform accrual over 23 months
dr <- -log(0.95) / 12
scn <- ms_scenario(ms_fleischer(median_pfs = 3, median_os = 9, dropout = dr),
                   ms_fleischer(median_pfs = 4.3, median_os = 12.9, dropout = dr),
                   accrual = repro_unif(23))
n_pub <- 482
t_os <- ms_event_time(scn, n = n_pub, events = d_os, endpoint = "os")
ev_pfs <- ms_events(scn, n = n_pub, calendar_time = t_os, endpoint = "pfs")$events
n_35 <- stats::uniroot(function(n) {
  ms_events(scn, n = n, calendar_time = 35)$events - 336
}, c(300, 700), tol = 1e-8)$root
pw_pfs <- ms_power(scn, n = n_pub, looks = ms_looks("os", fraction = c(0.75, 1)),
                   test = ms_test(alpha = a, endpoints = "pfs", sf_eff = "Pocock"),
                   events = 336)
nt <- "uniform accrual over 23 months assumed"
rows[["Vergote2024-n"]] <- repro_row("Vergote2024-n", n_35,
                                     note = paste0(nt, "; sample size with 336 OS events at 35 months"))
rows[["Vergote2024-time-IA"]] <- repro_row("Vergote2024-time-IA", t_os[1], note = nt)
rows[["Vergote2024-time-FA"]] <- repro_row("Vergote2024-time-FA", t_os[2], note = nt)
rows[["Vergote2024-PFS-events-IA"]] <- repro_row("Vergote2024-PFS-events-IA", ev_pfs[1], note = nt)
rows[["Vergote2024-PFS-events-FA"]] <- repro_row("Vergote2024-PFS-events-FA", ev_pfs[2], note = nt)
rows[["Vergote2024-PFS-power"]] <- repro_row(
  "Vergote2024-PFS-power", 100 * pw_pfs$power,
  note = paste0(nt, sprintf("; PFS events at the OS analyses %.1f and %.1f",
                           pw_pfs$looks$events[1], pw_pfs$looks$events[2])))

result <- repro_finish(src, rows)
