# Reproduction of Shih (1995), Sample size calculation for complex clinical
# trials with survival endpoints, Controlled Clinical Trials 16:395-407: the
# examples computed with SIZE.
#
# Time is measured in periods (6 months in Examples 1 and 3, 3 months in
# Example 2a, 1 month in Example 2b, 1 year in the osteosarcoma study). The
# hazards are piecewise constant by period: a probability P over m periods
# with hazard ratios r_1, ..., r_m between the periods gives the hazard
# -log(1 - P) r_j / sum(r) in period j. Switching from control to the
# experimental treatment (din) and back (dout) is represented with at most
# two switches, the hazard after a switch being measured from randomization
# (clock = "forward"); loss is dropout. Two-sided 5% (one-sided 2.5%), power
# 80%. All patients enter at time 0 unless stated otherwise.
#
# size_calc() in helpers.R is the SIZE algorithm (each period split into
# sbdv = 20 intervals; the hazard ratio of an interval is the ratio of the
# interval event probabilities).

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "Shih1995"
pp <- repro_per_period
za <- stats::qnorm(0.975)
zb <- stats::qnorm(0.8)
a <- 0.025

pkg_case <- function(hc, k, t_end, acc = NULL, alloc = 1, din = 0, dout = 0, loss = 0) {
  arms <- repro_switch_arms(hc, k, din = din, dout = dout, loss = loss)
  scn <- ms_scenario(arms$control, arms$treatment, allocation = alloc,
                     accrual = if (is.null(acc)) repro_simult(t_end) else acc)
  n <- repro_n_fixed(scn, t_end = t_end, alpha = a, power = 0.8)
  ev <- ms_events(scn, n = 1000, calendar_time = t_end)
  p <- scn$prob
  list(n = n, pc = ev$events_control / (1000 * p[[1L]]),
       pe = ev$events_treatment / (1000 * p[[2L]]), scn = scn)
}
size_n <- function(..., sbdv, max_cross = Inf) {
  size_calc(..., sbdv = sbdv, max_cross = max_cross, za = za, zb = zb)
}
rows <- list()
helpers <- list()
add_helper <- function(case, sbdv, max_cross, r) {
  helpers[[length(helpers) + 1L]] <<- data.frame(case = case, sbdv = sbdv,
                                                 max_cross = max_cross, n = r$n,
                                                 PC = r$PC, PE = r$PE)
}

# Example 1: 30% control event probability over 4 periods with hazard ratios
# 1, 2, 2, 2 between the periods; k = 0.65
hc1 <- pp(0.30, c(1, 2, 2, 2))
din1 <- pp(0.10, rep(1, 4))
dout1 <- pp(0.10, c(1, 1, 2, 2))
loss1 <- pp(0.10, rep(1, 4))

# 1a: no switching, no loss
r1a <- pkg_case(hc1, 0.65, 4)
c1a20 <- size_n(hc1, 0.65, sbdv = 20)
c1a2k <- size_n(hc1, 0.65, sbdv = 2000)
add_helper("Ex1a", 20, Inf, c1a20)
add_helper("Ex1a", 2000, Inf, c1a2k)
rows[["Shih1995-Ex1a-n"]] <- repro_row(
  "Shih1995-Ex1a-n", r1a$n, evidence = c1a20$n,
  evidence_extra = abs(r1a$n / c1a2k$n - 1) <= 1 / 2000,
  note = sprintf("SIZE algorithm: %.2f (sbdv 20), %.2f (sbdv 2000)", c1a20$n, c1a2k$n))

# 1c: switching in both directions and loss
r1c <- pkg_case(hc1, 0.65, 4, din = din1, dout = dout1, loss = loss1)
c1c20 <- size_n(hc1, 0.65, din = din1, dout = dout1, loss = loss1, sbdv = 20)
c1c2k <- size_n(hc1, 0.65, din = din1, dout = dout1, loss = loss1, sbdv = 2000, max_cross = 2)
c1c20m <- size_n(hc1, 0.65, din = din1, dout = dout1, loss = loss1, sbdv = 20, max_cross = 2)
add_helper("Ex1c", 20, Inf, c1c20)
add_helper("Ex1c", 2000, 2, c1c2k)
add_helper("Ex1c", 20, 2, c1c20m)
rows[["Shih1995-Ex1c-n"]] <- repro_row(
  "Shih1995-Ex1c-n", r1c$n, evidence = c1c20$n,
  evidence_extra = abs(r1c$n / c1c2k$n - 1) <= 1 / 2000,
  note = sprintf("SIZE algorithm: %.2f (sbdv 20), %.2f (sbdv 20, at most 2 switches), %.2f (sbdv 2000, at most 2 switches)",
                 c1c20$n, c1c20m$n, c1c2k$n))
rows[["Shih1995-Ex1c-pc"]] <- repro_row(
  "Shih1995-Ex1c-pc", r1c$pc, evidence = c1c20$PC,
  evidence_extra = abs(r1c$pc / c1c2k$PC - 1) <= 1 / 2000,
  note = sprintf("SIZE algorithm: %.5f (sbdv 20), %.5f (sbdv 2000)", c1c20$PC, c1c2k$PC))
rows[["Shih1995-Ex1c-pe"]] <- repro_row(
  "Shih1995-Ex1c-pe", r1c$pe, evidence = c1c20$PE,
  evidence_extra = abs(r1c$pe / c1c2k$PE - 1) <= 1 / 2000,
  note = sprintf("SIZE algorithm: %.5f (sbdv 20), %.5f (sbdv 2000)", c1c20$PE, c1c2k$PE))

# 1c with 2:1 allocation (control:treatment) (INFO)
r1c21 <- pkg_case(hc1, 0.65, 4, alloc = 0.5, din = din1, dout = dout1, loss = loss1)
c1c21 <- size_n(hc1, 0.65, din = din1, dout = dout1, loss = loss1, ratio_ce = 2, sbdv = 20)
add_helper("Ex1c 2:1", 20, Inf, c1c21)
rows[["Shih1995-Ex1c21-n"]] <- repro_row(
  "Shih1995-Ex1c21-n", r1c21$n, note = sprintf("SIZE algorithm: %.2f (sbdv 20)", c1c21$n))

# 1d (INFO): 5 periods, accrual uniform over the first 2 periods
hc1d <- -log(0.7) / 7 * c(1, 2, 2, 2, 2)
din1d <- rep(-log(0.9) / 4, 5)
dout1d <- -log(0.9) / 6 * c(1, 1, 2, 2, 2)
loss1d <- rep(-log(0.9) / 4, 5)
r1d <- pkg_case(hc1d, 0.65, 5, acc = repro_unif(2), din = din1d, dout = dout1d, loss = loss1d)
c1d <- size_n(hc1d, 0.65, din = din1d, dout = dout1d, loss = loss1d,
              accrual = c(0.5, 0.5, 0, 0, 0), sbdv = 20)
add_helper("Ex1d", 20, Inf, c1d)
rows[["Shih1995-Ex1d-n"]] <- repro_row(
  "Shih1995-Ex1d-n", r1d$n, note = sprintf("SIZE algorithm: %.2f (sbdv 20)", c1d$n))

# Example 2a: 8 periods of 3 months; hazard ratio 0.82 in the second period
# (input of Table 2; the text derives 0.825)
h2a <- rep(-log(0.7) / 8, 8)
k2a <- c(1, 0.82, rep(0.65, 6))
k2a5 <- c(1, 0.825, rep(0.65, 6))
r2a <- pkg_case(h2a, k2a, 8)
r2a5 <- pkg_case(h2a, k2a5, 8)
c2a20 <- size_n(h2a, k2a, sbdv = 20)
c2a2k <- size_n(h2a, k2a, sbdv = 2000)
add_helper("Ex2a", 20, Inf, c2a20)
add_helper("Ex2a", 2000, Inf, c2a2k)
rows[["Shih1995-Ex2a-n"]] <- repro_row(
  "Shih1995-Ex2a-n", r2a$n, evidence = c2a20$n,
  evidence_extra = abs(r2a$n / c2a2k$n - 1) <= 1 / 2000,
  note = sprintf("SIZE algorithm: %.2f (sbdv 20), %.2f (sbdv 2000); with 0.825: package %.2f",
                 c2a20$n, c2a2k$n, r2a5$n))

# Example 2b: 24 periods of 1 month
h2b <- rep(-log(0.7) / 24, 24)
k2b <- c(1, 0.9417, 0.8834, 0.825, 0.7667, 0.7084, rep(0.65, 18))
r2b <- pkg_case(h2b, k2b, 24)
rows[["Shih1995-Ex2b-n"]] <- repro_row("Shih1995-Ex2b-n", r2b$n)

# Example 3 (INFO): Example 1c extended to 12 periods, accrual of 200
# patients per period
ext <- function(v, n) c(v, rep(v[length(v)], n - length(v)))
hc3 <- ext(hc1, 12)
arms3 <- repro_switch_arms(hc3, 0.65, din = ext(din1, 12), dout = ext(dout1, 12),
                           loss = ext(loss1, 12))
acc3 <- ms_accrual(time = 0, rate = 200)
scn3 <- ms_scenario(arms3$control, arms3$treatment, accrual = acc3)
# 3a: 800 patients (4 periods of accrual); duration for power 80%
s3a <- ms_size(scn3, looks = ms_looks("os"), test = ms_test(alpha = a), power = 0.8,
               solve_for = "followup", n = 800)
rows[["Shih1995-Ex3a-duration"]] <- repro_row("Shih1995-Ex3a-duration",
                                              s3a$accrual_duration + s3a$solution)
# 3b: accrual continues until the end of the study
s3b <- ms_size(scn3, looks = ms_looks("os"), test = ms_test(alpha = a), power = 0.8,
               solve_for = "n", followup = 0)
rows[["Shih1995-Ex3b-duration"]] <- repro_row(
  "Shih1995-Ex3b-duration", s3b$solution / 200,
  note = sprintf("sample size %.1f", s3b$solution))

# Table 4: unconditional power over a two-point distribution of k (each with
# probability 0.5), 600 patients, 4 periods, 20% control event probability,
# 10% din, dout and loss
hc4 <- pp(0.20, rep(1, 4))
d4 <- pp(0.10, rep(1, 4))
pw4 <- function(k) {
  arms <- repro_switch_arms(hc4, k, din = d4, dout = d4, loss = d4)
  repro_power_fixed(ms_scenario(arms$control, arms$treatment, accrual = repro_simult(4)),
                    n = 600, t = 4, alpha = a)
}
k4 <- list(c(0.5, 0.5), c(0.55, 0.45), c(0.6, 0.4), c(0.65, 0.35), c(0.7, 0.3),
           c(0.75, 0.25), c(0.8, 0.2), c(0.85, 0.15), c(0.9, 0.1))
for (ks in k4) {
  id <- sprintf("Shih1995-T4-k%.2f-%.2f", ks[1], ks[2])
  rows[[id]] <- repro_row(id, 100 * mean(vapply(ks, pw4, numeric(1))))
}

# Osteosarcoma study: 45% event probability over 5 years, k = log(0.7) / log(0.55), 400 patients
hco <- pp(0.45, rep(1, 5))
ko <- log(0.7) / log(0.55)
arms_o <- repro_switch_arms(hco, ko)
rows[["Shih1995-Osteo-power"]] <- repro_row(
  "Shih1995-Osteo-power",
  100 * repro_power_fixed(ms_scenario(arms_o$control, arms_o$treatment, accrual = repro_simult(5)),
                          n = 400, t = 5, alpha = a),
  note = sprintf("k = %.4f", ko))

result <- repro_finish(src, rows, helpers = do.call(rbind, helpers))
