# Reproduction of Xia, George and Wang (2016), A multi-state model for
# designing clinical trials for testing overall survival allowing for
# crossover after progression, Statistics in Biopharmaceutical Research
# 8:12-21: the CALGB 30607 example (Section 3.2, Tables 3 and 4).
#
# Times are in months. Constant hazards: PFS medians 2.5 (control) and 3.9
# (experimental) months, a proportion 0.1 of the PFS events being deaths
# before progression; survival after progression with median 6.5 months for
# control without crossover, 9.0 months for control after crossover at
# progression (probability p) and 9.0 months for the experimental group.
# Accrual 4.8 patients per month, two-sided 5%. Apart from the ratio of the
# OS medians, the published values come from simulations, so all items are
# INFO.

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "Xia2016"
pv <- repro_published(src)
pub <- function(id) pv$value[pv$id == id]
l2 <- log(2)
a <- 0.025
za <- stats::qnorm(0.975)
zb <- stats::qnorm(0.9)

arm <- function(mpfs, om, mpost, p = 0, mswitch = NULL) {
  lp <- l2 / mpfs
  ms_model(death = om * lp, prog = (1 - om) * lp, death_prog = l2 / mpost,
           switch_prob = p, death_switch = if (is.null(mswitch)) NULL else l2 / mswitch,
           clock = "forward")
}
ctl <- function(p) arm(2.5, 0.1, 6.5, p, 9.0)
trt <- arm(3.9, 0.1, 9.0)
acc <- ms_accrual(time = 0, rate = 4.8)
os_median <- function(m) {
  stats::uniroot(function(t) ms_surv(m, t, "os")$surv - 0.5, c(0.1, 100), tol = 1e-10)$root
}
# Equation (2) of Appendix A (OS distribution function)
f_closed <- function(t, med, om, m3, m4, p) {
  l <- l2 / med
  la1 <- (1 - om) * l
  la2 <- om * l
  la3 <- l2 / m3
  la4 <- l2 / m4
  s <- la1 + la2
  1 - p * la1 / (s - la4) * exp(-la4 * t) - (1 - p) * la1 / (s - la3) * exp(-la3 * t) +
    exp(-s * t) / s * (p * la1 * la4 / (s - la4) + (1 - p) * la1 * la3 / (s - la3) - la2)
}
power_at_d <- function(p, d) {
  n <- round(1.1 * d)
  scn <- ms_scenario(ctl(p), trt, accrual = acc)
  t <- ms_event_time(scn, n = n, events = d)
  c(time = t, power = repro_power_fixed(scn, n, t, a))
}
plab <- c("0", "0.25", "0.50", "0.75", "1")
pval <- c(0, 0.25, 0.5, 0.75, 1)
rows <- list()

m_e <- os_median(trt)
m_e_closed <- stats::uniroot(function(t) f_closed(t, 3.9, 0.1, 9, 9, 0) - 0.5, c(0.1, 100),
                             tol = 1e-10)$root
for (i in seq_along(pval)) {
  p <- pval[i]
  scn <- ms_scenario(ctl(p), trt, accrual = acc)
  m_c <- os_median(ctl(p))
  m_c_closed <- stats::uniroot(function(t) f_closed(t, 2.5, 0.1, 6.5, 9, p) - 0.5,
                               c(0.1, 100), tol = 1e-10)$root
  ratio <- m_c / m_e
  id <- sprintf("Xia2016-median-ratio-p%s", plab[i])
  rows[[id]] <- repro_row(id, ratio,
                          note = sprintf("medians %.4f and %.4f; equation (2): %.4f and %.4f",
                                         m_c, m_e, m_c_closed, m_e_closed))
  # Table 3: 232 patients, analysis at 222 PFS events and 36 months later
  t1 <- ms_event_time(scn, n = 232, events = 222, endpoint = "pfs")
  for (k in 1:2) {
    tt <- if (k == 1) t1 else t1 + 36
    sfx <- if (k == 1) "PFS" else "3y"
    ev <- ms_events(scn, n = 232, calendar_time = tt)$events
    id <- sprintf("Xia2016-T3-p%s-deaths-%s", plab[i], sfx)
    rows[[id]] <- repro_row(id, ev, note = sprintf("analysis at %.2f months", tt))
    id <- sprintf("Xia2016-T3-p%s-power-%s", plab[i], sfx)
    rows[[id]] <- repro_row(id, repro_power_fixed(scn, 232, tt, a))
  }
  # Table 4: D* of the proportional hazards approach (equation (3)) and of the
  # multi-state model; sample size round(1.1 D)
  id <- sprintf("Xia2016-T4-p%s-PH-Dstar", plab[i])
  rows[[id]] <- repro_row(id, 4 * (za + zb)^2 / log(ratio)^2,
                          note = "equation (3) with the package's ratio of the OS medians")
  for (ap in c("PH", "MS")) {
    d_pub <- pub(sprintf("Xia2016-T4-p%s-%s-Dstar", plab[i], ap))
    r <- power_at_d(p, d_pub)
    id <- sprintf("Xia2016-T4-p%s-%s-power", plab[i], ap)
    rows[[id]] <- repro_row(id, r[["power"]], note = sprintf("at the published D = %d", d_pub))
    id <- sprintf("Xia2016-T4-p%s-%s-ET", plab[i], ap)
    rows[[id]] <- repro_row(id, r[["time"]], note = sprintf("at the published D = %d", d_pub))
  }
  d_ms <- pub(sprintf("Xia2016-T4-p%s-MS-Dstar", plab[i]))
  dstar <- stats::uniroot(function(d) power_at_d(p, d)[["power"]] - 0.9,
                          c(0.7 * d_ms, 1.3 * d_ms), tol = 0.5)$root
  id <- sprintf("Xia2016-T4-p%s-MS-Dstar", plab[i])
  rows[[id]] <- repro_row(id, dstar,
                          note = "package D (power 90% with 1.1 D patients), search tolerance 0.5")
}

result <- repro_finish(src, rows)
