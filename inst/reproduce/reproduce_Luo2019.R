# Reproduction of Luo et al. (2019), Statistics in Medicine 38:192-209
# (PWEALL), Examples 1-3.
#
# Times are in months. The power uses the variance of the log-rank score
# under the alternative (variance = "alternative"). Rows of the time at which
# a target power is reached are judged by the power at the published time;
# the group-wise numbers of events are compared at the published time. The
# computed times at which the target power is reached are given in the
# notes.

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "Luo2019"
pv <- repro_published(src)
pub <- function(id) pv$value[pv$id == id]
pow <- function(scn, n, t, alpha) {
  repro_power_fixed(scn, n, t, alpha, variance = "alternative")
}
t_target <- function(scn, n, alpha, target, lo, hi) {
  stats::uniroot(function(x) pow(scn, n, x, alpha) - target, c(lo, hi), tol = 1e-8)$root
}
rows <- list()

# Example 1 -----------------------------------------------------------------
tm <- c(0, 6, 12, 24)
lam0 <- c(0.09, 0.09, 0.07, 0.05) / 12
drop <- c(0, 0, 0.015, 0.015) / 12
disc <- c(0.04, 0.04, 0.03, 0.03) / 12
acc1 <- ms_accrual(time = 0:13, rate = c(100, 200, 300, 300, rep(400, 8), 300, 200))
ctl1 <- ms_model(time = tm, death = lam0, dropout = drop)
hr_d <- c(1, 0.75, 0.75, 0.75)
trt1 <- list(
  "constHR" = ms_model(time = tm, death = 0.8 * lam0, dropout = drop),
  "delay-w1" = ms_model(time = tm, death = hr_d * lam0, dropout = drop),
  "delay-w0" = ms_model(time = tm, death = hr_d * lam0, dropout = drop, prog = disc,
                        death_prog = lam0, time_post = tm, clock = "forward")
)
n1 <- 4600
a1 <- 0.025
for (lb in names(trt1)) {
  scn <- ms_scenario(ctl1, trt1[[lb]], accrual = acc1)
  t722 <- ms_event_time(scn, n = n1, events = 722)
  ev <- ms_events(scn, n = n1, calendar_time = t722)
  p <- sprintf("Luo2019-Ex1-%s-", lb)
  rows[[paste0(p, "t722")]] <- repro_row(paste0(p, "t722"), t722)
  rows[[paste0(p, "d722-C")]] <- repro_row(paste0(p, "d722-C"), ev$events_control)
  rows[[paste0(p, "d722-T")]] <- repro_row(paste0(p, "d722-T"), ev$events_treatment)
  if (lb != "constHR") {
    id <- paste0(p, "power722")
    rows[[id]] <- repro_row(id, 100 * pow(scn, n1, t722, a1))
    id <- paste0(p, "t85")
    tp <- pub(id)
    tc <- t_target(scn, n1, a1, 0.85, 25, 120)
    rows[[id]] <- repro_row(id, pow(scn, n1, tp, a1),
                            note = sprintf("85%% power reached at %.3f months", tc),
                            computed_perturbed = pow(scn, n1, 1.1 * tp, a1))
    ev <- ms_events(scn, n = n1, calendar_time = tp)
    rows[[paste0(p, "d85-C")]] <- repro_row(paste0(p, "d85-C"), ev$events_control)
    rows[[paste0(p, "d85-T")]] <- repro_row(paste0(p, "d85-T"), ev$events_treatment)
    m <- lr_moments(scn, n = n1, calendar_time = tp, variance = "alternative")
    nt <- "proxy exp(E(U) / V) and sqrt(var(U)) / V at the published time; not a Cox estimate"
    rows[[paste0(p, "HR85")]] <- repro_row(paste0(p, "HR85"), exp(m$score_mean / m$info),
                                           note = nt)
    rows[[paste0(p, "SE85")]] <- repro_row(paste0(p, "SE85"), sqrt(m$z_var / m$info),
                                           note = nt)
  }
}

# Example 2 -----------------------------------------------------------------
tm2 <- c(0, 1, 2, 3, 4, 5, 6, 12)
lam2 <- c(0.9, 0.75, 0.6, 0.55, 0.5, 0.4, 0.3, 0.1) / 12
acc2 <- ms_accrual(time = c(0, 3, 6, 9, 12), rate = c(10, 30, 60, 100, 100))
ctl2 <- ms_model(time = tm2, death = lam2, dropout = 0.03 / 12)
trt2 <- list(
  "nodil" = ms_model(time = tm2, death = 0.6 * lam2, dropout = 0.03 / 12),
  "dil" = ms_model(time = tm2, death = 0.6 * lam2, dropout = 0.03 / 12, prog = 0.2 / 12,
                   death_prog = lam2, time_post = tm2, clock = "forward")
)
n2 <- 900
a2 <- 0.005
for (lb in names(trt2)) {
  scn <- ms_scenario(ctl2, trt2[[lb]], allocation = 2, accrual = acc2)
  for (tg in c(80, 85, 90)) {
    p <- sprintf("Luo2019-Ex2-%s-p%d-", lb, tg)
    id <- paste0(p, "time")
    tp <- pub(id)
    tc <- t_target(scn, n2, a2, tg / 100, 8, 80)
    rows[[id]] <- repro_row(id, pow(scn, n2, tp, a2),
                            note = sprintf("%d%% power reached at %.3f months", tg, tc),
                            computed_perturbed = pow(scn, n2, 1.1 * tp, a2))
    ev <- ms_events(scn, n = n2, calendar_time = tp)
    # printed heading C+T: the first number is read as the control group;
    # evidence with the groups in the order T, C
    rows[[paste0(p, "first")]] <- repro_row(paste0(p, "first"), ev$events_control,
                                            evidence = ev$events_treatment)
    rows[[paste0(p, "second")]] <- repro_row(paste0(p, "second"), ev$events_treatment,
                                             evidence = ev$events_control)
  }
}

# Example 3 -----------------------------------------------------------------
# Hybrid crossover of equation (5): after progression in the control group
# the hazard is 0.7 times the hazard after crossover plus 0.3 times the
# control hazard (hazard mixture, measured from progression). The level is
# not given in the article; two-sided 5% is assumed.
acc3 <- ms_accrual(time = c(0, 3, 6, 9, 12), rate = c(10, 20, 30, 30, 30))
l0 <- 1.5 / 12
dp <- c(0.7 * 0.8 * l0 + 0.3 * l0, 0.7 * 0.5 * l0 + 0.3 * l0)
ctl3 <- ms_model(death = l0, dropout = 0.3 / 12, prog = 0.8 / 12, death_prog = dp,
                 time_post = c(0, 3), clock = "reset")
trt3 <- ms_model(time = c(0, 3), death = c(0.8, 0.5) * l0, dropout = 0.3 / 12)
scn3 <- ms_scenario(ctl3, trt3, accrual = acc3)
n3 <- 360
a3 <- 0.025
p60 <- pow(scn3, n3, 60, a3)
for (tg in c(80, 85, 90)) {
  p <- sprintf("Luo2019-Ex3-p%d-", tg)
  id <- paste0(p, "time")
  tp <- pub(id)
  pw <- pow(scn3, n3, tp, a3)
  pw_p <- pow(scn3, n3, 1.1 * tp, a3)
  rows[[id]] <- repro_row(id, pw, evidence = pw,
                          evidence_criterion = sprintf("power_at_time:%.2f:0.01", tg / 100),
                          note = sprintf("power at 60 months %.4f", p60),
                          computed_perturbed = pw_p, evidence_perturbed = pw_p)
  ev <- ms_events(scn3, n = n3, calendar_time = tp)
  rows[[paste0(p, "C")]] <- repro_row(paste0(p, "C"), ev$events_control,
                                      evidence = ev$events_control, evidence_criterion = "abs:2")
  rows[[paste0(p, "T")]] <- repro_row(paste0(p, "T"), ev$events_treatment,
                                      evidence = ev$events_treatment, evidence_criterion = "abs:2")
}

result <- repro_finish(src, rows)
