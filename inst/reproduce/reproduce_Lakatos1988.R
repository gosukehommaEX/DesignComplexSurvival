# Reproduction of Lakatos (1988), Sample sizes based on the log-rank
# statistic in complex clinical trials, Biometrics 44:229-241: the example of
# Section 3 and the log-rank columns of Tables 2 and 3.
#
# Times are in years. Yearly rates of loss, non-compliance and drop-in are
# yearly probabilities p, with hazards -log(1 - p). Switching between the
# treatments (non-compliance and drop-in) is represented with at most two
# switches, the hazard after a switch being measured from randomization.
# Two-sided 5% (one-sided 2.5%), power 90%. All patients enter at time 0
# unless stated otherwise.

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "Lakatos1988"
hp <- function(p) -log(1 - p)
a <- 0.025
za <- stats::qnorm(0.975)
zb <- stats::qnorm(0.9)
rows <- list()
helpers <- list()

# Section 3 example (Gail 1985): hazards 1 and 1/2 per year, loss 3%,
# non-compliance 4%, drop-in 5% per year, 2 years
arms <- repro_switch_arms(c(1, 1), 0.5, din = hp(0.05), dout = hp(0.04), loss = hp(0.03))
scn <- ms_scenario(arms$control, arms$treatment, accrual = repro_simult(2))
n <- repro_n_fixed(scn, t_end = 2, alpha = a, power = 0.9)
ev <- ms_events(scn, n = n, calendar_time = 2)
d <- ev$events
pc <- ev$events_control / (n / 2)
pe <- ev$events_treatment / (n / 2)
for (th in c("hazard", "prob")) {
  for (mc in c(Inf, 2)) {
    r <- size_calc(c(1, 1), 0.5, sbdv = 10, din = hp(0.05), dout = hp(0.04), loss = hp(0.03),
                   za = za, zb = zb, theta = th, max_cross = mc)
    helpers[[length(helpers) + 1L]] <- data.frame(case = "Section 3", theta = th,
                                                  sbdv = 10, max_cross = mc, d = r$d,
                                                  PC = r$PC, PE = r$PE)
  }
}
hz <- helpers[[1]]
hz2 <- helpers[[2]]
rows[["Lakatos1988-Ex-deaths"]] <- repro_row(
  "Lakatos1988-Ex-deaths", d,
  note = sprintf("SIZE algorithm (theta as hazard ratio, K = 10): %.3f, %.3f with at most 2 switches",
                 hz$d, hz2$d))
rows[["Lakatos1988-Ex-n"]] <- repro_row(
  "Lakatos1988-Ex-n", 2 * ceiling(d - 1e-9) / (pc + pe),
  note = sprintf("PC %.4f, PE %.4f at 2 years; package sample size %.2f", pc, pe, n))

# Tables 2 and 3, log-rank column, Adjust = No ---------------------------
accs <- function(trial) {
  if (trial == "CVD") {
    w <- c(0.3, 0.4, 0.6, 0.8, 4.0)
    list(simultaneous = list(t = 5, acc = repro_simult(5)),
         uniform = list(t = 6, acc = repro_unif(2)),
         nonuniform = list(t = 6, acc = ms_accrual(time = c(0, 0.25, 0.5, 0.75, 1),
                                                   prop = w / sum(w), duration = 2)))
  } else {
    w <- c(0.4, 0.6, 0.8, 1.0)
    list(simultaneous = list(t = 1.5, acc = repro_simult(1.5)),
         uniform = list(t = 2, acc = repro_unif(1)),
         nonuniform = list(t = 2, acc = ms_accrual(time = c(0, 0.25, 0.5, 0.75),
                                                   prop = w / sum(w), duration = 1)))
  }
}
for (trial in c("CVD", "cancer")) {
  if (trial == "CVD") {
    hc <- hp(0.016)
    he <- hp(0.0096)
    lag <- 2
  } else {
    hc <- 1
    he <- 0.5
    lag <- 1
  }
  ac <- accs(trial)
  for (model in c("PH", "LM2")) {
    trt <- if (model == "PH") ms_model(death = he) else ms_model(time = c(0, lag), death = c(hc, he))
    for (en in names(ac)) {
      id <- sprintf("Lakatos1988-%s-%s-No-%s", trial, model, en)
      scn <- ms_scenario(ms_model(death = hc), trt, accrual = ac[[en]]$acc)
      rows[[id]] <- repro_row(id, repro_n_fixed(scn, t_end = ac[[en]]$t, alpha = a, power = 0.9))
    }
  }
}

# Tables 2 and 3, log-rank column, proportional hazards, Adjust = Yes ------
loss <- hp(c(0.03, 0.032, 0.034, 0.036, 0.038))
nonc <- hp(c(0.07, 0.035, 0.035, 0.035, 0.035))
din <- hp(c(0.09, 0.045, 0.050, 0.055, 0.060))
arms_c <- repro_switch_arms(rep(hp(0.016), 5), hp(0.0096) / hp(0.016), din = din,
                            dout = nonc, loss = loss)
scn_c <- ms_scenario(arms_c$control, arms_c$treatment, accrual = repro_simult(5))
rows[["Lakatos1988-CVD-PH-Yes-simultaneous"]] <- repro_row(
  "Lakatos1988-CVD-PH-Yes-simultaneous", repro_n_fixed(scn_c, t_end = 5, alpha = a, power = 0.9))
arms_k <- repro_switch_arms(c(1, 1), 0.5, din = din[1:2], dout = nonc[1:2], loss = loss[1:2])
ac <- accs("cancer")
for (en in names(ac)) {
  id <- sprintf("Lakatos1988-cancer-PH-Yes-%s", en)
  scn <- ms_scenario(arms_k$control, arms_k$treatment, accrual = ac[[en]]$acc)
  rows[[id]] <- repro_row(id, repro_n_fixed(scn, t_end = ac[[en]]$t, alpha = a, power = 0.9))
}

result <- repro_finish(src, rows, helpers = do.call(rbind, helpers))
