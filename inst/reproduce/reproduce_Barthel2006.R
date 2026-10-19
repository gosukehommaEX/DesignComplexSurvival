# Reproduction of Barthel et al. (2006), Statistics in Medicine 25:2521-2542,
# Tables II-V: sample sizes computed with ART.
#
# Accrual uniform over 2 years, analysis at 4 years, 1:1 allocation,
# exponential survival with a control median of 1 year, two-sided 5% level
# (one-sided 2.5%), power 90%. Loss to follow-up and crossover are given as
# proportions by the end of the study (4 years), so the intensities are
# -log(1 - L) / 4. After a crossover the hazard is that of the other group.

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "Barthel2006"
lc <- log(2)
dh <- function(x) -log(1 - x) / 4
acc <- repro_unif(2)
n_req <- function(ctl, trt) {
  repro_n_fixed(ms_scenario(ctl, trt, accrual = acc), t_end = 4, alpha = 0.025,
                power = 0.9)
}
rows <- list()

# Table II: loss to follow-up (L1 control, L2 experimental)
t2 <- data.frame(hr = c(0.7, 0.8, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.8, 0.8, 0.8),
                 l1 = c(0, 0, 0.05, 0.2, 0.05, 0.2, 0.3, 0.4, 0.5, 0.3, 0.4, 0.5),
                 l2 = c(0, 0, 0.2, 0.05, 0.05, 0.2, 0.3, 0.4, 0.5, 0.3, 0.4, 0.5))
for (i in seq_len(nrow(t2))) {
  n <- n_req(ms_model(death = lc, dropout = dh(t2$l1[i])),
             ms_model(death = t2$hr[i] * lc, dropout = dh(t2$l2[i])))
  id <- sprintf("Barthel2006-T2-r%02d", i)
  if (i == 2L) {
    # evidence: the same design is printed as 1015 in Tables II-IV (Unadjusted)
    rows[[id]] <- repro_row(id, n, evidence = n,
                            evidence_ref_id = "Barthel2006-Unadj-HR0.8")
  } else {
    rows[[id]] <- repro_row(id, n)
  }
}

# Unadjusted sample sizes (no loss, no crossover)
for (hr in c(0.6, 0.7, 0.8, 0.9)) {
  id <- sprintf("Barthel2006-Unadj-HR%.1f", hr)
  rows[[id]] <- repro_row(id, n_req(ms_model(death = lc), ms_model(death = hr * lc)))
}

# Table III: hazard ratio HR1 in years 0-2 and HR2 from year 2
t3 <- data.frame(h1 = c(0.6, 0.6, 0.6, 0.7, 0.8, 0.8),
                 h2 = c(0.9, 0.8, 0.7, 0.8, 0.7, 0.6))
for (i in seq_len(nrow(t3))) {
  id <- sprintf("Barthel2006-T3-r%02d", i)
  trt <- ms_model(time = c(0, 2), death = c(t3$h1[i], t3$h2[i]) * lc)
  rows[[id]] <- repro_row(id, n_req(ms_model(death = lc), trt))
}

# Table IV: crossover (C1 control to experimental, C2 experimental to control)
t4 <- data.frame(hr = c(0.6, 0.6, 0.6, 0.6, 0.7, 0.8, 0.9, 0.7, 0.7, 0.7, 0.7, 0.7),
                 c1 = c(0, 0, 0, 0, 0, 0, 0, 0.1, 0.2, 0.3, 0.2, 0.3),
                 c2 = c(0.05, 0.1, 0.2, 0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1, 0.2, 0.3))
for (i in seq_len(nrow(t4))) {
  id <- sprintf("Barthel2006-T4-r%02d", i)
  hr <- t4$hr[i]
  ctl <- if (t4$c1[i] > 0) {
    ms_model(death = lc, prog = dh(t4$c1[i]), death_prog = hr * lc, clock = "forward")
  } else {
    ms_model(death = lc)
  }
  trt <- ms_model(death = hr * lc, prog = dh(t4$c2[i]), death_prog = lc,
                  clock = "forward")
  rows[[id]] <- repro_row(id, n_req(ctl, trt))
}

# Table V: combinations. The experimental hazard ratio is HR1 in years 0-2
# and HR2 from year 2. For crossover from control (rows 9-11), the hazard
# ratio pattern is measured from the time of the switch (clock = "reset");
# the values with the pattern measured from randomization (clock =
# "forward") are given in the note.
t5 <- data.frame(h1 = c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.8, 0.8, 0.6, 0.6, 0.6),
                 h2 = c(0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 0.6, 0.6, 0.8, 0.8, 0.8),
                 l = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2),
                 c1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0.1, 0.2),
                 c2 = c(0.2, 0.3, 0.2, 0.3, 0.2, 0.3, 0.2, 0.3, 0.1, 0.2, 0.2))
for (i in seq_len(nrow(t5))) {
  id <- sprintf("Barthel2006-T5-r%02d", i)
  he <- c(t5$h1[i], t5$h2[i]) * lc
  trt <- ms_model(time = c(0, 2), death = he, prog = dh(t5$c2[i]), death_prog = lc,
                  dropout = dh(t5$l[i]), clock = "forward")
  if (t5$c1[i] > 0) {
    ctl <- ms_model(death = lc, prog = dh(t5$c1[i]), death_prog = he,
                    time_post = c(0, 2), dropout = dh(t5$l[i]), clock = "reset")
    ctl_fwd <- ms_model(time = c(0, 2), death = lc, prog = dh(t5$c1[i]),
                        death_prog = he, time_post = c(0, 2), dropout = dh(t5$l[i]),
                        clock = "forward")
    n_fwd <- n_req(ctl_fwd, trt)
    rows[[id]] <- repro_row(id, n_req(ctl, trt),
                            note = sprintf("clock forward for the control group: %.2f", n_fwd))
  } else {
    ctl <- ms_model(death = lc, dropout = dh(t5$l[i]))
    rows[[id]] <- repro_row(id, n_req(ctl, trt))
  }
}

result <- repro_finish(src, rows)
