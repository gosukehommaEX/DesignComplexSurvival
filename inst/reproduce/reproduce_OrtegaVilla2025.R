# Reproduction of Ortega-Villa et al. (2025), Futility monitoring in clinical
# trials, Statistics in Medicine 44:e70157: the conditional powers of
# Section 3.1 (Figure 1), Section 3.2 (LUME-Lung 2) and Section 3.3 (reverse
# conditional power example).
#
# The article computes the conditional power with formula (3), with one
# remaining analysis at the bound 1.96: given the interim z at the
# information fraction t, the final Z has the mean z sqrt(t) + theta (1 - t)
# and the variance 1 - t, where theta is the mean of the final Z, either the
# originally hypothesized drift (theta = 1.96 + qnorm(power), 3.00 for 85%
# power) or the current trend z / sqrt(t).
#
# Design of the package: 1:1, exponential survival, mean = "schoenfeld",
# variance = "null", a final analysis at D = 400 OS events and an interim
# analysis at t D events. The spending function "user" with
# param_eff = c(0, 1) spends no alpha at the interim analysis, so that the
# final efficacy bound is qnorm(0.975) as in formula (3) (diagnostic). With
# the hazard ratio exp(-theta / 10) the mean of the final Z is
# -log(HR) sqrt(D / 4) = theta (diagnostic). The originally hypothesized
# effect is ms_cond_power(under = "hr") with this hazard ratio, whose mean
# and correlation sqrt(t) of the Brownian motion are those of formula (3);
# the current trend is under = "trend". The value of under = "design", which
# takes the correlation from the information under the design scenario
# (variance = "null") instead of sqrt(t) and differs from the conditional
# power of formula (3) by 0.001 to 0.002, is given in the note.
#
# The conditional power 0.103 under the current trend in Section 3.2 was
# taken by the authors from Hanna et al. and used to deduce the interim
# z = 0.733, so it is checked for consistency only (INFO).

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "OrtegaVilla2025"
rows <- list()
diag <- list()

ov_design <- function(theta, t) {
  hr <- exp(-theta / 10)
  lam <- log(2) / 12
  scn <- ms_scenario(ms_model(death = lam), ms_model(death = hr * lam),
                     accrual = ms_accrual(rate = 20))
  x <- ms_power(scn, n = 1000, looks = ms_looks("os", fraction = c(t, 1)),
                test = ms_test(alpha = 0.025, sf_eff = "user", param_eff = c(0, 1),
                               mean = "schoenfeld", variance = "null"),
                events = 400)
  list(x = x, hr = hr)
}

cases <- list(
  list(id = "S3.1", theta = 3, t = 0.5, z = 0.3),
  list(id = "S3.2", theta = 3.24, t = 0.485, z = 0.733),
  list(id = "S3.3-85", theta = 3, t = 0.5, z = 0),
  list(id = "S3.3-80", theta = 2.8, t = 0.5, z = 0)
)
for (cs in cases) {
  d <- ov_design(cs$theta, cs$t)
  x <- d$x
  diag[[length(diag) + 1L]] <- repro_diag(
    sprintf("%s-%s-final-bound", src, cs$id), x$looks$z_upper[2], stats::qnorm(0.975), 1e-6,
    note = "final efficacy bound against 1.96 of formula (3)")
  diag[[length(diag) + 1L]] <- repro_diag(
    sprintf("%s-%s-final-mean", src, cs$id), x$looks$z_mean[2], cs$theta, 1e-10,
    note = "mean of the final Z against theta")
  cp_hr <- ms_cond_power(x, 1, cs$z, under = "hr", hr = d$hr)$cond_power
  cp_design <- ms_cond_power(x, 1, cs$z, under = "design")$cond_power
  cp_trend <- ms_cond_power(x, 1, cs$z, under = "trend")$cond_power
  note_d <- sprintf("under = \"design\": %.5f", cp_design)
  if (cs$id %in% c("S3.1", "S3.2")) {
    id <- sprintf("%s-%s-original", src, cs$id)
    rows[[id]] <- repro_row(id, cp_hr, note = note_d)
    id <- sprintf("%s-%s-trend", src, cs$id)
    rows[[id]] <- repro_row(id, cp_trend)
  } else {
    id <- sprintf("%s-%s", src, cs$id)
    rows[[id]] <- repro_row(id, cp_hr, note = note_d)
  }
}

result <- repro_finish(src, rows, diagnostics = do.call(rbind, diag))
