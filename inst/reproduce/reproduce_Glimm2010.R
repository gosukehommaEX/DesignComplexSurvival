# Reproduction of Glimm, Maurer and Bretz (2010), Hierarchical testing of
# multiple endpoints in group-sequential trials, Statistics in Medicine
# 29:219-228: the Pocock bounds of Section 2, the sharp upper bound and the
# Pocock constant of Section 3.2, the parameters of the Kim-DeMets spending
# function in Section 4, and the numerical example of Section 5.
#
# Section 5: one interim analysis at half of the information, one-sided
# level 0.025, classical O'Brien-Fleming bounds for the primary endpoint
# and classical Pocock bounds for the secondary endpoint. The drifts are
# delta_p = z_0.975 + z_0.92 and delta_s = z_0.975 + z_0.86 (powers of 92
# and 86 per cent for a test without interim analysis), and the
# statistics have the distribution (7) with correlation 0 or 0.8 between
# the endpoints. Strategy (a) is gs_prob_mtp() with hierarchical =
# "stagewise", strategy (b) with "overall", and strategy (d) (coequal) with
# mtp = "bonferroni" (each hypothesis tested at its own bounds until both
# are rejected). Strategy (c) (partially hierarchical) is not provided by
# the package and is not transcribed; neither is the probability of
# "losing" a significant interim result (0.09 per cent), which is not a
# rejection probability, nor the lower bounds 82 and 73 per cent computed
# from the normal distribution function alone.
#
# The text prints the final secondary bound as 2.17 in P(Zs,2 >= 2.17); the
# probabilities 80.6 and 80.7 per cent follow from the Pocock bound 2.178.

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "Glimm2010"
rows <- list()
r2 <- function(t1) matrix(c(1, sqrt(t1), sqrt(t1), 1), 2)

# Section 2: Pocock bounds for alpha = 0.025 and t = 0.5
po <- repro_classical2("PO", 0.025, 0.5)
sp <- gs_prob(c(-Inf, -Inf), po, 0, r2(0.5))$prob_upper
rows[["Glimm2010-S2-PO-u"]] <- repro_row("Glimm2010-S2-PO-u", po[1])
rows[["Glimm2010-S2-PO-nominal"]] <- repro_row("Glimm2010-S2-PO-nominal",
                                               stats::pnorm(po[1], lower.tail = FALSE))
rows[["Glimm2010-S2-PO-spent1"]] <- repro_row("Glimm2010-S2-PO-spent1", sp[1])
rows[["Glimm2010-S2-PO-spent2"]] <- repro_row("Glimm2010-S2-PO-spent2", sp[2])

# Section 3.2: sharp upper bound (9) of the naive strategy (us,1 = us,2 =
# z_0.95, alpha = 0.05, t = 0.5), and the constant bound with (9) = 0.05.
# Evidence (see explained.csv): agreement within one unit in the last digit
# printed, and the probability (9) at the printed bound 1.8755 within 0.001
# of 0.05; the printed nominal level is that of the printed bound.
z <- stats::qnorm(0.95)
prob9 <- function(u) sum(gs_prob(c(-Inf, -Inf), c(u, u), 0, r2(0.5))$prob_upper)
rows[["Glimm2010-S32-naive"]] <- repro_row("Glimm2010-S32-naive", prob9(z))
u05 <- repro_classical2("PO", 0.05, 0.5)[1]
ok9 <- abs(prob9(1.8755) - 0.05) <= 0.001
rows[["Glimm2010-S32-us"]] <- repro_row(
  "Glimm2010-S32-us", u05, evidence = u05, evidence_criterion = "abs:0.0001",
  evidence_extra = ok9, evidence_perturbed = u05,
  note = sprintf("probability (9) at the printed bound: %.6f", prob9(1.8755)))
lev <- stats::pnorm(1.8755, lower.tail = FALSE)
rows[["Glimm2010-S32-us-level"]] <- repro_row(
  "Glimm2010-S32-us-level", stats::pnorm(u05, lower.tail = FALSE),
  evidence = lev, evidence_perturbed = lev)

# Section 4: parameters of the Kim-DeMets spending function alpha t^gamma
# whose bounds at t = 0.5 match the Pocock and the O'Brien-Fleming bounds
# (alpha = 0.025; the alpha spent at the interim analysis is matched)
match_power <- function(type) {
  b <- repro_classical2(type, 0.025, 0.5)
  a1 <- stats::pnorm(b[1], lower.tail = FALSE)
  stats::uniroot(function(g) gs_spend(0.5, 0.025, "power", g) - a1, c(0.05, 20),
                 tol = 1e-12)$root
}
rows[["Glimm2010-S4-gamma-PO"]] <- repro_row("Glimm2010-S4-gamma-PO", match_power("PO"))
rows[["Glimm2010-S4-gamma-OF"]] <- repro_row("Glimm2010-S4-gamma-OF", match_power("OF"))

# Section 5
of <- repro_classical2("OF", 0.025, 0.5)
rows[["Glimm2010-S5-up1"]] <- repro_row("Glimm2010-S5-up1", of[1])
rows[["Glimm2010-S5-up2"]] <- repro_row("Glimm2010-S5-up2", of[2])
rows[["Glimm2010-S5-us"]] <- repro_row("Glimm2010-S5-us", po[1])
dp <- stats::qnorm(0.975) + stats::qnorm(0.92)
ds <- stats::qnorm(0.975) + stats::qnorm(0.86)
ex <- function(rho, mtp = "fixed_sequence", hierarchical = "overall") {
  r <- r2(0.5)
  sigma <- rbind(cbind(r, rho * r), cbind(rho * r, r))
  gs_prob_mtp(list(of, po), mean = c(sqrt(0.5) * dp, dp, sqrt(0.5) * ds, ds),
              sigma = sigma, mtp = mtp, hierarchical = hierarchical)
}
for (rho in c(0, 0.8)) {
  tag <- if (rho == 0) "rho0" else "rho0.8"
  sa <- ex(rho, hierarchical = "stagewise")
  p <- sa$looks$prob_upper
  if (rho == 0) {
    rows[["Glimm2010-S5-Pp1"]] <- repro_row("Glimm2010-S5-Pp1", 100 * p[1])
  }
  rows[[paste0("Glimm2010-S5-cond1-", tag)]] <- repro_row(paste0("Glimm2010-S5-cond1-", tag),
                                                         100 * p[3] / p[1])
  rows[[paste0("Glimm2010-S5-cond2-", tag)]] <- repro_row(paste0("Glimm2010-S5-cond2-", tag),
                                                         100 * p[4] / p[2])
  # evidence for the row with correlation 0.8 (see explained.csv)
  rows[[paste0("Glimm2010-S5-a-", tag)]] <- repro_row(
    paste0("Glimm2010-S5-a-", tag), 100 * sa$reject[["second"]],
    evidence = 100 * sa$reject[["second"]], evidence_criterion = "abs:0.1",
    evidence_perturbed = 100 * sa$reject[["second"]])
  rows[[paste0("Glimm2010-S5-b-", tag)]] <- repro_row(
    paste0("Glimm2010-S5-b-", tag), 100 * ex(rho)$reject[["second"]])
}
rows[["Glimm2010-S5-d"]] <- repro_row("Glimm2010-S5-d",
                                      100 * ex(0, mtp = "bonferroni")$reject[["second"]])
rows[["Glimm2010-S5-us-t0.8"]] <- repro_row("Glimm2010-S5-us-t0.8",
                                            repro_classical2("PO", 0.025, 0.8)[1])

result <- repro_finish(src, rows)
