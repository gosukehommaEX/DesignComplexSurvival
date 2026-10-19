# Reproduction of Porcher, Levy and Chevret (2002), Sample size correction
# for treatment crossovers in randomized clinical trials with a survival
# endpoint, Controlled Clinical Trials 23:650-661: Tables 1 and 2 and the
# MAG90 example.
#
# Control group: death hazard lambda1 = -log(pi1) / T and relapse hazard
# lambda3 = -log(1 - F3(T)) / T; relapsed patients switch at once to the
# experimental treatment and then have the experimental hazard
# lambda2 = -log(pi2) / T. All patients enter at time 0 (approximated by
# uniform accrual over 1e-6 T) and are analysed at T. The inflation
# coefficient gamma is n / n0, where n is the sample size per group of the
# package (mean = "direct", variance = "null") and n0 that of Freedman's
# formula (equations (1) and (2)); gamma does not depend on the error rates.
#
# porcher_gamma() is the article's formula with K intervals (phi at the start
# of each interval, theta the mean hazard ratio of the interval), transcribed
# from Planning/phase3_feasibility/porcher_discrete.py. It is not part of the
# package.

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "Porcher2002"
za <- stats::qnorm(0.975)
zb <- stats::qnorm(0.9)

freedman_n0 <- function(p1, p2, za, zb) {
  th <- log(p1) / log(p2)
  (th + 1)^2 / (th - 1)^2 * (za + zb)^2 / (2 - p1 - p2)
}

porcher_gamma <- function(p1, p2, f3, k_int, t_end = 1) {
  l1 <- -log(p1) / t_end
  l2 <- -log(p2) / t_end
  l3 <- -log(1 - f3) / t_end
  s1 <- function(t) (l3 * exp(-l2 * t) + (l1 - l2) * exp(-(l1 + l3) * t)) / (l1 + l3 - l2)
  s2 <- function(t) exp(-l2 * t)
  k <- seq_len(k_int)
  a <- (k - 1) * t_end / k_int
  b <- k * t_end / k_int
  phi <- s1(a) / s2(a)
  theta <- log(s1(a) / s1(b)) / (l2 * t_end / k_int)
  rho <- (s1(a) - s1(b) + s2(a) - s2(b)) / ((1 - s1(t_end)) + (1 - s2(t_end)))
  zeta <- phi * theta / (1 + phi * theta) - phi / (1 + phi)
  eta <- phi / (1 + phi)^2
  n <- sum(rho * eta) / sum(rho * zeta)^2 / (2 - s1(t_end) - s2(t_end))
  n / freedman_n0(p1, p2, 0, 1)
}

porcher_scn <- function(p1, p2, f3, t_end = 1, eps = 1e-6) {
  l1 <- -log(p1) / t_end
  l2 <- -log(p2) / t_end
  l3 <- -log(1 - f3) / t_end
  ctl <- if (l3 > 0) ms_model(death = l1, prog = l3, death_prog = l2) else ms_model(death = l1)
  ms_scenario(ctl, ms_model(death = l2), accrual = repro_simult(t_end, eps))
}

gamma_pkg <- function(p1, p2, f3, t_end = 1, eps = 1e-6) {
  n <- repro_n_fixed(porcher_scn(p1, p2, f3, t_end, eps), t_end = t_end, alpha = 0.025,
                     power = 0.9) / 2
  n / freedman_n0(p1, p2, za, zb)
}

rows <- list()
helpers <- list()

# Table 1 (log-rank)
t1 <- rbind(data.frame(p2 = 0.50, p1 = 0.35, f3 = seq(0, 0.65, by = 0.05)),
            data.frame(p2 = 0.70, p1 = 0.60, f3 = seq(0, 0.40, by = 0.05)),
            data.frame(p2 = 0.70, p1 = 0.65, f3 = seq(0, 0.35, by = 0.05)))
for (i in seq_len(nrow(t1))) {
  p1 <- t1$p1[i]
  p2 <- t1$p2[i]
  f3 <- t1$f3[i]
  id <- sprintf("Porcher2002-T1-%.2f-%.2f-%.2f", p2, p1, f3)
  g <- gamma_pkg(p1, p2, f3)
  g50 <- porcher_gamma(p1, p2, f3, 50)
  g1e4 <- porcher_gamma(p1, p2, f3, 10000)
  helpers[[length(helpers) + 1L]] <- data.frame(case = id, p1 = p1, p2 = p2, f3 = f3,
                                                k_int = c(50, 10000), gamma = c(g50, g1e4))
  note <- sprintf("formula with K = 50: %.4f; K = 10000: %.5f", g50, g1e4)
  if (id == "Porcher2002-T1-0.50-0.35-0.20") {
    rows[[id]] <- repro_row(id, g, note = note, evidence = g50,
                            evidence_extra = abs(g / g1e4 - 1) <= 1e-4)
  } else {
    rows[[id]] <- repro_row(id, g, note = note)
  }
}

# MAG90 example: pi1 0.30, pi2 0.50, T 60 months, two-sided 5%, power 80%,
# n0 = 91 per group (Freedman's formula rounded up)
za2 <- stats::qnorm(0.975)
zb2 <- stats::qnorm(0.8)
n0 <- ceiling(freedman_n0(0.30, 0.50, za2, zb2) - 1e-9)
for (f3 in c(0.5, 0.6, 0.7, 0.8)) {
  scn <- porcher_scn(0.30, 0.50, f3, t_end = 60)
  n_arm <- repro_n_fixed(scn, t_end = 60, alpha = 0.025, power = 0.8) / 2
  g <- n_arm / freedman_n0(0.30, 0.50, za2, zb2)
  id <- sprintf("Porcher2002-MAG90-F%.2f-n", f3)
  # evidence: gamma rounded to 2 decimals times n0, rounded up
  rows[[id]] <- repro_row(id, n_arm, evidence = ceiling(round(g, 2) * n0 - 1e-9),
                          note = sprintf("gamma %.4f; n0 %d", g, n0))
  id <- sprintf("Porcher2002-MAG90-F%.2f-power", f3)
  rows[[id]] <- repro_row(id, 100 * repro_power_fixed(scn, 2 * n0, 60, 0.025),
                          note = sprintf("power with %d per group", n0))
}

# Table 2 (INFO): lambda3 = ratio * lambda1, with the printed pi1
t2 <- data.frame(n0 = c(100, 100, 500, 500), p1 = c(0.07, 0.37, 0.13, 0.50),
                 p2 = c(0.20, 0.60, 0.20, 0.60), lab = c("0.20", "0.60", "0.20", "0.6"))
ratios <- c(0.1, 0.2, 0.33, 0.5, 1)
for (i in seq_len(nrow(t2))) {
  p1 <- t2$p1[i]
  p2 <- t2$p2[i]
  # pi1 implied by n0 (power 90%, two-sided 5%)
  p1x <- stats::uniroot(function(x) freedman_n0(x, p2, za, zb) - t2$n0[i],
                        c(1e-4, p2 - 1e-4), tol = 1e-12)$root
  for (r in ratios) {
    id <- sprintf("Porcher2002-T2-n%d-pi2%s-r%s", t2$n0[i], t2$lab[i], format(r))
    g <- gamma_pkg(p1, p2, 1 - p1^r)
    gx <- gamma_pkg(p1x, p2, 1 - p1x^r)
    rows[[id]] <- repro_row(id, g, note = sprintf("pi1 from n0 (%.4f): %.4f", p1x, gx))
  }
}

# Diagnostics of the approximation of simultaneous entry
diag <- list()
scn0 <- porcher_scn(0.35, 0.50, 0)
ev0 <- ms_events(scn0, n = 1000, calendar_time = 1)$events
diag[[1]] <- repro_diag("Porcher2002-simultaneous-events-F0", ev0,
                        500 * (1 - 0.35) + 500 * (1 - 0.50), 1e-5,
                        "expected events at T with accrual over 1e-6 T versus the closed form")
g6 <- gamma_pkg(0.35, 0.50, 0.5, eps = 1e-6)
g5 <- gamma_pkg(0.35, 0.50, 0.5, eps = 1e-5)
g7 <- gamma_pkg(0.35, 0.50, 0.5, eps = 1e-7)
diag[[2]] <- repro_diag("Porcher2002-simultaneous-gamma-eps1e-5", g6, g5, 1e-4,
                        "gamma (pi2 0.50, pi1 0.35, F3 0.50) with accrual over 1e-6 T versus 1e-5 T")
diag[[3]] <- repro_diag("Porcher2002-simultaneous-gamma-eps1e-7", g6, g7, 1e-4,
                        "gamma (pi2 0.50, pi1 0.35, F3 0.50) with accrual over 1e-6 T versus 1e-7 T")

result <- repro_finish(src, rows, diagnostics = do.call(rbind, diag),
                       helpers = do.call(rbind, helpers))
