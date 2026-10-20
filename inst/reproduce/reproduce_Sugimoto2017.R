# Reproduction of Sugimoto et al. (2017), Statistics in Medicine 36:1363-1382,
# doi:10.1002/sim.7225: Tables III and VI, one fatal outcome (semi-competing
# risks) in the composite setting (PFS and OS), rows with rho = 0.
#
# With rho = 0 the times to progression (TTP) and to death (OS) of a patient
# are independent and exponential, so that PFS is exponential with the sum of
# the two hazards and the hazard of death does not change at progression (the
# model of ms_fleischer()). The hazards follow from the survival rates at
# tau = tau_a + tau_f = 5 in the control group (S_TTP(tau) = 0.3, and
# S_OS(tau) = 0.4 in Table III and 0.5 in Table VI) and the hazard ratios
# control / test of 1.7 for TTP and 1.3, 1.5 or 1.7 for OS. The censoring time
# U(0, tau_a) + tau_f corresponds to uniform accrual over tau_a = 2 and one
# analysis at tau = 5. No dropout, equal allocation, one-sided level 0.025,
# power 0.8. The article rounds the sample sizes up to even numbers.
#
# Table III (co-primary endpoints): N_cp is the sample size with which both
# hypotheses, each tested at level 0.025, are rejected with probability 0.8.
# It is solved with ms_size(target = "both") under the fixed sequence: with
# one analysis, rejecting both hypotheses under the fixed sequence is
# rejecting both at level 0.025. Table VI (multiple primary endpoints,
# weighted Bonferroni procedure): N_mp is the smallest sample size over the
# weights omega = 0, 0.05, ..., 1 of PFS (omega_0); omega = 1 and omega = 0
# are the designs for PFS alone and for OS alone. The search over omega uses
# the mean and covariance of lr_cov(), which do not depend on n apart from
# the factor sqrt(n) of the mean. As diagnostics, N_cp, N_PFS, N_OS and N_mp
# at omega = 0.95 are solved with ms_size() and compared with the solutions
# from lr_cov().
#
# The package uses variance = "alternative" (Luo et al. 2019). Formula (2) of
# the article takes the variances and the covariance of the scores with the
# weights H_k fixed (explained.csv, Sugimoto2017-variance); it is computed
# below as the evidence, and its values are written to
# results/helpers_Sugimoto2017.csv.
#
# Not reproduced: the rows with rho > 0 (dependence by the Clayton and Gumbel
# copulas, which the multistate models of the package do not represent), the
# non-composite setting (TTP is not an endpoint of the package) and Table V
# (ICON7), whose hazard ratios are printed rounded (1.31 and 1.23); with these
# values formula (2) does not give the printed sample sizes. The Gumbel rows
# with rho = 0 print the same values from formulas (8) and (11) as the
# Clayton rows and are not repeated.

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "Sugimoto2017"
a <- 0.025
pw <- 0.8
tau_a <- 2
tau_f <- 3
tau <- tau_a + tau_f
looks <- ms_looks("calendar", time = tau)
w_grid <- round(seq(0.05, 0.95, by = 0.05), 2)
pv <- repro_published(src)
pub <- function(id) pv$value[pv$id == id]
rows <- list()
diag <- list()
helpers <- list()

# Hazards of progression (prog) and death, control and test group
sugi_rates <- function(s_os, hr_os) {
  list(prog = -log(0.3) / tau * c(1, 1 / 1.7),
       death = -log(s_os) / tau * c(1, 1 / hr_os))
}

sugi_scenario <- function(r) {
  ms_scenario(ms_model(prog = r$prog[1], death = r$death[1], death_prog = r$death[1]),
              ms_model(prog = r$prog[2], death = r$death[2], death_prog = r$death[2]),
              accrual = repro_unif(tau_a))
}

# Z = (Z_PFS, Z_OS) is normal with mean sqrt(n) d and covariance s; the
# hypothesis of endpoint k is rejected when Z_k exceeds z[k].
sugi_p_one <- function(n, d, s, k, z) {
  stats::pnorm((sqrt(n) * d[k] - z) / sqrt(s[k, k]))
}

# P(Z_PFS > z[1], Z_OS > z[2]) as an integral over the standardized Z_PFS
# (the normal density is negligible outside (-10, 10))
sugi_p_both <- function(n, d, s, z) {
  m <- sqrt(n) * d
  sd <- sqrt(diag(s))
  rho <- s[1, 2] / (sd[1] * sd[2])
  f <- function(u) {
    stats::dnorm(u) * stats::pnorm((m[2] + rho * sd[2] * u - z[2]) / (sd[2] * sqrt(1 - rho^2)))
  }
  lo <- max((z[1] - m[1]) / sd[1], -10)
  if (lo >= 10) {
    return(0)
  }
  stats::integrate(f, lo, 10, rel.tol = 1e-10, abs.tol = 1e-14)$value
}

sugi_p_any <- function(n, d, s, z) {
  sugi_p_one(n, d, s, 1L, z[1]) + sugi_p_one(n, d, s, 2L, z[2]) - sugi_p_both(n, d, s, z)
}

# Sample size with probability pw (the probabilities increase with n)
sugi_solve <- function(f) {
  stats::uniroot(function(n) f(n) - pw, c(10, 1e5), tol = 1e-10)$root
}

# N_PFS, N_OS and N_cp, and with weights, N_mp over omega = 0, weights, 1
sugi_sizes <- function(ds, weights = NULL) {
  z <- stats::qnorm(1 - a)
  out <- list(n_pfs = sugi_solve(function(n) sugi_p_one(n, ds$d, ds$s, 1L, z)),
              n_os = sugi_solve(function(n) sugi_p_one(n, ds$d, ds$s, 2L, z)),
              n_cp = sugi_solve(function(n) sugi_p_both(n, ds$d, ds$s, c(z, z))))
  if (!is.null(weights)) {
    n_w <- vapply(weights, function(w) {
      zw <- stats::qnorm(1 - c(w, 1 - w) * a)
      sugi_solve(function(n) sugi_p_any(n, ds$d, ds$s, zw))
    }, numeric(1))
    om <- c(0, weights, 1)
    n_all <- c(out$n_os, n_w, out$n_pfs)
    i <- which.min(n_all)
    out$n_w <- n_w
    out$omega0 <- om[i]
    out$n_mp <- n_all[i]
  }
  out
}

# Package: mean per sqrt(n) and covariance of (Z_PFS, Z_OS) at tau
sugi_package_design <- function(scn) {
  n0 <- 1000
  jd <- lr_cov(scn, n = n0, calendar_time = tau, endpoint = c("pfs", "os"),
               variance = "alternative")
  i <- match(c("pfs", "os"), jd$endpoint)
  list(d = jd$mean[i] / sqrt(n0), s = jd$sigma[i, i])
}

# Package: sample size (unrounded) and power of the design with one analysis
size_pkg <- function(scn, endpoints, target = "first", ...) {
  tst <- ms_test(alpha = a, endpoints = endpoints, variance = "alternative", ...)
  ms_size(scn, looks = looks, test = tst, power = pw, solve_for = "n",
          target = target)$solution
}

power_pkg <- function(scn, n, endpoints, ...) {
  tst <- ms_test(alpha = a, endpoints = endpoints, variance = "alternative", ...)
  ms_power(scn, n = n, looks = looks, test = tst)$power
}

# Formula (2) of the article for independent exponential TTP and OS: means mu,
# null variances v0 and variances v of the scores of PFS and OS, and their
# covariance v12, with the weights H_k = G r_0 r_1 S_0 S_1 / (r_0 S_0 + r_1 S_1)
# fixed (S_j is the survival function of the endpoint in group j, and G that
# of the censoring time, here uniform on (tau_f, tau)). The hazard of death
# does not change at progression, so the counting process martingales of PFS
# and OS of a patient have the predictable covariation int Y_PFS b_j dt (a
# death before progression is a jump of both), and
# v12 = sum_j int H_PFS H_OS b_j / (r_j G S_OS,j) dt. mu is taken positive
# for a benefit of the test group.
sugi_formula <- function(r, alloc = c(0.5, 0.5)) {
  G <- function(t) ifelse(t <= tau_f, 1, pmax(tau - t, 0) / tau_a)
  lam <- list(pfs = r$prog + r$death, os = r$death)
  S <- function(t, j, ep) exp(-lam[[ep]][j] * t)
  hw <- function(t, ep) {
    s0 <- S(t, 1L, ep)
    s1 <- S(t, 2L, ep)
    alloc[1] * alloc[2] * s0 * s1 / (alloc[1] * s0 + alloc[2] * s1)
  }
  int <- function(f) {
    stats::integrate(f, 0, tau_f, rel.tol = 1e-10, abs.tol = 1e-15)$value +
      stats::integrate(f, tau_f, tau, rel.tol = 1e-10, abs.tol = 1e-15)$value
  }
  out <- list()
  for (ep in c("pfs", "os")) {
    l <- lam[[ep]]
    out[[paste0("mu_", ep)]] <- int(function(t) G(t) * hw(t, ep) * (l[1] - l[2]))
    out[[paste0("v0_", ep)]] <- int(function(t) {
      G(t) * hw(t, ep)^2 * (l[1] / (alloc[2] * S(t, 2L, ep)) + l[2] / (alloc[1] * S(t, 1L, ep)))
    })
    out[[paste0("v_", ep)]] <- int(function(t) {
      G(t) * hw(t, ep)^2 * (l[1] / (alloc[1] * S(t, 1L, ep)) + l[2] / (alloc[2] * S(t, 2L, ep)))
    })
  }
  out$v12 <- sum(vapply(1:2, function(j) {
    int(function(t) G(t) * hw(t, "pfs") * hw(t, "os") * r$death[j] / (alloc[j] * S(t, j, "os")))
  }, numeric(1)))
  out
}

sugi_formula_design <- function(f) {
  s12 <- f$v12 / sqrt(f$v0_pfs * f$v0_os)
  list(d = c(f$mu_pfs / sqrt(f$v0_pfs), f$mu_os / sqrt(f$v0_os)),
       s = matrix(c(f$v_pfs / f$v0_pfs, s12, s12, f$v_os / f$v0_os), 2))
}

sugi_helper_row <- function(case, s_os, hr_os, r, f, ev) {
  data.frame(case = case, s_os = s_os, hr_os = hr_os, prog0 = r$prog[1], prog1 = r$prog[2],
             death0 = r$death[1], death1 = r$death[2], mu_pfs = f$mu_pfs, v0_pfs = f$v0_pfs,
             v_pfs = f$v_pfs, mu_os = f$mu_os, v0_os = f$v0_os, v_os = f$v_os, v12 = f$v12,
             n_pfs = ev$n_pfs, n_os = ev$n_os, n_cp = ev$n_cp,
             omega0 = if (is.null(ev$omega0)) NA_real_ else ev$omega0,
             n_mp = if (is.null(ev$n_mp)) NA_real_ else ev$n_mp, stringsAsFactors = FALSE)
}

dnote <- "ms_size() against the solution from lr_cov()"
hr_os <- c(1.3, 1.5, 1.7)

# Table III: co-primary endpoints, S_OS(tau) = 0.4
for (h in hr_os) {
  id <- function(q) sprintf("%s-T3-OS%.1f-%s", src, h, q)
  r <- sugi_rates(0.4, h)
  scn <- sugi_scenario(r)
  f <- sugi_formula(r)
  ev <- sugi_sizes(sugi_formula_design(f))
  pk <- sugi_sizes(sugi_package_design(scn))
  n_cp <- size_pkg(scn, c("pfs", "os"), target = "both", mtp = "fixed_sequence")
  n_pfs <- size_pkg(scn, "pfs")
  n_os <- size_pkg(scn, "os")
  diag[[length(diag) + 1L]] <- repro_diag(id("Ncp-lr_cov"), n_cp, pk$n_cp, 1e-6, dnote)
  diag[[length(diag) + 1L]] <- repro_diag(id("NPFS-lr_cov"), n_pfs, pk$n_pfs, 1e-6, dnote)
  diag[[length(diag) + 1L]] <- repro_diag(id("NOS-lr_cov"), n_os, pk$n_os, 1e-6, dnote)
  n_pub <- pub(id("Ncp"))
  p_cp <- power_pkg(scn, n_pub, c("pfs", "os"), mtp = "fixed_sequence")[["both"]]
  rows[[id("Nsim")]] <- repro_row(id("Nsim"), n_cp, note = "N_cp of the package (unrounded)")
  rows[[id("Ncp")]] <- repro_row(id("Ncp"), n_cp, evidence = ev$n_cp)
  rows[[id("Pcp")]] <- repro_row(id("Pcp"), 100 * p_cp,
                                 note = sprintf("power to reject both at the printed N_cp = %g", n_pub))
  rows[[id("NPFS")]] <- repro_row(id("NPFS"), n_pfs, evidence = ev$n_pfs)
  rows[[id("NOS")]] <- repro_row(id("NOS"), n_os, evidence = ev$n_os)
  helpers[[length(helpers) + 1L]] <- sugi_helper_row(sprintf("T3-OS%.1f", h), 0.4, h, r, f, ev)
}

# Table VI: multiple primary endpoints (weighted Bonferroni), S_OS(tau) = 0.5
for (h in hr_os) {
  id <- function(q) sprintf("%s-T6-OS%.1f-%s", src, h, q)
  r <- sugi_rates(0.5, h)
  scn <- sugi_scenario(r)
  f <- sugi_formula(r)
  ev <- sugi_sizes(sugi_formula_design(f), weights = w_grid)
  pk <- sugi_sizes(sugi_package_design(scn), weights = w_grid)
  n_pfs <- size_pkg(scn, "pfs")
  n_os <- size_pkg(scn, "os")
  n_095 <- size_pkg(scn, c("pfs", "os"), target = "any", mtp = "bonferroni",
                    weights = c(0.95, 0.05))
  diag[[length(diag) + 1L]] <- repro_diag(id("NPFS-lr_cov"), n_pfs, pk$n_pfs, 1e-6, dnote)
  diag[[length(diag) + 1L]] <- repro_diag(id("NOS-lr_cov"), n_os, pk$n_os, 1e-6, dnote)
  diag[[length(diag) + 1L]] <- repro_diag(id("Nmp-w0.95-lr_cov"), n_095,
                                          pk$n_w[abs(w_grid - 0.95) < 1e-9], 1e-6,
                                          paste(dnote, "(Bonferroni, omega = 0.95, target any)"))
  w0 <- pk$omega0
  n_mp <- if (w0 == 1) {
    n_pfs
  } else if (w0 == 0) {
    n_os
  } else {
    size_pkg(scn, c("pfs", "os"), target = "any", mtp = "bonferroni", weights = c(w0, 1 - w0))
  }
  n_pub <- pub(id("Nmp"))
  p_mp <- if (w0 == 1) {
    power_pkg(scn, n_pub, "pfs")
  } else if (w0 == 0) {
    power_pkg(scn, n_pub, "os")
  } else {
    power_pkg(scn, n_pub, c("pfs", "os"), mtp = "bonferroni", weights = c(w0, 1 - w0))[["any"]]
  }
  rows[[id("w0")]] <- repro_row(
    id("w0"), w0,
    note = sprintf("package N_mp %.2f at omega_0, %.2f at omega = 0.95; formula (2): omega_0 = %.2f",
                   n_mp, pk$n_w[abs(w_grid - 0.95) < 1e-9], ev$omega0))
  rows[[id("Nmp")]] <- repro_row(id("Nmp"), n_mp, evidence = ev$n_mp,
                                 note = sprintf("omega_0 = %.2f", w0))
  rows[[id("Pmp")]] <- repro_row(id("Pmp"), 100 * p_mp,
                                 note = sprintf("power at the printed N_mp = %g with omega_0 = %.2f", n_pub, w0))
  rows[[id("NPFS")]] <- repro_row(id("NPFS"), n_pfs, evidence = ev$n_pfs)
  rows[[id("NOS")]] <- repro_row(id("NOS"), n_os, evidence = ev$n_os)
  helpers[[length(helpers) + 1L]] <- sugi_helper_row(sprintf("T6-OS%.1f", h), 0.5, h, r, f, ev)
}

result <- repro_finish(src, rows, diagnostics = do.call(rbind, diag),
                       helpers = do.call(rbind, helpers))
