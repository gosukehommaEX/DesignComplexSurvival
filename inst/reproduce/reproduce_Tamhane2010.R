# Reproduction of Tamhane, Mehta and Liu (2010), Testing a primary and a
# secondary endpoint in a group sequential design, Biometrics 66:1174-1184:
# Table 1 (critical constants d of the secondary boundary and the values of
# Delta_1 that maximize the familywise error rate), the constants of
# Section 3.2, and the bounds of the example in Section 5.
#
# Two looks with equal information (gamma = sqrt(1/2)) and one-sided level
# 0.05. The trial stops at the look at which the primary hypothesis is
# rejected, and the secondary hypothesis is tested only at that look
# (gs_prob_mtp() with hierarchical = "stagewise"). Under H2 the familywise
# error rate is the probability of rejecting the secondary hypothesis, with
# means Delta_1 and Delta_1 / gamma of the primary statistics, means 0 of the
# secondary statistics and the correlations of equation (3).
#
# For rho = 1 the covariance matrix is singular. The secondary statistics
# then equal the primary statistics minus their means, and the familywise
# error rate P(X1 > max(c1, d1 + Delta_11)) + P(X1 <= c1,
# X2 > max(c2, d2 + Delta_12)) is computed from the primary statistics with
# pnorm() and gs_prob().
#
# The critical constant d solves max over Delta_1 >= 0 of the familywise
# error rate = 0.05. The maximum is located on a grid of Delta_1 from 0 to 8
# (step 0.1) and refined by optimize(). When the largest value is at the end
# of the grid, the familywise error rate increases in Delta_1 and its
# supremum is the limit P(Y1 > d1) as Delta_1 tends to infinity (printed as
# infinity in Table 1). The values of Delta_1 printed as infinity, and those
# printed as 0 (rho = 1, OF1-OF2 and PO1-PO2, attained at the end of the
# range), are not transcribed, because the perturbation check (published
# value times 1.1) cannot be applied to them.
#
# Table 1 is computed with the bounds printed in its notes: the primary
# bounds (1.678 sqrt(2), 1.678) and (1.876, 1.876), and z_0.05 = 1.645 in the
# ad hoc (AH) secondary boundary (d, z_0.05). The printed d is the root
# rounded up to 3 decimals, that is, the smallest d on a grid of 0.001 at
# which the maximum is at most 0.05 (criterion ceil_to:0.001). The printed
# Delta_1 maximizes the familywise error rate of the tabled boundary and is
# computed at the printed d. The constants themselves (Section 3.2) are
# computed exactly.
#
# Evidence for the registered explanations (explained.csv):
# - Tamhane2010-T1-OF1OF2-rho0.4-d and -delta: the printed d = 1.476 is a
#   misprint of 1.456; the maximizer at 1.456 (the root rounded up) is
#   compared with the printed Delta_1 of the cell, and the maximum at the
#   printed d is below 0.0495.
# - Tamhane2010-T1-OF1AH2-rho0.4-delta: 0.949 is a misprint of 1.949; the
#   computed value minus 1 is within 0.001 of 0.949 and the familywise
#   error rate at 1.949 is within 1e-7 of the maximum.
# - Tamhane2010-T1-flat: the maximum is flat; the package value is within
#   0.001 of the printed value and the familywise error rate at the printed
#   Delta_1 is within 1e-7 of the maximum.

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "Tamhane2010"
a <- 0.05
gam <- sqrt(0.5)
rows <- list()
r2 <- matrix(c(1, gam, gam, 1), 2)

prim <- list(OF = c(1.678 / gam, 1.678), PO = c(1.876, 1.876))
z_a <- 1.645
sec_bound <- function(type, d) {
  switch(type, OF = c(d / gam, d), PO = c(d, d), AH = c(d, z_a))
}

fwer <- function(d1, rho, cb, db) {
  m <- c(d1, d1 / gam)
  if (rho >= 1) {
    p_first <- stats::pnorm(max(cb[1], db[1] + m[1]) - m[1], lower.tail = FALSE)
    p_second <- gs_prob(c(-Inf, -Inf), c(cb[1], max(cb[2], db[2] + m[2])), m, r2)$prob_upper[2]
    return(p_first + p_second)
  }
  sigma <- rbind(cbind(r2, rho * r2), cbind(rho * r2, r2))
  gs_prob_mtp(list(cb, db), mean = c(m, 0, 0), sigma = sigma,
              hierarchical = "stagewise")$reject[["second"]]
}

max_fwer <- function(rho, cb, db) {
  grid <- seq(0, 8, by = 0.1)
  v <- vapply(grid, fwer, numeric(1), rho = rho, cb = cb, db = db)
  i <- which.max(v)
  if (i == length(grid)) {
    return(c(value = stats::pnorm(db[1], lower.tail = FALSE), delta = Inf))
  }
  o <- stats::optimize(fwer, grid[c(max(i - 1L, 1L), i + 1L)], maximum = TRUE,
                       rho = rho, cb = cb, db = db, tol = 1e-7)
  if (o$objective < v[i]) {
    return(c(value = v[i], delta = grid[i]))
  }
  c(value = o$objective, delta = o$maximum)
}

solve_d <- function(p, s, rho) {
  f <- function(d) max_fwer(rho, prim[[p]], sec_bound(s, d))[["value"]] - a
  stats::uniroot(f, c(1, 3), tol = 1e-8)$root
}
round_up3 <- function(x) ceiling(x * 1000 - 1e-9) / 1000

# Table 1
combos <- list(c("OF", "OF"), c("OF", "PO"), c("OF", "AH"), c("PO", "PO"),
               c("PO", "OF"), c("PO", "AH"))
pv <- repro_published(src)
pub_value <- function(id) pv$value[pv$id == id]
for (cmb in combos) {
  cb <- prim[[cmb[1]]]
  for (rho in c(0, 0.2, 0.4, 0.6, 0.8, 1)) {
    id <- sprintf("Tamhane2010-T1-%s1%s2-rho%.1f", cmb[1], cmb[2], rho)
    id_d <- paste0(id, "-d")
    id_delta <- paste0(id, "-delta")
    d_root <- solve_d(cmb[1], cmb[2], rho)
    if (id == "Tamhane2010-T1-OF1OF2-rho0.4") {
      ev <- max_fwer(rho, cb, sec_bound(cmb[2], round_up3(d_root)))[["delta"]]
      low <- max_fwer(rho, cb, sec_bound(cmb[2], pub_value(id_d)))[["value"]] < 0.0495
      rows[[id_d]] <- repro_row(id_d, d_root, evidence = ev, evidence_criterion = "abs:0.0005",
                                evidence_ref_id = id_delta, evidence_extra = low,
                                evidence_perturbed = ev)
    } else {
      rows[[id_d]] <- repro_row(id_d, d_root)
    }
    if (id_delta %in% pv$id) {
      db <- sec_bound(cmb[2], pub_value(id_d))
      mx <- max_fwer(rho, cb, db)
      if (id == "Tamhane2010-T1-OF1OF2-rho0.4") {
        ev <- max_fwer(rho, cb, sec_bound(cmb[2], round_up3(d_root)))[["delta"]]
        rows[[id_delta]] <- repro_row(id_delta, mx[["delta"]], evidence = ev,
                                      evidence_perturbed = ev)
      } else if (id == "Tamhane2010-T1-OF1AH2-rho0.4") {
        loss <- mx[["value"]] - fwer(1 + pub_value(id_delta), rho, cb, db)
        rows[[id_delta]] <- repro_row(id_delta, mx[["delta"]], evidence = mx[["delta"]] - 1,
                                      evidence_criterion = "abs:0.001",
                                      evidence_extra = loss < 1e-7,
                                      evidence_perturbed = mx[["delta"]] - 1,
                                      note = sprintf(paste("evidence: the computed value minus 1 (reading 1.949);",
                                                           "FWER at 1.949 below the maximum by %.1e"), loss))
      } else {
        loss <- mx[["value"]] - fwer(pub_value(id_delta), rho, cb, db)
        rows[[id_delta]] <- repro_row(id_delta, mx[["delta"]], evidence = mx[["delta"]],
                                      evidence_criterion = "abs:0.001",
                                      evidence_extra = loss < 1e-7,
                                      evidence_perturbed = mx[["delta"]],
                                      note = sprintf("FWER at the printed value below the maximum by %.1e",
                                                     loss))
      }
    }
  }
}

# Section 3.2: the classical constants (alpha = 0.05, printed rounded up)
# and the nominal level of the Pocock constant; evidence: the level of the
# printed constant 1.876
c_of <- repro_classical2("OF", a, 0.5)
c_po <- repro_classical2("PO", a, 0.5)
rows[["Tamhane2010-S32-OF-c"]] <- repro_row("Tamhane2010-S32-OF-c", c_of[2])
rows[["Tamhane2010-S32-PO-c"]] <- repro_row("Tamhane2010-S32-PO-c", c_po[1])
lev_printed <- stats::pnorm(1.876, lower.tail = FALSE)
rows[["Tamhane2010-S32-alpha-star"]] <- repro_row(
  "Tamhane2010-S32-alpha-star", stats::pnorm(c_po[1], lower.tail = FALSE),
  note = "nominal level of the exact Pocock constant",
  evidence = lev_printed, evidence_perturbed = lev_printed)

# Section 3.2, observation (4): PO primary boundary (1.876, 1.876) and OF
# secondary boundary (1.678 sqrt(2), 1.678), rho = 1
mx <- max_fwer(1, prim$PO, sec_bound("OF", 1.678))
rows[["Tamhane2010-S32-obs4-delta"]] <- repro_row("Tamhane2010-S32-obs4-delta", mx[["delta"]])
rows[["Tamhane2010-S32-obs4-fwer"]] <- repro_row("Tamhane2010-S32-obs4-fwer", mx[["value"]])

# Section 5: interim analysis at 1050 of 1400 subjects, one-sided 0.025;
# Lan-DeMets spending functions (5) and (6), and the classical bounds
b_of <- gs_bounds(c(0.75, 1), alpha = 0.025, sf_eff = "OF")$upper
b_po <- gs_bounds(c(0.75, 1), alpha = 0.025, sf_eff = "Pocock")$upper
rows[["Tamhane2010-S5-c1"]] <- repro_row("Tamhane2010-S5-c1", b_of[1])
rows[["Tamhane2010-S5-c2"]] <- repro_row("Tamhane2010-S5-c2", b_of[2])
rows[["Tamhane2010-S5-d1"]] <- repro_row("Tamhane2010-S5-d1", b_po[1])
rows[["Tamhane2010-S5-d2"]] <- repro_row("Tamhane2010-S5-d2", b_po[2])
cl_of <- repro_classical2("OF", 0.025, 0.75)
rows[["Tamhane2010-S5-classical-OF-c1"]] <- repro_row("Tamhane2010-S5-classical-OF-c1", cl_of[1])
rows[["Tamhane2010-S5-classical-OF-c2"]] <- repro_row("Tamhane2010-S5-classical-OF-c2", cl_of[2])
rows[["Tamhane2010-S5-classical-PO"]] <- repro_row("Tamhane2010-S5-classical-PO",
                                                   repro_classical2("PO", 0.025, 0.75)[1])

result <- repro_finish(src, rows)
