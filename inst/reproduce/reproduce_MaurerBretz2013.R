# Reproduction of Maurer and Bretz (2013), Multiple testing in group
# sequential trials using graphical approaches, Statistics in
# Biopharmaceutical Research 5:311-320: the nominal significance levels of
# Table 1 and of the text of Section 4, and the repeated p-values of Table 2
# and of the text of Section 4 (computed in the article with ADDPLAN 6.0.1).
#
# Case study of Section 4: O'Brien-Fleming type spending function at the
# information fractions 1/3, 2/3 and 1. The nominal level of a hypothesis at
# level gamma and look t is the nominal p-value of the efficacy bound of
# gs_bounds() at level gamma (the bounds after recycling of ms_power() are
# computed in the same way). The repeated p-value at look t is the level
# gamma at which the nominal level equals the nominal p-value at look t
# (Section 3.3); it is found by root finding over gamma. The levels of the
# case study are 0.0125 (alpha / 2, initial level of a primary hypothesis),
# 0.01875 and 0.00625 (after the rejection of H1) and 0.025 (H4 alone).
#
# Table 2 prints the nominal p-values with one or two significant digits.
# The repeated p-values computed from them differ from four printed
# repeated p-values in the last digit; the evidence (see explained.csv) is
# the range of the repeated p-value over the nominal p-values that round to
# the printed one. The probabilities of the procedure (power, familywise
# error rate) are not given in the article.

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "MaurerBretz2013"
rows <- list()
info <- c(1, 2, 3) / 3

nominal <- function(gamma, look) {
  gs_bounds(info, alpha = gamma, sf_eff = "OF")$nominal_p[look]
}

# Table 1 and the text of p. 317
lev <- list(
  "MaurerBretz2013-T1-a1-0.0125" = c(0.0125, 1),
  "MaurerBretz2013-T1-a2-0.0125" = c(0.0125, 2),
  "MaurerBretz2013-T1-a2-0.01875" = c(0.01875, 2),
  "MaurerBretz2013-T1-a2-0.00625" = c(0.00625, 2),
  "MaurerBretz2013-T1-a3-0.025" = c(0.025, 3),
  "MaurerBretz2013-S4-a1-0.025" = c(0.025, 1))
for (id in names(lev)) {
  rows[[id]] <- repro_row(id, nominal(lev[[id]][1], lev[[id]][2]))
}

# Table 2 and the text of p. 318: repeated p-values
repeated_p <- function(p, look) {
  stats::uniroot(function(g) nominal(g, look) - p, c(1e-6, 0.9), tol = 1e-12)$root
}
rp <- list(
  "MaurerBretz2013-T2-r1-1" = c("0.0062", 1),
  "MaurerBretz2013-T2-r2-1" = c("0.017", 1),
  "MaurerBretz2013-T2-r3-1" = c("0.009", 1),
  "MaurerBretz2013-T2-r4-1" = c("0.13", 1),
  "MaurerBretz2013-T2-r1-2" = c("0.0002", 2),
  "MaurerBretz2013-T2-r2-2" = c("0.0035", 2),
  "MaurerBretz2013-T2-r3-2" = c("0.0020", 2),
  "MaurerBretz2013-T2-r4-2" = c("0.06", 2),
  "MaurerBretz2013-S4-r2-2" = c("0.0035", 2))
pv <- repro_published(src)
for (id in names(rp)) {
  p <- as.numeric(rp[[id]][1])
  look <- as.integer(rp[[id]][2])
  # half a unit of the last digit of the printed nominal p-value
  h <- 0.5 * 10^-(nchar(rp[[id]][1]) - 2L)
  rng <- c(repeated_p(p - h, look), repeated_p(p + h, look))
  pub <- pv$value[pv$id == id]
  rows[[id]] <- repro_row(
    id, repeated_p(p, look),
    evidence = min(max(pub, rng[1L]), rng[2L]),
    evidence_perturbed = min(max(1.1 * pub, rng[1L]), rng[2L]),
    note = sprintf(paste("nominal p-value %s; repeated p-value from %.6f to %.6f for nominal",
                         "p-values from %s to %s"),
                   rp[[id]][1], rng[1L], rng[2L], format(p - h), format(p + h)))
}

result <- repro_finish(src, rows)
