# Reproduction of Kim and Tsiatis (1990), Study duration for clinical trials
# with survival response and early stopping rule, Biometrics 46:81-92: the
# ECOG example of Section 5 (Tables 1-6 and the values in Sections 5 and 6).
#
# Design: one-sided level 0.05 and power 0.90 for the log hazard ratio
# -0.41, exponential survival with the failure rate 1.02 per year in the
# control group, 100 patients per year accrued uniformly for 2.5 years (250
# patients), no loss to follow-up, and five analyses at the information
# times t(1) = (.2, .4, .6, .8, 1), t(2) = (.3, .6, .8, .9, 1) or
# t(3) = (.1, .2, .3, .6, 1), with the use function of O'Brien-Fleming type
# (sf_eff = "OF") or of Pocock type (sf_eff = "Pocock"). The drift of the
# article, eta = |theta| sqrt(E1(Sa + Sf) / 4), is the mean of the final Z
# with mean = "schoenfeld" and equal allocation.
#
# Convention of the tables (Section 4): the follow-up Sf is the printed one,
# and the final analysis is at S = Sa + Sf under both hypotheses, with
# F1 = E1(S) failures under H1 and F0 = E0(S) under H0. The k-th analysis is
# at the real time s_k with E0(s_k) = t_k F0, and the expected values under
# both hypotheses use these times: the real time, the patients a(s_k), the
# failures t_k F0 under H0 and t_k F1 under H1, and the probabilities of
# stopping under H0 and H1 (at the last analysis one minus the probability
# of stopping earlier). The bounds and the means of the statistics come from
# ms_power() with followup = Sf, the probabilities from gs_prob() with the
# correlations sqrt(t_j / t_k) of the Brownian motion of the article (the
# information is proportional to the number of failures, V(s) = E(s) / 4),
# and the real times, failures and patients from ms_event_time() and
# ms_events() under ms_null(). ms_power() itself takes the correlations
# under the design scenario from the information (variance = "null"), whose
# fractions differ from t here (0.205 at t = 0.2 for the design of Table 1),
# so that its power differs from that of the article's model in the fourth
# decimal; it is given in the notes. The expected values of the package
# (ms_power()$expected), which use under each hypothesis its own analysis
# times (under H0 at the numbers of failures of the design scenario) and the
# correlations of the information, are given as INFO rows next to the
# printed values (ids ending in -pkg).
#
# The numbers of failures and patients are printed rounded up (criterion
# ceil). The follow-up is judged by the power at the printed follow-up
# (criterion power_at_time; Sections 3 and 5 choose Sf for the power 0.90).

if (!exists("repro_finish", mode = "function")) {
  source(file.path("inst", "reproduce", "helpers.R"))
}

src <- "KimTsiatis1990"
rows <- list()
n_kt <- 250
scn_kt <- ms_scenario(ms_model(death = 1.02), ms_model(death = exp(-0.41) * 1.02),
                      accrual = repro_unif(2.5))
scn0_kt <- ms_null(scn_kt)
times_kt <- list(t1 = c(0.2, 0.4, 0.6, 0.8, 1), t2 = c(0.3, 0.6, 0.8, 0.9, 1),
                 t3 = c(0.1, 0.2, 0.3, 0.6, 1))
followup_kt <- list(OF = c(t1 = 1.22, t2 = 1.29, t3 = 1.10),
                    Pocock = c(t1 = 3.87, t2 = 3.92, t3 = 2.77))
table_kt <- c(OF.t1 = "T1", Pocock.t1 = "T2", OF.t2 = "T3", Pocock.t2 = "T4", OF.t3 = "T5",
              Pocock.t3 = "T6")

kt_power <- function(sf, tt, followup) {
  ms_power(scn_kt, n = n_kt, looks = ms_looks("os", fraction = tt),
           test = ms_test(alpha = 0.05, sf_eff = sf, mean = "schoenfeld"), followup = followup)
}

# probabilities of crossing the efficacy bound at each look under H1 and H0
# with the correlations of the Brownian motion in the information time t
kt_probs <- function(pw, tt) {
  sig <- sqrt(outer(tt, tt, pmin) / outer(tt, tt, pmax))
  lo <- rep(-Inf, length(tt))
  list(h1 = gs_prob(lo, pw$looks$z_upper, pw$looks$z_mean, sig)$prob_upper,
       h0 = gs_prob(lo, pw$looks$z_upper, 0, sig)$prob_upper)
}

kt_power_bm <- function(sf, tt, followup) {
  sum(kt_probs(kt_power(sf, tt, followup), tt)$h1)
}

kt_design <- function(sf, tt, followup) {
  pw <- kt_power(sf, tt, followup)
  pr <- kt_probs(pw, tt)
  k <- length(tt)
  s_end <- pw$study_duration
  f1 <- ms_events(scn_kt, n = n_kt, calendar_time = s_end)$events
  f0 <- ms_events(scn0_kt, n = n_kt, calendar_time = s_end)$events
  s0 <- c(ms_event_time(scn0_kt, n = n_kt, events = tt[-k] * f0), s_end)
  pat <- ms_events(scn0_kt, n = n_kt, calendar_time = s0)$enrolled
  st <- function(p) c(p[-k], 1 - sum(p[-k]))
  st1 <- st(pr$h1)
  st0 <- st(pr$h0)
  list(pw = pw, power = sum(pr$h1), tt = tt, f0 = f0, f1 = f1, s0 = s0, pat = pat, st0 = st0, st1 = st1,
       e_dur = c(H0 = sum(st0 * s0), H1 = sum(st1 * s0)),
       e_fail = c(H0 = sum(st0 * tt * f0), H1 = sum(st1 * tt * f1)),
       e_pat = c(H0 = sum(st0 * pat), H1 = sum(st1 * pat)))
}

pkg_expected <- function(d, quantity, hyp) {
  ex <- d$pw$expected
  ex[[quantity]][ex$hypothesis == if (hyp == "H1") "design" else "null"]
}

for (sf in c("OF", "Pocock")) {
  for (tk in names(times_kt)) {
    fu <- followup_kt[[sf]][[tk]]
    d <- kt_design(sf, times_kt[[tk]], fu)
    tab <- table_kt[[paste(sf, tk, sep = ".")]]
    # Tables 1-6
    for (k in seq_along(d$tt)) {
      b <- sprintf("%s-%s-%.1f", src, tab, d$tt[k])
      rows[[paste0(b, "-time")]] <- repro_row(paste0(b, "-time"), d$s0[k])
      rows[[paste0(b, "-bound")]] <- repro_row(paste0(b, "-bound"), d$pw$looks$z_upper[k],
                                               evidence = d$pw$looks$z_upper[k],
                                               evidence_criterion = "abs:0.01")
      # the same increment of the Pocock-type use function (from .6 to .8)
      # is printed in Table 2
      rows[[paste0(b, "-PH0")]] <- if (paste0(b, "-PH0") == "KimTsiatis1990-T4-0.8-PH0") {
        repro_row(paste0(b, "-PH0"), d$st0[k], evidence = d$st0[k],
                  evidence_ref_id = "KimTsiatis1990-T2-0.8-PH0")
      } else {
        repro_row(paste0(b, "-PH0"), d$st0[k])
      }
      rows[[paste0(b, "-PH1")]] <- repro_row(paste0(b, "-PH1"), d$st1[k], evidence = d$st1[k],
                                             evidence_criterion = "abs:0.002")
      rows[[paste0(b, "-F0")]] <- repro_row(paste0(b, "-F0"), d$tt[k] * d$f0)
      rows[[paste0(b, "-F1")]] <- repro_row(paste0(b, "-F1"), d$tt[k] * d$f1)
      rows[[paste0(b, "-patients")]] <- repro_row(paste0(b, "-patients"), d$pat[k])
    }
    # Section 6: follow-up (power at the printed follow-up and at 1.1 times
    # it) and expected study duration under H1
    h <- sprintf("%s-S6-%s-%s", src, sf, tk)
    pw_pert <- kt_power_bm(sf, d$tt, 1.1 * fu)
    rows[[paste0(h, "-followup")]] <- repro_row(
      paste0(h, "-followup"), d$power, computed_perturbed = pw_pert,
      note = sprintf(paste("power at the printed follow-up; ms_power()$power %.4f; ms_size() for the",
                           "power 0.90: follow-up %.4f"), d$pw$power,
                     ms_size(scn_kt, looks = ms_looks("os", fraction = d$tt),
                             test = ms_test(alpha = 0.05, sf_eff = sf, mean = "schoenfeld"),
                             power = 0.9, solve_for = "followup", n = n_kt)$solution))
    rows[[paste0(h, "-Edur-H1")]] <- repro_row(paste0(h, "-Edur-H1"), d$e_dur[["H1"]])
    rows[[paste0(h, "-Edur-H1-pkg")]] <- repro_row(
      paste0(h, "-Edur-H1-pkg"), pkg_expected(d, "study_duration", "H1"),
      note = "ms_power()$expected: analysis times of the design scenario")
    if (tk == "t1") {
      # Section 6: power of the design and of the same design at other
      # information times
      rows[[paste0(h, "-power")]] <- repro_row(paste0(h, "-power"), d$power,
                                               evidence = d$power,
                                               evidence_criterion = "abs:0.002",
                                               note = sprintf("ms_power()$power %.4f", d$pw$power))
      for (lk in list(c(0.15, 0.4, 0.7, 1), c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1))) {
        id <- sprintf("%s-power-%dlooks", h, length(lk))
        p <- kt_power_bm(sf, lk, fu)
        rows[[id]] <- repro_row(id, p, evidence = p, evidence_criterion = "abs:0.002")
      }
      # Section 5: expected study duration under H0 (evidence: with the
      # smallest follow-up that rounds to the printed one)
      h5 <- sprintf("%s-S5-%s-t1", src, sf)
      dm <- kt_design(sf, d$tt, fu - 0.005)
      rows[[paste0(h5, "-Edur-H0")]] <- repro_row(
        paste0(h5, "-Edur-H0"), d$e_dur[["H0"]], evidence = dm$e_dur[["H0"]],
        note = sprintf("with the follow-up %.3f years: %.4f", fu - 0.005, dm$e_dur[["H0"]]))
      rows[[paste0(h5, "-Edur-H0-pkg")]] <- repro_row(
        paste0(h5, "-Edur-H0-pkg"), pkg_expected(d, "study_duration", "H0"),
        note = "ms_power()$expected: analysis times of the null scenario")
      if (sf == "OF") {
        for (hy in c("H0", "H1")) {
          id <- sprintf("%s-Efail-%s", h5, hy)
          rows[[id]] <- repro_row(id, d$e_fail[[hy]])
          rows[[paste0(id, "-pkg")]] <- repro_row(paste0(id, "-pkg"),
                                                  pkg_expected(d, "events", hy),
                                                  note = "ms_power()$expected")
          id <- sprintf("%s-Epat-%s", h5, hy)
          # evidence for H0: the expected number cannot exceed the 250
          # patients accrued (Table 1, last row)
          rows[[id]] <- if (hy == "H0") {
            repro_row(id, d$e_pat[[hy]], evidence = d$e_pat[[hy]],
                      evidence_ref_id = "KimTsiatis1990-T1-1.0-patients")
          } else {
            repro_row(id, d$e_pat[[hy]])
          }
          rows[[paste0(id, "-pkg")]] <- repro_row(paste0(id, "-pkg"),
                                                  pkg_expected(d, "enrolled", hy),
                                                  note = "ms_power()$expected")
        }
      }
    }
  }
}

# Section 5: fixed-duration design (one analysis), power at the printed
# follow-up 1.02 years and at 1.1 times it
fixed_power <- function(fu) {
  ms_power(scn_kt, n = n_kt, looks = ms_looks("os"),
           test = ms_test(alpha = 0.05, mean = "schoenfeld"), followup = fu)$power
}
rows[["KimTsiatis1990-S5-fixed-followup"]] <- repro_row(
  "KimTsiatis1990-S5-fixed-followup", fixed_power(1.02),
  computed_perturbed = fixed_power(1.1 * 1.02))

result <- repro_finish(src, rows)
