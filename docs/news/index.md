# Changelog

## DesignComplexSurvival 0.3.3

- Documentation: the introductory vignette covers the workflow of a
  design with one endpoint, and five new vignettes describe the
  multistate model and treatment switching, group sequential designs and
  the covariance of the statistics, PFS and OS tested together,
  switching from a calendar time, and the validation of the package. The
  pkgdown site lists them under Articles.
- Validation: in `inst/validation/validate_fastsurvival.R` (V1-V4), the
  death that defines a look driven by OS deaths was lost in some trials,
  because the follow-up time cut - entry could be smaller than its time
  by rounding (on average, about 0.1 to 0.3 deaths fewer than the target
  at these looks). The comparisons now use the cut times (1 + 1e-12),
  and the script checks that every trial has the target number of deaths
  and that the log-rank Z agrees with
  [`survival::survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html).
  The package code is unchanged.

## DesignComplexSurvival 0.3.2

- Documentation:
  [`ms_power()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)
  notes that the normal approximation does not describe the skewness of
  the log-rank statistic with unequal allocation and few events (by
  simulation, a type I error of about 0.028 instead of 0.025 with 2:1
  allocation and looks at 75 and 150 deaths).
- Validation: larger simulations of V2 and V7 under the null hypothesis
  (`inst/validation/python/simulate_type1.py`): the excess type I error
  of V2 comes from this skewness and that of V7 is a chance deviation
  (`inst/validation/README.md`).

## DesignComplexSurvival 0.3.1

- Faster designs whose final analysis is fixed by `followup`
  ([`ms_power()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)
  with `followup`, and the searches of
  [`ms_size()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md)):
  the time of the final analysis, `followup` after the end of accrual,
  is no longer also searched for from its expected number of events. The
  results are unchanged. The gain is largest with switching from a
  calendar time, where these times come from Newton’s method.
- Tests of
  [`ms_size()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md)
  with recycling of the level and with switching from a calendar time,
  with expected values from
  `inst/validation/python/expected_size_stage3.py`.

## DesignComplexSurvival 0.3.0

- Stage 3 of the development plan: further models of switching and
  censoring, switching that starts at a calendar time, and recycling of
  the level between PFS and OS. With the default values of the new
  arguments the results are unchanged: the outputs of `inst/reproduce/`
  and `inst/validation/` of version 0.2.2 are reproduced without any
  numerical difference.
- [`ms_model()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md)
  gains `switch_clock`: with `"switch"`, the death hazard after a switch
  is a function of the time since the switch instead of the
  post-progression clock (`"same"`, the default).
- [`ms_model()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md)
  gains `switch_af`, an accelerated failure time effect of the switch:
  the remaining time after the switch is the remaining time without
  switch divided by `switch_af` (the rank preserving structural failure
  time model of Robins and Tsiatis 1991).
- [`ms_model()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md)
  gains switching before progression: `switch_pre` (hazard of switching
  while progression-free), `prog_switch` and `death_switch_pre` (hazards
  of progression and death after such a switch). After a later
  progression these patients die with the hazard `death_switch` (with
  `switch_af`, `switch_af * death_prog`).
  [`ms_prob()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_prob.md)
  then adds the column `progression_free_switched`.
- [`ms_scenario()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md)
  gains `pfs_switch`: PFS is followed up to progression or death
  (`"event"`, the default) or censored at a switch before progression
  (`"censor"`).
- [`ms_model()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md)
  gains `dropout_pfs`, a censoring of PFS only (for example the end of
  tumour assessments). The PFS dropouts of
  [`ms_events()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_events.md)
  include it and, with `pfs_switch = "censor"`, the censoring at a
  switch.
- [`ms_scenario()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md)
  gains `switch_start`, a calendar time before which nobody switches
  (for example when switching is allowed only after an interim
  analysis). The distribution of OS then depends on the entry time, and
  the moments and covariances of the OS statistic integrate over the
  entry time. Under the null hypothesis the score no longer has
  independent increments, and the correlations between the looks come
  from the per-patient influence functions. The OS statistic alone with
  `mean = "direct"` is covered, at calendar times or at looks driven by
  OS events; PFS is not affected. Two endpoints together give an error.
- [`ms_test()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md)
  gains `transition` for `mtp = "bonferroni"`: the level of a rejected
  hypothesis is passed on to the other (graphical procedure of Maurer
  and Bretz 2013, Algorithm 1), and
  [`gs_prob_mtp()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob_mtp.md)
  gains `upper_recycled`. The bounds after recycling are checked not to
  exceed the bounds at the initial level at every look. With recycling,
  the table of
  [`ms_power()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)
  gains the columns `z_upper_recycled` and `p_upper_recycled`.
- Reproduction scripts added for Maurer and Bretz (2013) and
  Torres-Saavedra et al. (2026) (`inst/reproduce/`).
- Validation: Python simulations of the new options (switch clock and
  accelerated failure time effect; switching before progression and
  censoring of PFS only, with PFS and OS together; switching from a
  calendar time), the Bonferroni procedure with recycling in the
  simulations of PFS and OS together, and a comparison of the rejection
  probabilities with the decisions of
  [`lrstat::fseqbon()`](https://kaifenglu.github.io/lrstat/reference/fseqbon.html)
  and of the bounds with
  [`lrstat::getBound()`](https://kaifenglu.github.io/lrstat/reference/getBound.html)
  (`inst/validation/`).

## DesignComplexSurvival 0.2.2

- Bug fix: rectangle probabilities of nine or more coordinates with two
  finite limits for some coordinate (designs with futility bounds and
  nine or more looks of an endpoint) were wrong with the Miwa algorithm
  of mvtnorm, which leaves out the terms of nine or more coordinates of
  the inclusion-exclusion over the limits (errors of 7.5e-4 with nine
  and 6.9e-3 with ten coordinates under the correlation of independent
  increments, so also with the default `variance = "null"`). They now
  use the Genz-Bretz algorithm with a fixed seed. Rectangles with
  one-sided limits only (no futility bounds) and designs with at most
  eight looks per endpoint are unchanged.

## DesignComplexSurvival 0.2.1

- Bug fix: with `variance = "alternative"` or `"studentized"` and four
  or more looks of an endpoint, the stopping probabilities, the power,
  the futility bounds by beta spending and the familywise error rates
  under `null_scenarios` were computed with the Miwa algorithm, which is
  inaccurate for covariance matrices that are not those of statistics
  with independent increments (errors of about 1e-3 with four looks and
  above 0.1 with six or more looks, found by comparison with exact
  values). Rectangle probabilities of four or more coordinates with such
  a covariance matrix now use the Genz-Bretz algorithm with a fixed
  seed, in
  [`gs_prob()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob.md),
  [`gs_bounds()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_bounds.md)
  and
  [`gs_prob_mtp()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob_mtp.md)
  as well. The default `variance = "null"`, the efficacy bounds, the
  type I error and designs with at most three looks per endpoint are
  unchanged.
- The option for the absolute error tolerance of the Genz-Bretz
  algorithm is renamed `DesignComplexSurvival.genz_abseps` (formerly
  `DesignComplexSurvival.mtp_abseps`); the default remains 1e-5.
- With binding futility bounds, a futility bound that reaches the
  efficacy bound at an interim look stops every trial there, so that the
  efficacy bounds of the later looks cannot be computed; the error
  message now names that look.

## DesignComplexSurvival 0.2.0

- Stage 2 of the development plan: designs that test progression-free
  survival (PFS) and overall survival (OS) together.
- [`ms_test()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md)
  gains `mtp` (`"none"`, `"bonferroni"` or `"fixed_sequence"`),
  `weights` and `hierarchical` (`"overall"` or `"stagewise"`,
  procedures (b) and (a) of Glimm et al. 2010), added after the existing
  arguments. `endpoints` may name both endpoints in the order of
  testing, and `sf_eff` and `param_eff` may differ between them. Each
  endpoint has its own efficacy bounds; futility bounds are placed on
  the first endpoint.
- [`ms_looks()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_looks.md)
  accepts one trigger per look (`"os"`, `"pfs"` or `"calendar"`) and,
  through `test`, the endpoints tested at each look.
- [`ms_power()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)
  returns, for two endpoints, the power of each endpoint, the
  probability of rejecting both and at least one hypothesis, and the
  familywise error rate under the global null hypothesis and under
  user-given partial null scenarios (`null_scenarios`). `events` takes
  the final numbers of events named by endpoint, for example
  `c(pfs = 300, os = 220)`.
- [`ms_size()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md)
  gains `target` (`"first"`, `"each"`, `"both"` or `"any"`) and
  `null_scenarios`; with looks driven by the events of both endpoints,
  `solve_for = "events"` rounds up the final numbers of events of both.
- [`lr_cov()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md)
  returns the joint distribution of the PFS and OS statistics
  (`endpoint = c("pfs", "os")`). Under the null hypothesis their
  covariance is that of Danzer et al. (2025, Theorem 1).
- New option `variance = "studentized"` in
  [`lr_moments()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_moments.md),
  [`lr_cov()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md)
  and
  [`ms_test()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md):
  the asymptotic covariance of the standardized statistic itself,
  including the variability of the variance estimator and, through
  `trigger`, of the analysis times driven by observed numbers of events.
  The default remains `variance = "null"`, and the results of the
  existing options are unchanged.
- New function
  [`gs_prob_mtp()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob_mtp.md):
  rejection probabilities of the two hypotheses under the Bonferroni and
  fixed-sequence procedures. Rectangles of four or more coordinates over
  both endpoints use the Genz-Bretz algorithm with a fixed seed
  (absolute error set by
  `options(DesignComplexSurvival.mtp_abseps = ...)`, default 1e-5); the
  other probabilities use the Miwa algorithm as before. `mvtnorm` (\>=
  1.2-0) is now required.
- [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  handle two endpoints.
- Reproduction scripts added for Glimm et al. (2010), Tamhane et
  al. (2010), Sugimoto et al. (2017) and Erdmann et al. (2025)
  (`inst/reproduce/`).
- Validation: a Python simulation of PFS and OS analysed together in
  four scenarios, with designs of two endpoints, and the covariance
  estimator of Danzer et al. (2025) applied to large simulated trials
  under the null hypothesis (`inst/validation/`). The Python simulation
  of single endpoints now counts the event at an event-driven analysis
  time in every trial.

## DesignComplexSurvival 0.1.0

- First version (stage 1 of the development plan): designs that test one
  endpoint, progression-free survival (PFS) or overall survival (OS).
- Model: a four-state model per group (progression-free, progressed,
  switched, dead) with piecewise constant intensities, immediate
  switching at progression and switching after progression, and a choice
  of the clock of the post-progression intensities (`clock = "reset"` or
  `"forward"`)
  ([`ms_model()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md),
  [`ms_fleischer()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_fleischer.md));
  piecewise uniform accrual by rates or by proportions
  ([`ms_accrual()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_accrual.md));
  piecewise exponential dropout by group.
- Endpoint distributions and state occupation probabilities
  ([`ms_surv()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_surv.md),
  [`ms_prob()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_prob.md)),
  expected events and dropouts over calendar time
  ([`ms_events()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_events.md),
  [`ms_event_time()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_event_time.md)).
- Mean, variance and covariance of the log-rank statistic over calendar
  time
  ([`lr_moments()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_moments.md),
  [`lr_cov()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md)),
  with the mean by the direct formula or the Schoenfeld approximation
  and the covariance with unit variances or under the alternative (Luo
  et al. 2019).
- Group sequential designs: alpha spending for efficacy; no futility,
  fixed futility bounds (Z or hazard ratio scale) or beta spending;
  binding or non-binding futility
  ([`gs_spend()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_spend.md),
  [`gs_bounds()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_bounds.md),
  [`gs_prob()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob.md)).
- Power and type I error
  ([`ms_power()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md))
  and the number of events, sample size, accrual duration or follow-up
  for a target power
  ([`ms_size()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md)),
  with event-driven or calendar-driven looks
  ([`ms_looks()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_looks.md),
  [`ms_test()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md)).
- Tests compare the results with an independent Python implementation
  (`inst/validation/python/`) and, when installed, with gsDesign,
  gsDesign2, rpact and lrstat.
- Scripts that reproduce published sample sizes, numbers of events,
  analysis times, bounds and powers (`inst/reproduce/`) and compare the
  package with patient-level simulations (FastSurvival and Python) and
  with PWEALL (`inst/validation/`). Each folder writes its judgements to
  `results/summary.csv` and `results/summary.md`.
