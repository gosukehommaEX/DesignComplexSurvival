# Changelog

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
