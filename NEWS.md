# DesignComplexSurvival 0.1.0

* First version (stage 1 of the development plan): designs that test one
  endpoint, progression-free survival (PFS) or overall survival (OS).
* Model: a four-state model per group (progression-free, progressed,
  switched, dead) with piecewise constant intensities, immediate switching
  at progression and switching after progression, and a choice of the clock
  of the post-progression intensities (`clock = "reset"` or `"forward"`)
  (`ms_model()`, `ms_fleischer()`); piecewise uniform accrual by rates or by
  proportions (`ms_accrual()`); piecewise exponential dropout by group.
* Endpoint distributions and state occupation probabilities (`ms_surv()`,
  `ms_prob()`), expected events and dropouts over calendar time
  (`ms_events()`, `ms_event_time()`).
* Mean, variance and covariance of the log-rank statistic over calendar time
  (`lr_moments()`, `lr_cov()`), with the mean by the direct formula or the
  Schoenfeld approximation and the covariance with unit variances or under
  the alternative (Luo et al. 2019).
* Group sequential designs: alpha spending for efficacy; no futility, fixed
  futility bounds (Z or hazard ratio scale) or beta spending; binding or
  non-binding futility (`gs_spend()`, `gs_bounds()`, `gs_prob()`).
* Power and type I error (`ms_power()`) and the number of events, sample
  size, accrual duration or follow-up for a target power (`ms_size()`), with
  event-driven or calendar-driven looks (`ms_looks()`, `ms_test()`).
* Tests compare the results with an independent Python implementation
  (`inst/validation/python/`) and, when installed, with gsDesign, gsDesign2,
  rpact and lrstat.
