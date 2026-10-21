# Package index

## Model and scenario

- [`ms_model()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md)
  : Multistate model for one treatment group
- [`ms_fleischer()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_fleischer.md)
  : Exponential model for PFS and OS of Fleischer et al. (2009)
- [`ms_accrual()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_accrual.md)
  : Piecewise uniform accrual
- [`ms_scenario()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md)
  : Two-group trial scenario
- [`ms_null()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_null.md)
  : Null scenario

## Endpoint distributions

- [`ms_surv()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_surv.md)
  : Survival function, hazard and density of an endpoint
- [`ms_prob()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_prob.md)
  : State occupation probabilities

## Events and the log-rank statistic

- [`ms_events()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_events.md)
  : Expected numbers of patients, events and dropouts over calendar time
- [`ms_event_time()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_event_time.md)
  : Calendar time at which a number of events is expected
- [`lr_moments()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_moments.md)
  : Mean and variance of the log-rank statistic over calendar time
- [`lr_cov()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md)
  : Joint distribution of the log-rank statistic at several calendar
  times

## Group sequential bounds and probabilities

- [`gs_spend()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_spend.md)
  : Error spending functions
- [`gs_bounds()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_bounds.md)
  : Bounds of a group sequential design
- [`gs_prob()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob.md)
  : Stopping probabilities of a group sequential design
- [`gs_prob_mtp()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob_mtp.md)
  : Rejection probabilities of a group sequential design with two
  endpoints

## Design

- [`ms_looks()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_looks.md)
  : Timing of the analyses
- [`ms_test()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md)
  : Test and stopping rules
- [`ms_power()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)
  : Power of a group sequential log-rank design
- [`ms_size()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md)
  : Number of events, sample size, accrual duration or follow-up for a
  target power

## Methods

- [`print(`*`<ms_model>`*`)`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/print.ms_model.md)
  : Print an ms_model object
- [`print(`*`<ms_accrual>`*`)`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/print.ms_accrual.md)
  : Print an ms_accrual object
- [`print(`*`<ms_scenario>`*`)`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/print.ms_scenario.md)
  : Print an ms_scenario object
- [`print(`*`<ms_looks>`*`)`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/print.ms_looks.md)
  : Print an ms_looks object
- [`print(`*`<ms_test>`*`)`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/print.ms_test.md)
  : Print an ms_test object
- [`print(`*`<ms_power>`*`)`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/print.ms_power.md)
  : Print an ms_power object
- [`summary(`*`<ms_power>`*`)`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/summary.ms_power.md)
  : Summarize an ms_power object
- [`print(`*`<summary.ms_power>`*`)`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/print.summary.ms_power.md)
  : Print a summary.ms_power object
- [`as.data.frame(`*`<ms_power>`*`)`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/as.data.frame.ms_power.md)
  : Per-look results of an ms_power object as a data frame
- [`print(`*`<ms_size>`*`)`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/print.ms_size.md)
  : Print an ms_size object
