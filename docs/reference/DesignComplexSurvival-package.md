# DesignComplexSurvival: Simulation-Free Design of Complex Time-to-Event Trials

Computes the operating characteristics of two-group time-to-event trials
analysed with the log-rank test, without patient-level simulation. Each
group is described by a four-state model (progression-free, progressed,
switched, dead) with piecewise constant intensities, so that
progression-free survival, overall survival and treatment switching
after progression are handled within one model. Accrual is piecewise
uniform and dropout is piecewise exponential.

## Details

Expected numbers of events, the mean, variance and covariance of the
log-rank statistic over calendar time, and the probabilities of group
sequential designs with efficacy and futility bounds are obtained from
closed-form expressions, Gauss-Legendre quadrature and multivariate
normal probabilities.

A design may test one endpoint, progression-free survival (PFS) or
overall survival (OS), or both with a fixed-sequence or Bonferroni
procedure. With two endpoints, the analyses may be driven by the events
of either endpoint or by calendar time, and the joint distribution of
the PFS and OS statistics gives the power of each endpoint, the
probability of rejecting both or at least one hypothesis, and the
familywise error rate.

## Workflow

1.  Define each group with
    [`ms_model`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md)
    (or
    [`ms_fleischer`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_fleischer.md))
    and the accrual with
    [`ms_accrual`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_accrual.md),
    and combine them with
    [`ms_scenario`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md).

2.  Define the analysis times (and, with two endpoints, the endpoints
    tested at each analysis) with
    [`ms_looks`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_looks.md),
    and the test and the multiple testing procedure with
    [`ms_test`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md).

3.  Compute the power of a given design with
    [`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)
    or solve for the number of events, the sample size, the accrual
    duration or the follow-up with
    [`ms_size`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md).

The building blocks
[`ms_surv`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_surv.md),
[`ms_prob`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_prob.md),
[`ms_events`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_events.md),
[`ms_event_time`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_event_time.md),
[`lr_moments`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_moments.md),
[`lr_cov`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md),
[`gs_spend`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_spend.md),
[`gs_bounds`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_bounds.md),
[`gs_prob`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob.md)
and
[`gs_prob_mtp`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob_mtp.md)
are exported so that each step can be checked separately.

## See also

Useful links:

- <https://github.com/gosukehommaEX/DesignComplexSurvival>

- Report bugs at
  <https://github.com/gosukehommaEX/DesignComplexSurvival/issues>

## Author

**Maintainer**: Gosuke Homma <my.name.is.gosuke@gmail.com>
([ORCID](https://orcid.org/0000-0002-6854-5627))

Authors:

- Gosuke Homma <my.name.is.gosuke@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-6854-5627))
