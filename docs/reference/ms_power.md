# Power of a group sequential log-rank design

Computes the analysis times, expected events, bounds, stopping
probabilities, power and type I error of a design with a given sample
size, without simulation.

## Usage

``` r
ms_power(
  scenario,
  n,
  looks,
  test,
  events = NULL,
  followup = NULL,
  null_scenarios = NULL
)
```

## Arguments

- scenario:

  An `ms_scenario` object.

- n:

  Total sample size.

- looks:

  An `ms_looks` object.

- test:

  An `ms_test` object.

- events:

  Final number of events of the trigger endpoint (for event-driven
  looks), or final numbers of events named by endpoint (`"os"`,
  `"pfs"`).

- followup:

  Time from the end of accrual to the final analysis (for an
  event-driven final look, instead of `events`).

- null_scenarios:

  For two endpoints, an optional list of `ms_scenario` objects with the
  allocation and accrual of `scenario`, each named by the endpoint
  (`"os"` or `"pfs"`) whose null hypothesis is true in the scenario, for
  example `list(os = s_os)`.

## Value

An object of class `ms_power`, a list with

- looks:

  A data frame with one row per look: `time`, `enrolled`, the expected
  events of the tested endpoint (`events`, `events_control`,
  `events_treatment`), `info_frac`, the efficacy bound (`z_upper`,
  nominal one-sided p-value `p_upper`, hazard ratio `hr_upper`), the
  futility bound (`z_lower`, `hr_lower`), the mean and standard
  deviation of the statistic (`z_mean`, `z_sd`), the probabilities of
  stopping for efficacy and futility at the look (`prob_upper`,
  `prob_lower`), the cumulative power `cum_power`, the cumulative alpha
  spent `alpha_spent` and the cumulative non-binding type I error
  `cum_alpha`. When the looks are driven by another endpoint,
  `events_trigger` gives the events of that endpoint.

- power:

  Overall power.

- type1:

  Type I error, non-binding and binding.

- n, accrual_duration, study_duration, events:

  Sample size, accrual duration, calendar time of the final analysis and
  expected events of the tested endpoint at the final analysis.

- settings, scenario, looks_def, test:

  The inputs.

With two endpoints, `looks` has one row per endpoint and look at which
the endpoint is tested, with the columns above and `endpoint`, `look`,
`trigger` and `events_trigger` (the events of the trigger endpoint at
the look, `NA` for calendar-driven looks); `prob_upper` is the
probability of rejecting the hypothesis of the endpoint at the look
under the multiple testing procedure (with recycling of the level,
`z_upper_recycled` and `p_upper_recycled` give the efficacy bound and
its nominal p-value at the level after the other hypothesis has been
rejected, after `hr_upper`), `prob_lower` the probability of stopping
for futility at the look (first endpoint), `cum_power` the cumulative
sum of `prob_upper`, `alpha_spent` the cumulative alpha of the bounds of
the endpoint and `cum_alpha` the cumulative non-binding probability of
rejecting the hypothesis of the endpoint under the global null
hypothesis. `power` has the powers of the two endpoints (named by
endpoint), the probability of rejecting both hypotheses (`both`) and of
rejecting at least one (`any`); `type1` is the familywise error rate
under the global null hypothesis, non-binding and binding; `fwer` is a
data frame with the familywise error rate under the global null
hypothesis and under each scenario of `null_scenarios`; `events` has the
expected events of each endpoint at its last look; `final_events` the
final numbers of events of the trigger endpoints; and `mean` and `sigma`
are the mean vector and covariance matrix of the statistics in the order
of the rows of `looks`.

## Details

For event-driven looks, the final analysis takes place when the expected
number of events of the trigger endpoint reaches `events`, or at
`followup` after the end of accrual. The interim analyses take place
when the expected number of events reaches the fractions given in
`looks`. For calendar-driven looks, `events` and `followup` are not
used.

At each look, the standardized log-rank statistics have the mean and
covariance of
[`lr_cov`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md)
with the approximations chosen in `test`. The bounds are computed by
[`gs_bounds`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_bounds.md)
and the stopping probabilities by
[`gs_prob`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob.md).
The power is the probability of crossing an efficacy bound when the
futility bounds are applied. The type I error is computed under
[`ms_null`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_null.md)
applied to the scenario, with event-driven looks at the same numbers of
events of the trigger endpoint; the non-binding value ignores the
futility bounds and the binding value applies them.

The power and the type I error rest on the normal approximation of the
log-rank statistics. With unequal allocation and few events the
statistic is skewed, so that the probability of rejecting in favor of
the treatment exceeds the computed type I error when the treatment group
is the larger group, and falls below it when the control group is: by
simulation (`inst/validation` of the source package), about 0.028
instead of 0.025 with 2:1 allocation and looks at 75 and 150 deaths.

The hazard ratios at the bounds are `exp(-z / sqrt(p0 p1 d))`, where d
is the expected number of events of the tested endpoint.

When the looks have different triggers (see
[`ms_looks`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_looks.md)),
each trigger endpoint has its own final number of events, and a look
driven by it takes place when its expected number of events reaches the
fraction of the look times that number. The final numbers are given by
`events`, named by endpoint when the looks are driven by the events of
two endpoints, for example `c(pfs = 300, os = 250)`. With `followup`, or
without `events` when the final look is calendar-driven, they are the
expected numbers of events at the final analysis; the looks driven by an
endpoint other than that of the final look then need fractions below 1.

With two endpoints (`test` from
[`ms_test`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md)
with two endpoints), the statistics of the two endpoints at the looks at
which they are tested have the joint distribution of
[`lr_cov`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md),
with the observed events of the trigger endpoints determining the
event-driven looks for `variance = "studentized"`. Each endpoint has its
own bounds, and the rejection probabilities under the multiple testing
procedure are computed by
[`gs_prob_mtp`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob_mtp.md);
with recycling of the level (`transition` of
[`ms_test`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md)),
each endpoint also has bounds at its level after recycling. The
familywise error rate is computed under
[`ms_null`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_null.md)
applied to the scenario (the global null hypothesis), where the
correlations between the endpoints are those of Danzer et al. (2025,
Theorem 1), and under each scenario of `null_scenarios`, where the
hypothesis of one endpoint is true and that of the other is not (a
partial null hypothesis). Such scenarios cannot be derived from the
multistate model of the design scenario and are given by the user, for
example a model with the same OS in both groups and a difference in PFS.
Under each null scenario the event-driven looks take place at the same
numbers of events of their trigger endpoints.

With switching from a calendar time (`switch_start` of
[`ms_scenario`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md)),
only a test of one endpoint with `mean = "direct"` for OS is covered,
and an OS test cannot use looks driven by PFS events. For OS the
covariance of the statistics comes from the per-patient influence
functions (see
[`lr_cov`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md));
the bounds still come from the information or event fractions, and the
type I error under
[`ms_null`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_null.md)
uses the exact correlations under the null hypothesis, so it can differ
from `alpha` when the switching changes the death hazard.

## References

Danzer MF, Rufibach K, Beyersmann J, Schmidt R (2025). Exhausting the
type I error level in event-driven group-sequential designs with a
closed testing procedure for progression-free and overall survival.
arXiv:2512.08658.

## See also

[`ms_size`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md),
[`ms_looks`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_looks.md),
[`ms_test`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md),
[`gs_prob_mtp`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob_mtp.md)

## Examples

``` r
ctl <- ms_model(death = log(2) / 12)
trt <- ms_model(death = 0.7 * log(2) / 12)
scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
ms_power(scn, n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
         test = ms_test(), events = 300)
#> Group sequential log-rank design
#>   Tested endpoint: OS; looks driven by OS events
#>   Sample size: 400; accrual duration: 20; final analysis at: 39.67
#>   Power: 0.8689
#>   Type I error: 0.025 (non-binding), 0.025 (binding)
#>  look  time events info_frac z_upper  p_upper hr_upper z_lower hr_lower
#>     1 20.47    150       0.5   2.963 0.001525   0.6164    -Inf       NA
#>     2 39.67    300       1.0   1.969 0.024500   0.7967   1.969   0.7967
#>  prob_upper prob_lower cum_power
#>      0.2158     0.0000    0.2158
#>      0.6531     0.1311    0.8689

# PFS first, then OS: final PFS analysis at 300 PFS events with an OS
# interim analysis, final OS analysis at 220 OS events
fl <- ms_scenario(ms_fleischer(4, 12), ms_fleischer(6, 16),
                  accrual = ms_accrual(rate = 25))
lk <- ms_looks(c("pfs", "os"), fraction = c(1, 1), test = list(c("pfs", "os"), "os"))
ms_power(fl, n = 450, looks = lk,
         test = ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence"),
         events = c(pfs = 300, os = 220))
#> Group sequential log-rank design with two endpoints
#>   Endpoints: PFS, then OS (fixed sequence, overall); looks driven by PFS events and OS events
#>   Sample size: 450; accrual duration: 18; final analysis at: 23.05
#>   Power: PFS 0.9387, OS 0.5452; both 0.5452; at least one 0.9387
#>   FWER under global null: 0.025 (non-binding), 0.025 (binding)
#>  endpoint look  time events info_frac z_upper  p_upper hr_upper z_lower
#>       PFS    1 18.62    300    1.0000   1.960 0.025000   0.7975   1.960
#>        OS    1 18.62    163    0.7407   2.357 0.009207   0.6912    -Inf
#>        OS    2 23.05    220    1.0000   2.009 0.022243   0.7626   2.009
#>  hr_lower prob_upper prob_lower cum_power
#>    0.7975     0.9387          0    0.9387
#>        NA     0.2966          0    0.2966
#>    0.7626     0.2486          0    0.5452
```
