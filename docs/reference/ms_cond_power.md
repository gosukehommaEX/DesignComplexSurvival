# Conditional power at an interim look

Computes, for a design with one endpoint, the probability that the trial
crosses an efficacy bound at a later look, given the standardized
log-rank statistic observed at an interim look, under a chosen
distribution of the statistics.

## Usage

``` r
ms_cond_power(
  x,
  look,
  z,
  under = c("design", "null", "trend", "hr", "scenario"),
  hr = NULL,
  scenario = NULL
)
```

## Arguments

- x:

  An `ms_power` object with one endpoint (from
  [`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)
  or
  [`ms_size`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md)).

- look:

  Interim look at which `z` is observed (from 1 to the number of looks
  minus 1).

- z:

  Observed standardized log-rank statistic at `look` (positive when the
  treatment is better), or the observed statistics at the looks 1 to
  `look`. It must be below the efficacy bounds and, with binding
  futility bounds, not below the futility bounds.

- under:

  Distribution of the statistics: `"design"`, `"null"`, `"trend"`,
  `"hr"` or `"scenario"` (see Details).

- hr:

  Hazard ratio (treatment versus control) for `under = "hr"`.

- scenario:

  An `ms_scenario` object with the allocation and accrual of the design
  scenario, for `under = "scenario"`.

## Value

An object of class `ms_cond_power`, a list with

- cond_power:

  The conditional power.

- cond_power_nonbinding:

  The conditional power without the futility bounds of the later interim
  looks.

- looks:

  A data frame with one row per later look: `look`, `time` and `events`
  (calendar time and expected events of the tested endpoint under the
  distribution of `under`), the bounds `z_upper` and `z_lower`, the
  conditional mean and standard deviation of the statistic (`cond_mean`,
  `cond_sd`), the conditional probabilities of stopping for efficacy and
  for futility at the look (`prob_upper`, `prob_lower`) and the
  cumulative conditional power `cum_upper`.

- look, z, under:

  The inputs.

- hr:

  The hazard ratio of `"trend"` and `"hr"` (`NA` otherwise).

## Details

At the interim look j with the observed statistic z, the statistics of
the later looks given \\Z_j = z\\ have the conditional normal
distribution of the joint normal distribution of the statistics chosen
by `under`. The conditional power is the probability that, at some later
look k, the statistic crosses the efficacy bound after the statistics of
the looks between j and k have stayed between the futility and efficacy
bounds, as in `gsCP()` of gsDesign; `cond_power_nonbinding` leaves out
the futility bounds of the later interim looks (a trial that continues
despite non-binding futility bounds). The bounds are those of the design
(`x$looks`). When the events observed at the interim look differ from
those planned, the bounds of the design no longer apply: compute the
design again with
[`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)
with looks that give the observed events and the remaining planned
events, and use the new design.

With `z` of length `look`, the statistics of the looks 1 to `look` are
all conditioned on. For statistics with independent increments
(`variance = "null"` of
[`ms_test`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md),
except OS with switching from a calendar time) the result is the same as
with the statistic of the look alone; with `"alternative"` and
`"studentized"` it can differ. In both cases the conditioning does not
use the fact that the trial did not stop at the earlier looks.

The distributions (`under`):

- `"design"`:

  The mean and covariance of the statistics under the design scenario
  (`x$mean`, `x$sigma`).

- `"null"`:

  Mean zero and the correlations under the null scenario
  (`x$sigma_null`), at the looks of the type I error (`x$looks_null`);
  the result is the conditional type I error.

- `"trend"`:

  The current trend: a constant hazard ratio `exp(-z / sqrt(p0 p1 d_j))`
  estimated from z, where p0 and p1 are the allocation probabilities and
  d the expected events of the tested endpoint at the looks of the
  design, so that the mean of the statistic at look k is
  `z sqrt(d_k / d_j)`, with the correlations `sqrt(d_j / d_k)` used for
  the bounds.

- `"hr"`:

  A constant hazard ratio `hr`: the mean at look k is
  `-log(hr) sqrt(p0 p1 d_k)`, with the correlations `sqrt(d_j / d_k)`.

- `"scenario"`:

  The mean and covariance under `scenario`, with the approximations of
  the test of the design; the event-driven looks take place when the
  expected events of their trigger endpoint reach the event targets of
  the design (`x$targets`), also when the final look of the design is
  fixed by `followup`, and the calendar-driven looks keep their times.

`"trend"` and `"hr"` assume proportional hazards (the Schoenfeld
approximation) and do not use the shape of the hazard ratio of the
design scenario, such as a delayed effect; `"scenario"` does.

## References

Ortega-Villa AM, Grieco MC, Rubenstein K, Wang J, Proschan MA (2025).
Futility monitoring in clinical trials. Statistics in Medicine, 44,
e70157.

## See also

[`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md),
[`gs_prob`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob.md)

## Examples

``` r
ctl <- ms_model(death = log(2) / 12)
trt <- ms_model(death = 0.7 * log(2) / 12)
scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
des <- ms_power(scn, n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
                test = ms_test(), events = 300)
ms_cond_power(des, look = 1, z = 1.2)
#> Conditional power of a group sequential log-rank design
#>   Observed at look 1: z = 1.2
#>   Distribution: design scenario
#>   Conditional power: 0.7269
#>   Without the futility bounds of the later interim looks: 0.7269
#>  look  time events z_upper z_lower cond_mean cond_sd prob_upper prob_lower
#>     2 39.67    300   1.969   1.969     2.393  0.7032     0.7269     0.2731
#>  cum_upper
#>     0.7269
ms_cond_power(des, look = 1, z = 1.2, under = "trend")
#> Conditional power of a group sequential log-rank design
#>   Observed at look 1: z = 1.2
#>   Distribution: current trend (hazard ratio 0.822)
#>   Conditional power: 0.3505
#>   Without the futility bounds of the later interim looks: 0.3505
#>  look  time events z_upper z_lower cond_mean cond_sd prob_upper prob_lower
#>     2 39.67    300   1.969   1.969     1.697  0.7071     0.3505     0.6495
#>  cum_upper
#>     0.3505
ms_cond_power(des, look = 1, z = 1.2, under = "hr", hr = 0.8)
#> Conditional power of a group sequential log-rank design
#>   Observed at look 1: z = 1.2
#>   Distribution: hazard ratio 0.8
#>   Conditional power: 0.4139
#>   Without the futility bounds of the later interim looks: 0.4139
#>  look  time events z_upper z_lower cond_mean cond_sd prob_upper prob_lower
#>     2 39.67    300   1.969   1.969     1.815  0.7071     0.4139     0.5861
#>  cum_upper
#>     0.4139
```
