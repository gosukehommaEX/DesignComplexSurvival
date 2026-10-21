# Mean and variance of the log-rank statistic over calendar time

Computes, at each calendar time, the expected events, the mean of the
log-rank score U (observed minus expected events in the treatment
group), its information V (the limit of the usual variance estimator)
and the mean and variance of the standardized statistic Z = -U /
sqrt(V), which is positive when the treatment is better.

## Usage

``` r
lr_moments(
  scenario,
  n,
  calendar_time,
  endpoint = c("os", "pfs"),
  mean = c("direct", "schoenfeld"),
  variance = c("null", "alternative", "studentized"),
  trigger = NULL
)
```

## Arguments

- scenario:

  An `ms_scenario` object.

- n:

  Total sample size.

- calendar_time:

  Calendar times (positive).

- endpoint:

  `"os"` or `"pfs"`.

- mean:

  `"direct"` or `"schoenfeld"`.

- variance:

  `"null"`, `"alternative"` or `"studentized"`.

- trigger:

  Endpoint whose observed number of events determines each calendar time
  (`"os"` or `"pfs"`; one value or one per time), or `NULL` or `NA` for
  fixed calendar times. Used only with `variance = "studentized"`.

## Value

A data frame with `time`, `events_control`, `events_treatment`,
`events`, `score_mean`, `info`, `z_mean` and `z_var`.

## Details

Let y0 and y1 be the probabilities of being at risk at follow-up time s
in each group, multiplied by the allocation probabilities and by the
enrolled fraction at t - s. Then E(U) = n times the integral of y0 y1 /
(y0 + y1) (lambda1 - lambda0) and V = n times the integral of y0 y1 /
(y0 + y1)^2 (y0 lambda0 + y1 lambda1), from 0 to t.

The mean of Z is

- `mean = "direct"`: -E(U) / sqrt(V), the mean of Yung and Liu (2020,
  Section 2.3.3);

- `mean = "schoenfeld"`: -sqrt(p0 p1 D) times the event-weighted average
  of the log hazard ratio, where D is the expected number of events
  (Schoenfeld 1981).

The variance of Z is 1 for `variance = "null"`. For
`variance = "alternative"` it is the asymptotic variance of the score
under the alternative divided by V, from Luo et al. (2019), which
accounts for the variability of the numbers at risk. This variance
assumes simple randomization; Yung and Liu (2020) show that block
randomization gives a variance that is not larger. For
`variance = "studentized"` it is the asymptotic variance of Z itself,
which also accounts for the variability of V and, for times determined
by the observed number of events of an endpoint (`trigger`), of the
analysis time (see
[`lr_cov`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md)).

## References

Luo X, Mao X, Chen X, Qiu J, Bai S, Quan H (2019). Design and monitoring
of survival trials in complex scenarios. Statistics in Medicine
38:192-209.

Schoenfeld D (1981). The asymptotic properties of nonparametric tests
for comparing survival distributions. Biometrika 68:316-319.

Yung G, Liu Y (2020). Sample size and power for the weighted log-rank
test and Kaplan-Meier based tests with allowance for nonproportional
hazards. Biometrics 76:939-950.

## See also

[`lr_cov`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md),
[`ms_events`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_events.md)

## Examples

``` r
ctl <- ms_model(death = log(2) / 12)
trt <- ms_model(death = 0.7 * log(2) / 12)
scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
lr_moments(scn, n = 400, calendar_time = c(24, 36))
#>   time events_control events_treatment   events score_mean     info   z_mean
#> 1   24       105.8727         83.33015 189.2029  -16.77781 47.04842 2.446035
#> 2   36       152.9364        128.18129 281.1176  -24.89193 69.42161 2.987523
#>   z_var
#> 1     1
#> 2     1
```
