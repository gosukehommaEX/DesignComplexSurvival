# Joint distribution of the log-rank statistic at several calendar times

Returns the mean vector and covariance matrix of the standardized
log-rank statistics Z(t_1), ..., Z(t_K) of one endpoint, or of PFS and
OS together.

## Usage

``` r
lr_cov(
  scenario,
  n,
  calendar_time,
  endpoint = "os",
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

  Strictly increasing calendar times (positive).

- endpoint:

  `"os"`, `"pfs"` or both (for example `c("pfs", "os")`).

- mean:

  `"direct"` or `"schoenfeld"`.

- variance:

  `"null"`, `"alternative"` or `"studentized"`.

- trigger:

  Endpoint whose observed number of events determines each calendar time
  (`"os"` or `"pfs"`; one value or one per time), or `NULL` or `NA` for
  fixed calendar times. Used only with `variance = "studentized"`, and
  with switching from a calendar time.

## Value

A list with `time`, `endpoint`, `events`, `info`, `mean` (mean vector of
Z) and `sigma` (covariance matrix of Z), with one element (row) per
endpoint and time.

## Details

With `variance = "null"` the statistics have unit variances and
correlations `sqrt(V(t_j) / V(t_k))` for t_j \<= t_k, from the
independent increments of the score (Tsiatis 1981). With
`variance = "alternative"` the covariance of the score under the
alternative is that of Luo et al. (2019), computed here from the
per-patient influence functions of the score, and is divided by
`sqrt(V(t_j) V(t_k))`.

With `variance = "studentized"` the covariance is the asymptotic
covariance of Z = -U / sqrt(V-hat) itself, where V-hat is the usual
variance estimator. The per-patient influence function of Z adds to that
of the score a term for the variability of V-hat and, for times
determined by the observed number of events of an endpoint (`trigger`),
a term for the variability of the analysis time. Both terms vanish under
the null hypothesis.

For PFS and OS together (`endpoint = c("pfs", "os")`) the statistics are
ordered by endpoint and then by time. With `"studentized"` all
covariances come from the influence functions. With `"null"` the
statistics of one endpoint have the covariances described above and the
correlations between the endpoints are those of the studentized
statistics. With `"alternative"` the covariances between the endpoints
are the covariances of the scores divided by `sqrt(V_P V_O)`. Under the
null hypothesis the three choices agree, and the covariance between the
endpoints is that of Danzer et al. (2025, Theorem 1).

With switching from a calendar time (`switch_start` of
[`ms_scenario`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md)),
the hazard of OS depends on the entry time and the score has no
independent increments, even under the null hypothesis. The means,
information and covariances of OS then integrate over the entry time,
with the proportion in the risk set and the pooled hazard of each
calendar time. With `"null"` the statistics have unit variances and the
correlations of the studentized statistics (which, under the null
hypothesis, are the exact correlations), with `"alternative"` the
covariances of the scores divided by `sqrt(V_j V_k)`, and with
`"studentized"` the covariances of the studentized statistics; `trigger`
is used with all three. Only OS alone with `mean = "direct"` is covered;
PFS is not affected.

## References

Tsiatis AA (1981). The asymptotic joint distribution of the efficient
scores test for the proportional hazards model calculated over time.
Biometrika 68:311-315.

Luo X, Mao X, Chen X, Qiu J, Bai S, Quan H (2019). Design and monitoring
of survival trials in complex scenarios. Statistics in Medicine
38:192-209.

Danzer MF, Rufibach K, Beyersmann J, Schmidt R (2025). Exhausting the
type I error level in event-driven group-sequential designs with a
closed testing procedure for progression-free and overall survival.
arXiv:2512.08658.

## See also

[`lr_moments`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_moments.md),
[`gs_prob`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob.md)

## Examples

``` r
ctl <- ms_model(death = log(2) / 12)
trt <- ms_model(death = 0.7 * log(2) / 12)
scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
lr_cov(scn, n = 400, calendar_time = c(24, 36))
#> $time
#> [1] 24 36
#> 
#> $endpoint
#> [1] "os" "os"
#> 
#> $events
#> [1] 189.2029 281.1176
#> 
#> $info
#> [1] 47.04842 69.42161
#> 
#> $mean
#> [1] 2.446035 2.987523
#> 
#> $sigma
#>           [,1]      [,2]
#> [1,] 1.0000000 0.8232376
#> [2,] 0.8232376 1.0000000
#> 
fl <- ms_scenario(ms_fleischer(3, 9), ms_fleischer(4.3, 12.9),
                  accrual = ms_accrual(rate = 20))
lr_cov(fl, n = 400, calendar_time = 30, endpoint = c("pfs", "os"),
       variance = "studentized")
#> $time
#> [1] 30 30
#> 
#> $endpoint
#> [1] "pfs" "os" 
#> 
#> $events
#> [1] 383.8650 281.1558
#> 
#> $info
#> [1] 91.89643 69.35358
#> 
#> $mean
#> [1] 3.507546 3.015257
#> 
#> $sigma
#>          [,1]     [,2]
#> [1,] 0.983061 0.384386
#> [2,] 0.384386 0.999025
#> 
```
