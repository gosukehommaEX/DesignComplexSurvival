# Expected numbers of patients, events and dropouts over calendar time

Expected numbers of patients, events and dropouts over calendar time

## Usage

``` r
ms_events(scenario, n, calendar_time, endpoint = c("os", "pfs"))
```

## Arguments

- scenario:

  An `ms_scenario` object.

- n:

  Total sample size.

- calendar_time:

  Calendar times from the start of the study.

- endpoint:

  `"os"` or `"pfs"`.

## Value

A data frame with `time`, the numbers enrolled (`enrolled`,
`enrolled_control`, `enrolled_treatment`), the expected events
(`events_control`, `events_treatment`, `events`) and the expected
dropouts (`dropouts_control`, `dropouts_treatment`).

## Details

For group g with allocation probability p_g, the expected number of
events by calendar time t is `n p_g` times the integral from 0 to t of
A(t - s) f_g(s) G_g(s), where A is the enrolled fraction, f_g the
density of the endpoint and G_g the probability of no dropout. The
expected number of dropouts is obtained in the same way with the dropout
density.

## See also

[`ms_event_time`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_event_time.md),
[`lr_moments`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_moments.md)

## Examples

``` r
ctl <- ms_model(death = log(2) / 12, dropout = 0.005)
trt <- ms_model(death = 0.7 * log(2) / 12, dropout = 0.005)
scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
ms_events(scn, n = 400, calendar_time = c(12, 24, 36))
#>   time enrolled enrolled_control enrolled_treatment events_control
#> 1   12      240              120                120       32.85127
#> 2   24      400              200                200      102.49878
#> 3   36      400              200                200      145.65790
#>   events_treatment    events dropouts_control dropouts_treatment
#> 1         24.47062  57.32188         2.843661           3.026026
#> 2         80.49241 182.99119         8.872469           9.953657
#> 3        121.46823 267.12613        12.608396          15.020710
```
