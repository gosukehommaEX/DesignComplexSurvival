# Calendar time at which a number of events is expected

Solves for the calendar time at which the expected total number of
events of an endpoint reaches each target.

## Usage

``` r
ms_event_time(scenario, n, events, endpoint = c("os", "pfs"))
```

## Arguments

- scenario:

  An `ms_scenario` object.

- n:

  Total sample size.

- events:

  Target numbers of events (positive).

- endpoint:

  `"os"` or `"pfs"`.

## Value

Numeric vector of calendar times. An error is raised when a target
cannot be reached.

## See also

[`ms_events`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_events.md)

## Examples

``` r
ctl <- ms_model(death = log(2) / 12)
trt <- ms_model(death = 0.7 * log(2) / 12)
scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
ms_event_time(scn, n = 400, events = c(150, 300))
#> [1] 20.46883 39.66740
```
