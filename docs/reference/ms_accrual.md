# Piecewise uniform accrual

Defines the accrual pattern either by accrual rates or by the
proportions of patients enrolled in intervals of a given accrual
duration.

## Usage

``` r
ms_accrual(time = 0, rate = NULL, prop = NULL, duration = NULL)
```

## Arguments

- time:

  Interval start times (calendar time from the start of the study). Must
  start at 0 and be strictly increasing.

- rate:

  Accrual rates in each interval (patients per unit time).

- prop:

  Proportions of patients enrolled in each interval (summing to 1).

- duration:

  Accrual duration (with `prop`).

## Value

An object of class `ms_accrual`.

## Details

With `rate`, patients are enrolled at `rate[j]` patients per unit time
from `time[j]` on, and the accrual duration follows from the sample
size. The last rate applies until the sample size is reached.

With `prop` and `duration`, a proportion `prop[j]` of the patients is
enrolled uniformly in the interval starting at `time[j]`, and the last
interval ends at `duration`. The accrual rates then scale with the
sample size. When
[`ms_size`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md)
solves for the accrual duration, all interval limits are scaled by the
same factor.

## See also

[`ms_scenario`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md)

## Examples

``` r
# 10 patients per month in the first 3 months, then 30 per month
ms_accrual(time = c(0, 3), rate = c(10, 30))
#> Piecewise uniform accrual by rates (patients per unit time):
#>  time rate
#>     0   10
#>     3   30

# 30 percent of the patients in the first 6 months of an 18-month accrual
ms_accrual(time = c(0, 6), prop = c(0.3, 0.7), duration = 18)
#> Piecewise uniform accrual by proportions, duration 18:
#>  time prop
#>     0  0.3
#>     6  0.7
```
