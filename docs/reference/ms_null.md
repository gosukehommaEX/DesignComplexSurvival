# Null scenario

Replaces the model of the treatment group by the model of the control
group, keeping the allocation and the accrual.

## Usage

``` r
ms_null(scenario)
```

## Arguments

- scenario:

  An `ms_scenario` object.

## Value

An `ms_scenario` object under the null hypothesis.

## See also

[`ms_scenario`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md)

## Examples

``` r
ctl <- ms_model(death = log(2) / 12)
trt <- ms_model(death = 0.7 * log(2) / 12)
scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
lr_moments(ms_null(scn), n = 400, calendar_time = 30)
#>   time events_control events_treatment   events score_mean     info z_mean
#> 1   30        133.442          133.442 266.8839          0 66.72098      0
#>   z_var
#> 1     1
```
