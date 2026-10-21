# Two-group trial scenario

Combines the models of the control and treatment groups, the allocation
ratio and the accrual.

## Usage

``` r
ms_scenario(control, treatment, allocation = 1, accrual)
```

## Arguments

- control, treatment:

  `ms_model` objects of the control and treatment groups.

- allocation:

  Allocation ratio, treatment to control (for example 2 for 2:1).

- accrual:

  An `ms_accrual` object.

## Value

An object of class `ms_scenario`, a list with the two models,
`allocation`, the allocation probabilities `prob` (control, treatment)
and `accrual`.

## See also

[`ms_model`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md),
[`ms_accrual`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_accrual.md),
[`ms_null`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_null.md)

## Examples

``` r
ctl <- ms_model(death = log(2) / 12)
trt <- ms_model(death = 0.7 * log(2) / 12)
ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#> Two-group scenario, allocation 1:1 (treatment:control)
#> 
#> Control group. Four-state model for one group (clock after progression: "reset")
#> Before progression (time since randomization):
#>  time prog      death switch_prob dropout
#>     0    0 0.05776227           0       0
#> 
#> Treatment group. Four-state model for one group (clock after progression: "reset")
#> Before progression (time since randomization):
#>  time prog      death switch_prob dropout
#>     0    0 0.04043359           0       0
#> 
#> Piecewise uniform accrual by rates (patients per unit time):
#>  time rate
#>     0   20
```
