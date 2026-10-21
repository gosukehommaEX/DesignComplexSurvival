# Print an ms_scenario object

Print an ms_scenario object

## Usage

``` r
# S3 method for class 'ms_scenario'
print(x, ...)
```

## Arguments

- x:

  An `ms_scenario` object.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
ctl <- ms_model(death = log(2) / 12)
trt <- ms_model(death = 0.7 * log(2) / 12)
print(ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20)))
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
