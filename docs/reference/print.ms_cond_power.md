# Print an `ms_cond_power` object

Print an `ms_cond_power` object

## Usage

``` r
# S3 method for class 'ms_cond_power'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An `ms_cond_power` object.

- digits:

  Number of significant digits.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
ctl <- ms_model(death = log(2) / 12)
trt <- ms_model(death = 0.7 * log(2) / 12)
scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
des <- ms_power(scn, n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
                test = ms_test(), events = 300)
print(ms_cond_power(des, look = 1, z = 1.2, under = "trend"))
#> Conditional power of a group sequential log-rank design
#>   Observed at look 1: z = 1.2
#>   Distribution: current trend (hazard ratio 0.822)
#>   Conditional power: 0.3505
#>   Without the futility bounds of the later interim looks: 0.3505
#>  look  time events z_upper z_lower cond_mean cond_sd prob_upper prob_lower
#>     2 39.67    300   1.969   1.969     1.697  0.7071     0.3505     0.6495
#>  cum_upper
#>     0.3505
```
