# Print an ms_power object

Print an ms_power object

## Usage

``` r
# S3 method for class 'ms_power'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An `ms_power` object.

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
res <- ms_power(scn, n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
                test = ms_test(), events = 300)
print(res)
#> Group sequential log-rank design
#>   Tested endpoint: OS; looks driven by OS events
#>   Sample size: 400; accrual duration: 20; final analysis at: 39.67
#>   Power: 0.8689
#>   Type I error: 0.025 (non-binding), 0.025 (binding)
#>  look  time events info_frac z_upper  p_upper hr_upper z_lower hr_lower
#>     1 20.47    150       0.5   2.963 0.001525   0.6164    -Inf       NA
#>     2 39.67    300       1.0   1.969 0.024500   0.7967   1.969   0.7967
#>  prob_upper prob_lower cum_power
#>      0.2158     0.0000    0.2158
#>      0.6531     0.1311    0.8689
```
