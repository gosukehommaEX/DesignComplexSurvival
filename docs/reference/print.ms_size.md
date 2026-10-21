# Print an ms_size object

Print an ms_size object

## Usage

``` r
# S3 method for class 'ms_size'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An `ms_size` object.

- digits:

  Number of significant digits.

- ...:

  Passed to
  [`print.ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/print.ms_power.md).

## Value

`x`, invisibly.

## Examples

``` r
ctl <- ms_model(death = log(2) / 12)
trt <- ms_model(death = 0.7 * log(2) / 12)
scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
res <- ms_size(scn, looks = ms_looks("os"), test = ms_test(),
               power = 0.9, solve_for = "events", n = 500)
print(res)
#> Solved for the final number of events: 331.159 (target power 0.9)
#> 
#> Group sequential log-rank design
#>   Tested endpoint: OS; looks driven by OS events
#>   Sample size: 500; accrual duration: 25; final analysis at: 36.36
#>   Power: 0.9007
#>   Type I error: 0.025 (non-binding), 0.025 (binding)
#>  look  time events info_frac z_upper p_upper hr_upper z_lower hr_lower
#>     1 36.36    332         1    1.96   0.025   0.8064    1.96   0.8064
#>  prob_upper prob_lower cum_power
#>      0.9007    0.09927    0.9007
```
