# Per-look results of an ms_power object as a data frame

Per-look results of an ms_power object as a data frame

## Usage

``` r
# S3 method for class 'ms_power'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An `ms_power` object.

- row.names, optional:

  Not used.

- ...:

  Not used.

## Value

The data frame `x$looks`.

## Examples

``` r
ctl <- ms_model(death = log(2) / 12)
trt <- ms_model(death = 0.7 * log(2) / 12)
scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
res <- ms_power(scn, n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
                test = ms_test(), events = 300)
as.data.frame(res)
#>   look     time enrolled events events_control events_treatment info_frac
#> 1    1 20.46883      400    150       84.57556         65.42444       0.5
#> 2    2 39.66740      400    300      161.92098        138.07902       1.0
#>    z_upper     p_upper  hr_upper  z_lower  hr_lower   z_mean z_sd prob_upper
#> 1 2.962588 0.001525323 0.6164435     -Inf        NA 2.176151    1  0.2158058
#> 2 1.968596 0.024499771 0.7966707 1.968596 0.7966707 3.086900    1  0.6531038
#>   prob_lower cum_power alpha_spent   cum_alpha
#> 1  0.0000000 0.2158058 0.001525323 0.001525323
#> 2  0.1310905 0.8689095 0.025000000 0.025000000
```
