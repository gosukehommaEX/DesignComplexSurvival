# Summarize an ms_power object

Summarize an ms_power object

## Usage

``` r
# S3 method for class 'ms_power'
summary(object, ...)
```

## Arguments

- object:

  An `ms_power` object.

- ...:

  Not used.

## Value

An object of class `summary.ms_power`: a list with the design overview
(`overview`), the per-look table (`looks`) and the settings
(`settings`); for two endpoints also the familywise error rates
(`fwer`).

## Examples

``` r
ctl <- ms_model(death = log(2) / 12)
trt <- ms_model(death = 0.7 * log(2) / 12)
scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
res <- ms_power(scn, n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
                test = ms_test(), events = 300)
summary(res)
#> Group sequential log-rank design: summary
#>   Tested endpoint: OS; trigger: os; mean: "direct"; variance: "null"
#>   Futility: none; spending time: "events"; alpha: 0.025
#> 
#> Overview:
#>                n accrual_duration   study_duration           events 
#>         400.0000          20.0000          39.6700         300.0000 
#>            power type1_nonbinding    type1_binding 
#>           0.8689           0.0250           0.0250 
#> 
#> By look:
#>  look  time enrolled events events_control events_treatment info_frac z_upper
#>     1 20.47      400    150          84.58            65.42       0.5   2.963
#>     2 39.67      400    300         161.92           138.08       1.0   1.969
#>   p_upper hr_upper z_lower hr_lower z_mean z_sd prob_upper prob_lower cum_power
#>  0.001525   0.6164    -Inf       NA  2.176    1     0.2158     0.0000    0.2158
#>  0.024500   0.7967   1.969   0.7967  3.087    1     0.6531     0.1311    0.8689
#>  alpha_spent cum_alpha
#>     0.001525  0.001525
#>     0.025000  0.025000
```
