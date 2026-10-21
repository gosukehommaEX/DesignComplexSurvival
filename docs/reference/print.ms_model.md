# Print an ms_model object

Print an ms_model object

## Usage

``` r
# S3 method for class 'ms_model'
print(x, ...)
```

## Arguments

- x:

  An `ms_model` object.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
print(ms_fleischer(median_pfs = 6, median_os = 15))
#> Four-state model for one group (clock after progression: "reset")
#> Before progression (time since randomization):
#>  time       prog      death switch_prob dropout
#>     0 0.06931472 0.04620981           0       0
#> After progression (time since progression):
#>  time_post death_prog switch_rate death_switch
#>          0 0.04620981           0   0.04620981
```
