# Print an ms_accrual object

Print an ms_accrual object

## Usage

``` r
# S3 method for class 'ms_accrual'
print(x, ...)
```

## Arguments

- x:

  An `ms_accrual` object.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
print(ms_accrual(time = c(0, 3), rate = c(10, 30)))
#> Piecewise uniform accrual by rates (patients per unit time):
#>  time rate
#>     0   10
#>     3   30
```
