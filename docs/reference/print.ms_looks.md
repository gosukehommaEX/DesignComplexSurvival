# Print an ms_looks object

Print an ms_looks object

## Usage

``` r
# S3 method for class 'ms_looks'
print(x, ...)
```

## Arguments

- x:

  An `ms_looks` object.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
print(ms_looks("os", fraction = c(0.5, 1)))
#> 2 look(s) driven by OS events at fractions: 0.5, 1.0
print(ms_looks(c("pfs", "os"), fraction = c(1, 1), test = list(c("pfs", "os"), "os")))
#> 2 look(s):
#>   1: driven by PFS events at fraction 1; tests PFS and OS
#>   2: driven by OS events at fraction 1; tests OS
```
