# Print an ms_test object

Print an ms_test object

## Usage

``` r
# S3 method for class 'ms_test'
print(x, ...)
```

## Arguments

- x:

  An `ms_test` object.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
print(ms_test(futility = "beta"))
#> One-sided log-rank test of OS at level 0.025
#>   Efficacy: alpha spending "OF"
#>   Futility: beta spending "OF" with beta = 0.1 (non-binding)
#>   Mean: "direct"; variance: "null"; spending time: "events"
print(ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence"))
#> One-sided log-rank tests of PFS and OS at familywise level 0.025
#>   Fixed sequence: PFS, then OS (overall)
#>   Efficacy: alpha spending "OF" (PFS), "OF" (OS)
#>   Futility (PFS): none
#>   Mean: "direct"; variance: "null"; spending time: "events"
```
