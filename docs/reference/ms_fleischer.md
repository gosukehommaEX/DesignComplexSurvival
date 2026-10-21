# Exponential model for PFS and OS of Fleischer et al. (2009)

Creates the `ms_model` in which the time to progression and the time to
death are exponential and the hazard of death does not change at
progression. Then OS is exponential with median `median_os` and PFS is
exponential with median `median_pfs`.

## Usage

``` r
ms_fleischer(median_pfs, median_os, dropout = 0)
```

## Arguments

- median_pfs:

  Median PFS.

- median_os:

  Median OS (not smaller than `median_pfs`).

- dropout:

  Dropout hazard.

## Value

An object of class `ms_model`.

## Details

With death hazard `lambda_d = log(2) / median_os` and PFS hazard
`log(2) / median_pfs`, the progression hazard is
`log(2) / median_pfs - lambda_d`. This is the model of Section 2 of
Fleischer et al. (2009) with the hazard of death before and after
progression equal to `lambda_d`.

## References

Fleischer F, Gaschler-Markefski B, Bluhmki E (2009). A statistical model
for the dependence between progression-free survival and overall
survival. Statistics in Medicine 28:2669-2686.

## See also

[`ms_model`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md)

## Examples

``` r
m <- ms_fleischer(median_pfs = 6, median_os = 15)
ms_surv(m, time = c(6, 15), endpoint = "pfs")
#>   time      surv       haz       dens
#> 1    6 0.5000000 0.1155245 0.05776227
#> 2   15 0.1767767 0.1155245 0.02042204
ms_surv(m, time = c(6, 15), endpoint = "os")
#>   time      surv        haz       dens
#> 1    6 0.7578583 0.04620981 0.03502049
#> 2   15 0.5000000 0.04620981 0.02310491
```
