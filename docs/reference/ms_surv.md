# Survival function, hazard and density of an endpoint

Evaluates the distribution of progression-free survival or overall
survival implied by an `ms_model`, without dropout.

## Usage

``` r
ms_surv(model, time, endpoint = c("os", "pfs"))
```

## Arguments

- model:

  An `ms_model` object.

- time:

  Times since randomization (nonnegative).

- endpoint:

  `"os"` or `"pfs"`.

## Value

A data frame with `time`, the survival function `surv`, the hazard `haz`
(`NA` where `surv` is 0) and the density `dens`.

## Details

PFS has the hazard `prog + death`. OS combines death before progression
with death after progression, with or without switch; its hazard is not
piecewise constant when `prog` is positive. With a switch before
progression (`switch_pre`), PFS is the time to progression or death
whatever `pfs_switch` of
[`ms_scenario`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md),
and its hazard is not piecewise constant.

## See also

[`ms_model`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md),
[`ms_prob`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_prob.md)

## Examples

``` r
m <- ms_fleischer(median_pfs = 6, median_os = 15)
ms_surv(m, time = c(6, 15), endpoint = "os")
#>   time      surv        haz       dens
#> 1    6 0.7578583 0.04620981 0.03502049
#> 2   15 0.5000000 0.04620981 0.02310491
```
