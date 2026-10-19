
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DesignComplexSurvival

DesignComplexSurvival computes the operating characteristics of
two-group time-to-event trials analysed with the log-rank test, without
patient-level simulation:

- each group follows a four-state model (progression-free, progressed,
  switched, dead) with piecewise constant intensities, covering
  progression-free survival, overall survival and treatment switching at
  or after progression, with the post-progression intensities measured
  from progression or from randomization;
- accrual is piecewise uniform (by rates or by proportions) and dropout
  is piecewise exponential by group;
- expected events, the mean, variance and covariance of the log-rank
  statistic over calendar time, group sequential bounds by alpha and
  beta spending (binding or non-binding futility), power, type I error,
  and the number of events, sample size, accrual duration or follow-up
  for a target power.

This version tests one endpoint (PFS or OS). Designs that test both PFS
and OS with multiplicity adjustment are planned.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("gosukehommaEX/DesignComplexSurvival")
```

## Example

A delayed treatment effect: the hazard ratio is 1 for 6 months and 0.65
afterwards. The number of events for 90 percent power with interim
analyses at 50 and 75 percent of the events:

``` r
library(DesignComplexSurvival)
lam <- log(2) / 12
scn <- ms_scenario(
  control = ms_model(death = lam, dropout = 0.003),
  treatment = ms_model(time = c(0, 6), death = c(lam, 0.65 * lam), dropout = 0.003),
  accrual = ms_accrual(time = c(0, 3), rate = c(10, 30))
)
ms_size(scn, looks = ms_looks("os", fraction = c(0.5, 0.75, 1)),
        test = ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2),
        power = 0.9, solve_for = "events", n = 700)
#> Solved for the final number of events: 587.802 (target power 0.9)
#> 
#> Group sequential log-rank design
#>   Tested endpoint: OS; looks driven by OS events
#>   Sample size: 700; accrual duration: 25.33; final analysis at: 58.66
#>   Power: 0.9001
#>   Type I error: 0.025 (non-binding), 0.02482 (binding)
#>  look  time events info_frac z_upper  p_upper hr_upper z_lower hr_lower
#>     1 25.45    294      0.50   2.963 0.001525   0.7078 -0.5039   1.0605
#>     2 35.51    441      0.75   2.359 0.009162   0.7988  0.7166   0.9340
#>     3 58.66    588      1.00   2.014 0.022000   0.8469  2.0141   0.8469
#>  prob_upper prob_lower cum_power
#>     0.06203    0.02689   0.06203
#>     0.46398    0.02760   0.52601
#>     0.37413    0.04536   0.90014
```

## Validation

The tests compare the results with an independent Python implementation
(`inst/validation/python/`, using matrix exponentials and adaptive
quadrature) and, when the packages are installed, with gsDesign,
gsDesign2, rpact and lrstat.

The scripts in `inst/reproduce/` recompute sample sizes, numbers of
events, analysis times, bounds and powers printed in 10 sources
(articles and a statistical analysis plan). The scripts in
`inst/validation/` compare the package with patient-level simulations of
8 scenarios (FastSurvival and a Python simulation) and with the package
PWEALL in 8 designs. Each folder is run by its `run_all.R`, which judges
every value as PASS (the criterion is met), EXPLAINED (the criterion is
met after a registered explanation, such as an apparent misprint), FAIL,
or INFO (not judged: a published value computed by a method that the
package does not reproduce, or an asymptotic approximation compared with
simulation):

| Folder             | PASS | EXPLAINED | FAIL | INFO | Total |
|:-------------------|-----:|----------:|-----:|-----:|------:|
| `inst/reproduce/`  |  144 |        44 |    0 |  149 |   337 |
| `inst/validation/` |  737 |         0 |    0 |  838 |  1575 |

`run_all.R` also multiplies the compared values by 1.1 and checks that
the judgements change as expected. The criteria are described in the
`README.md` of each folder and the registered explanations in its
`explained.csv`. Every compared value is listed with its judgement in
`results/summary.csv`, and `results/summary.md` summarizes the results.
FastSurvival and PWEALL are needed only to rerun `inst/validation/`.
