# DesignComplexSurvival

DesignComplexSurvival computes the operating characteristics of
two-group time-to-event trials analyzed with the log-rank test, without
patient-level simulation:

- each group follows a multistate model (progression-free, progressed,
  switched, dead, and progression-free after a switch before
  progression) with piecewise constant intensities, covering
  progression-free survival, overall survival and treatment switching
  before, at or after progression, with the post-progression intensities
  measured from progression or from randomization, and the effect of a
  switch on the death hazard (measured on the post-progression clock or
  from the switch) or as an accelerated failure time effect;
- switching may start at a calendar time, for example after an interim
  analysis (for the OS statistic);
- accrual is piecewise uniform (by rates or by proportions), dropout is
  piecewise exponential by group, and PFS may have a censoring of its
  own and be censored at a switch before progression;
- expected events, the mean, variance and covariance of the log-rank
  statistic over calendar time, group sequential bounds by alpha and
  beta spending (binding or non-binding futility), power, type I error,
  and the number of events, sample size, accrual duration or follow-up
  for a target power;
- progression-free survival (PFS) and overall survival (OS) can be
  tested together with a group sequential fixed-sequence or Bonferroni
  procedure, the latter with recycling of the level after a rejection
  (Maurer and Bretz 2013), with analyses driven by the events of either
  endpoint or by calendar time; the joint distribution of the two
  log-rank statistics gives the power of each endpoint, the probability
  of rejecting both or at least one hypothesis, and the familywise error
  rate;
- for a design with one endpoint, the expected study duration and
  numbers of events and enrolled patients when the trial stops at its
  bounds, and the conditional power at an interim analysis under the
  design scenario, the null hypothesis, the current trend, a hazard
  ratio or another scenario;
- plots of the bounds, of the cumulative probabilities of stopping and
  of the survival functions of the groups.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("gosukehommaEX/DesignComplexSurvival")
```

## Documentation

The articles on the package website
(<https://gosukehommaEX.github.io/DesignComplexSurvival/>) cover the
workflow of a design, the multistate model and treatment switching,
group sequential designs (with the expected duration, the conditional
power and the plots of a design) and the covariance of the statistics,
PFS and OS tested together, switching from a calendar time, and the
validation of the package.

## Examples

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
des <- ms_size(scn, looks = ms_looks("os", fraction = c(0.5, 0.75, 1)),
               test = ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2),
               power = 0.9, solve_for = "events", n = 700)
des
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

The expected study duration, events and enrolled patients when the trial
stops at the bounds, and the conditional power at the first interim
analysis with the observed statistic 1.5 under the current trend:

``` r

des$expected
#>   hypothesis prob_efficacy_early prob_futility_early study_duration   events
#> 1     design          0.52601131          0.05449458       44.32643 489.5941
#> 2       null          0.00964929          0.76522251       34.26686 428.7172
#>   enrolled
#> 1 700.0000
#> 2 692.9973
ms_cond_power(des, look = 1, z = 1.5, under = "trend")
#> Conditional power of a group sequential log-rank design
#>   Observed at look 1: z = 1.5
#>   Distribution: current trend (hazard ratio 0.8395)
#>   Conditional power: 0.5711
#>   Without the futility bounds of the later interim looks: 0.5717
#>  look  time events z_upper z_lower cond_mean cond_sd prob_upper prob_lower
#>     2 35.51    441   2.359  0.7166     1.837  0.5774     0.1830    0.02614
#>     3 58.66    588   2.014  2.0141     2.121  0.7071     0.3881    0.40272
#>  cum_upper
#>     0.1830
#>     0.5711
```

PFS and OS tested together: PFS and OS are exponential in both groups
(the model of Fleischer et al. 2009), with medians of 4 and 12 months in
the control group and 6 and 16 months in the treatment group. PFS is
tested at 300 PFS events, where OS is also tested, and OS again at 220
OS events. OS is tested only after PFS has been rejected (fixed
sequence):

``` r

fl <- ms_scenario(
  control = ms_fleischer(median_pfs = 4, median_os = 12),
  treatment = ms_fleischer(median_pfs = 6, median_os = 16),
  accrual = ms_accrual(rate = 25)
)
ms_power(fl, n = 450,
         looks = ms_looks(c("pfs", "os"), fraction = c(1, 1),
                          test = list(c("pfs", "os"), "os")),
         test = ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence"),
         events = c(pfs = 300, os = 220))
#> Group sequential log-rank design with two endpoints
#>   Endpoints: PFS, then OS (fixed sequence, overall); looks driven by PFS events and OS events
#>   Sample size: 450; accrual duration: 18; final analysis at: 23.05
#>   Power: PFS 0.9387, OS 0.5452; both 0.5452; at least one 0.9387
#>   FWER under global null: 0.025 (non-binding), 0.025 (binding)
#>  endpoint look  time events info_frac z_upper  p_upper hr_upper z_lower
#>       PFS    1 18.62    300    1.0000   1.960 0.025000   0.7975   1.960
#>        OS    1 18.62    163    0.7407   2.357 0.009207   0.6912    -Inf
#>        OS    2 23.05    220    1.0000   2.009 0.022243   0.7626   2.009
#>  hr_lower prob_upper prob_lower cum_power
#>    0.7975     0.9387          0    0.9387
#>        NA     0.2966          0    0.2966
#>    0.7626     0.2486          0    0.5452
```

## Validation

The tests compare the results with an independent Python implementation
(`inst/validation/python/`, using matrix exponentials and adaptive
quadrature) and, when the packages are installed, with gsDesign,
gsDesign2, rpact and lrstat.

The scripts in `inst/reproduce/` recompute sample sizes, numbers of
events, analysis times, bounds, powers, expected study durations and
conditional powers printed in 18 sources (articles and a statistical
analysis plan). The scripts in `inst/validation/` compare the package
with patient-level simulations of 15 scenarios (FastSurvival and a
Python simulation; PFS and OS analyzed together in 6 of them, switching
from a calendar time in 3), with the covariance estimator of Danzer et
al. (2025) applied to large simulated trials under the null hypothesis,
with the package PWEALL in 8 designs, with the decisions of
[`lrstat::fseqbon()`](https://kaifenglu.github.io/lrstat/reference/fseqbon.html)
in 3 designs with recycling of the level, with rpact, lrstat and
gsDesign for the expected values when the trial stops at the bounds in 9
designs, with
[`gsDesign::gsCP()`](https://keaven.github.io/gsDesign//reference/gsCP.html)
for the conditional power in 4 designs, with the values printed in the
vignettes of rpact and gsDesign, and with a patient-level simulation of
the expected values in 3 designs. Each folder is run by its `run_all.R`,
which judges every value as PASS (the criterion is met), EXPLAINED (the
criterion is met after a registered explanation, such as an apparent
misprint), FAIL, or INFO (not judged: a published value computed by a
method that the package does not reproduce, or an asymptotic
approximation compared with simulation):

| Folder             | PASS | EXPLAINED | FAIL | INFO | Total |
|:-------------------|-----:|----------:|-----:|-----:|------:|
| `inst/reproduce/`  |  485 |        91 |    0 |  204 |   780 |
| `inst/validation/` | 1726 |         3 |    0 | 2217 |  3946 |

`run_all.R` also multiplies the compared values by 1.1 and checks that
the judgments change as expected. The criteria are described in the
`README.md` of each folder and the registered explanations in its
`explained.csv`. Every compared value is listed with its judgment in
`results/summary.csv`, and `results/summary.md` summarizes the results.
FastSurvival and PWEALL are needed only to rerun `inst/validation/`,
which also uses lrstat (0.3.4 or later), rpact and gsDesign.
