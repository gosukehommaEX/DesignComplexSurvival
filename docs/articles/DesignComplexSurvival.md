# Designing complex time-to-event trials without simulation

``` r

library(DesignComplexSurvival)
```

DesignComplexSurvival computes the operating characteristics of
two-group time-to-event trials analysed with the log-rank test, without
patient-level simulation. Each group follows a multistate model
(progression-free, progressed, switched, dead, and progression-free
after a switch before progression) with piecewise constant intensities.
A design is built in four steps:

1.  the scenario: the models of the two groups
    ([`ms_model()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md),
    or
    [`ms_fleischer()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_fleischer.md)
    for exponential progression-free survival (PFS) and overall survival
    (OS)), the allocation and the accrual
    ([`ms_accrual()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_accrual.md),
    [`ms_scenario()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md));
2.  the timing of the analyses
    ([`ms_looks()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_looks.md));
3.  the test: the endpoints, the spending functions, the futility bounds
    and the approximations of the distribution of the statistics
    ([`ms_test()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md));
4.  the evaluation: the power and type I error of a given design
    ([`ms_power()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)),
    or the number of events, sample size, accrual duration or follow-up
    for a target power
    ([`ms_size()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md)).

This article goes through the four steps for a single endpoint with a
delayed treatment effect. The other articles are listed at the end.

## A delayed treatment effect

The control group has exponential survival with a median of 12 months.
In the treatment group the hazard ratio is 1 for the first 6 months and
0.65 afterwards. Patients are enrolled at 10 per month for 3 months and
at 30 per month afterwards, and the dropout hazard is 0.003 per month in
both groups.

``` r

lam <- log(2) / 12
ctl <- ms_model(death = lam, dropout = 0.003)
trt <- ms_model(time = c(0, 6), death = c(lam, 0.65 * lam), dropout = 0.003)
acc <- ms_accrual(time = c(0, 3), rate = c(10, 30))
scn <- ms_scenario(ctl, trt, allocation = 1, accrual = acc)
ev <- ms_events(scn, n = 700, calendar_time = c(12, 24, 36))
ev[, c("time", "enrolled", "events_control", "events_treatment", "events")]
#>   time enrolled events_control events_treatment    events
#> 1   12      300       36.19604         34.58153  70.77756
#> 2   24      660      140.90234        124.55344 265.45579
#> 3   36      700      239.82054        206.40600 446.22654
```

With 700 patients, 446 events are expected by month 36.

The design has interim analyses at 50 and 75 percent of the events, an
efficacy bound of O’Brien-Fleming type and a non-binding futility bound
by beta spending.
[`ms_size()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md)
finds the number of events for a power of 90 percent.

``` r

looks <- ms_looks("os", fraction = c(0.5, 0.75, 1))
test <- ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2)
des <- ms_size(scn, looks = looks, test = test, power = 0.9,
               solve_for = "events", n = 700)
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

The accrual ends at month 25.3. The final analysis needs 588 events and
takes place at month 58.7. The type I error is 0.025 when the futility
bounds are ignored.

The mean of the standardized statistic can be computed directly from the
hazards (the default, following Yung and Liu 2020) or with the
approximation of Schoenfeld (1981), which uses the event-weighted
average of the log hazard ratio.

``` r

t_fin <- des$study_duration
rbind(direct = lr_moments(scn, 700, t_fin)$z_mean,
      schoenfeld = lr_moments(scn, 700, t_fin, mean = "schoenfeld")$z_mean)
#>                [,1]
#> direct     3.395909
#> schoenfeld 3.417165
```

## Summaries of a design

[`summary()`](https://rdrr.io/r/base/summary.html) gives an overview of
the design, and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) the table
of the looks:

``` r

summary(des)
#> Group sequential log-rank design: summary
#>   Tested endpoint: OS; trigger: os; mean: "direct"; variance: "null"
#>   Futility: beta (non-binding); spending time: "events"; alpha: 0.025
#> 
#> Overview:
#>                n accrual_duration   study_duration           events 
#>        700.00000         25.33000         58.66000        588.00000 
#>            power type1_nonbinding    type1_binding 
#>          0.90010          0.02500          0.02482 
#> 
#> By look:
#>  look  time enrolled events events_control events_treatment info_frac z_upper
#>     1 25.45      700    294          156.3            137.7      0.50   2.963
#>     2 35.51      700    441          237.0            204.0      0.75   2.359
#>     3 58.66      700    588          309.3            278.7      1.00   2.014
#>   p_upper hr_upper z_lower hr_lower z_mean z_sd prob_upper prob_lower cum_power
#>  0.001525   0.7078 -0.5039   1.0605  1.425    1    0.06203    0.02689   0.06203
#>  0.009162   0.7988  0.7166   0.9340  2.424    1    0.46398    0.02760   0.52601
#>  0.022000   0.8469  2.0141   0.8469  3.396    1    0.37413    0.04536   0.90014
#>  alpha_spent cum_alpha
#>     0.001525  0.001525
#>     0.009649  0.009649
#>     0.025000  0.025000
as.data.frame(des)[, c("time", "events", "z_upper", "hr_upper", "z_lower",
                       "cum_power")]
#>       time events  z_upper  hr_upper    z_lower  cum_power
#> 1 25.44525    294 2.962588 0.7078215 -0.5039193 0.06202802
#> 2 35.50958    441 2.359018 0.7987811  0.7165938 0.52601131
#> 3 58.66078    588 2.014084 0.8469456  2.0140837 0.90014124
```

[`ms_power()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)
evaluates a design with a given final number of events, a given
follow-up after the end of accrual, or analyses at calendar times. With
the same looks and test and the final analysis 12 months after the end
of accrual:

``` r

pw12 <- ms_power(scn, n = 700, looks = looks, test = test, followup = 12)
c(power = pw12$power, final_analysis = pw12$study_duration,
  events = pw12$events)
#>          power final_analysis         events 
#>      0.7066646     37.3333333    459.7884476
```

The final analysis then takes place at month 37.3 with 460 expected
events, and the power is 0.707.

## Further articles

- [`vignette("switching")`](https://gosukehommaEX.github.io/DesignComplexSurvival/articles/switching.md):
  the multistate model, switching at, after and before progression, the
  clocks of the post-progression intensities, and the censoring of PFS
  only.
- [`vignette("group-sequential")`](https://gosukehommaEX.github.io/DesignComplexSurvival/articles/group-sequential.md):
  the timing of the analyses, the spending functions, the futility
  bounds and the approximations of the covariance of the statistics.
- [`vignette("two-endpoints")`](https://gosukehommaEX.github.io/DesignComplexSurvival/articles/two-endpoints.md):
  PFS and OS tested together with a fixed-sequence or Bonferroni
  procedure, recycling of the level and the familywise error rate.
- [`vignette("calendar-switching")`](https://gosukehommaEX.github.io/DesignComplexSurvival/articles/calendar-switching.md):
  switching that starts at a calendar time, for example after an interim
  analysis.
- [`vignette("validation")`](https://gosukehommaEX.github.io/DesignComplexSurvival/articles/validation.md):
  the reproduction of published results and the comparison with
  patient-level simulations.

## References

Schoenfeld D (1981). The asymptotic properties of nonparametric tests
for comparing survival distributions. Biometrika 68:316-319.

Yung G, Liu Y (2020). Sample size and power for the weighted log-rank
test and Kaplan-Meier based tests with allowance for nonproportional
hazards. Biometrics 76:939-950.
