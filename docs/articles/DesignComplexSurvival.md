# Designing complex time-to-event trials without simulation

``` r

library(DesignComplexSurvival)
```

DesignComplexSurvival computes the operating characteristics of
two-group time-to-event trials analysed with the log-rank test, without
patient-level simulation. This article first describes what the package
computes and which trial settings it covers, and then goes through a
design with a single endpoint. The other articles are listed at the end.

## What the package computes

- The expected numbers of enrolled patients, events and dropouts over
  calendar time.
- The calendar times at which given numbers of events are expected, that
  is, the times of the interim and final analyses.
- The efficacy and futility bounds at each analysis and the
  probabilities of crossing them.
- The power and the type I error; with two endpoints, the power of each,
  the probabilities of rejecting both or at least one hypothesis, and
  the familywise error rate.
- The number of events, sample size, accrual duration or follow-up
  needed for a target power.

The analyses can be driven by the events of progression-free survival
(PFS) or overall survival (OS), take place at calendar times, or mix the
two. The computations use closed forms, numerical integration and
multivariate normal probabilities, so that they are fast and give the
same result every time; they rest on the normal approximation of the
log-rank statistics.

## The model of each group

Each group follows a multistate model with five states:

![](data:image/svg+xml;base64,PHN2ZyB2aWV3Ym94PSIwIDAgNzIwIDQyMCIgd2lkdGg9IjEwMCUiIHN0eWxlPSJtYXgtd2lkdGg6IDcyMHB4OyBoZWlnaHQ6IGF1dG87IGNvbG9yOiBpbmhlcml0OyIgcm9sZT0iaW1nIiBhcmlhLWxhYmVsPSJUaGUgZml2ZSBzdGF0ZXMgb2YgdGhlIG1vZGVsIG9mIG9uZSBncm91cCBhbmQgdGhlIHRyYW5zaXRpb25zIGJldHdlZW4gdGhlbSI+PGRlZnM+PG1hcmtlciBpZD0iZGNzLWFycm93IiB2aWV3Ym94PSIwIDAgMTAgMTAiIHJlZng9IjkiIHJlZnk9IjUiIG1hcmtlcndpZHRoPSI3IiBtYXJrZXJoZWlnaHQ9IjciIG9yaWVudD0iYXV0by1zdGFydC1yZXZlcnNlIj48cGF0aCBkPSJNMCwwIEwxMCw1IEwwLDEwIHoiIGZpbGw9ImN1cnJlbnRDb2xvciIgZmlsbC1vcGFjaXR5PSIwLjciIC8+PC9tYXJrZXI+PC9kZWZzPjxnIGZpbGw9Im5vbmUiIHN0cm9rZT0iY3VycmVudENvbG9yIiBzdHJva2Utb3BhY2l0eT0iMC42IiBzdHJva2Utd2lkdGg9IjEuNCI+PHBhdGggZD0iTTk1LDU4IEM5NSwyMCA2MjUsMjAgNjI1LDU4IiBtYXJrZXItZW5kPSJ1cmwoI2Rjcy1hcnJvdykiIC8+PGxpbmUgeDE9IjE4MCIgeTE9IjkwIiB4Mj0iMjczIiB5Mj0iOTAiIG1hcmtlci1lbmQ9InVybCgjZGNzLWFycm93KSI+PC9saW5lPjxsaW5lIHgxPSI0NDUiIHkxPSI5MCIgeDI9IjUzOCIgeTI9IjkwIiBtYXJrZXItZW5kPSJ1cmwoI2Rjcy1hcnJvdykiPjwvbGluZT48bGluZSB4MT0iMTMwIiB5MT0iMTIwIiB4Mj0iMzE4IiB5Mj0iMzA4IiBtYXJrZXItZW5kPSJ1cmwoI2Rjcy1hcnJvdykiPjwvbGluZT48bGluZSB4MT0iMzYwIiB5MT0iMTIwIiB4Mj0iMzYwIiB5Mj0iMzA4IiBtYXJrZXItZW5kPSJ1cmwoI2Rjcy1hcnJvdykiPjwvbGluZT48bGluZSB4MT0iNjEwIiB5MT0iMTIwIiB4Mj0iNDEwIiB5Mj0iMzA4IiBtYXJrZXItZW5kPSJ1cmwoI2Rjcy1hcnJvdykiPjwvbGluZT48bGluZSB4MT0iOTUiIHkxPSIxMjAiIHgyPSI5NSIgeTI9IjE5OCIgbWFya2VyLWVuZD0idXJsKCNkY3MtYXJyb3cpIj48L2xpbmU+PGxpbmUgeDE9IjE4MCIgeTE9IjI0NSIgeDI9IjI3MyIgeTI9IjMzNiIgbWFya2VyLWVuZD0idXJsKCNkY3MtYXJyb3cpIj48L2xpbmU+PHBhdGggZD0iTTk1LDI2MiBWNDAwIEg3MTUgVjkwIEg3MTIiIG1hcmtlci1lbmQ9InVybCgjZGNzLWFycm93KSIgLz48L2c+PGcgZmlsbD0iY3VycmVudENvbG9yIiBmaWxsLW9wYWNpdHk9IjAuNzUiIGZvbnQtZmFtaWx5PSJtb25vc3BhY2UiIGZvbnQtc2l6ZT0iMTEuNSI+PHRleHQgeD0iMzYwIiB5PSIxNCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+c3dpdGNoIGF0IHByb2dyZXNzaW9uIChwcm9nLCBzd2l0Y2hfcHJvYik8L3RleHQ+PHRleHQgeD0iMjI2IiB5PSI4MiIgdGV4dC1hbmNob3I9Im1pZGRsZSI+cHJvZzwvdGV4dD48dGV4dCB4PSI0OTEiIHk9IjExMCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+c3dpdGNoX3JhdGU8L3RleHQ+PHRleHQgeD0iMjIyIiB5PSIxOTYiPmRlYXRoPC90ZXh0Pjx0ZXh0IHg9IjM2OCIgeT0iMjI2Ij5kZWF0aF9wcm9nPC90ZXh0Pjx0ZXh0IHg9IjUxNiIgeT0iMjE0Ij5kZWF0aF9zd2l0Y2g8L3RleHQ+PHRleHQgeD0iMTAzIiB5PSIxNjQiPnN3aXRjaF9wcmU8L3RleHQ+PHRleHQgeD0iMjAwIiB5PSIzMTgiIHRleHQtYW5jaG9yPSJlbmQiPmRlYXRoX3N3aXRjaF9wcmU8L3RleHQ+PHRleHQgeD0iNDIwIiB5PSIzOTIiIHRleHQtYW5jaG9yPSJtaWRkbGUiPnByb2dyZXNzaW9uIGFmdGVyIGEgc3dpdGNoIChwcm9nX3N3aXRjaCk8L3RleHQ+PC9nPjxnIGZpbGw9Im5vbmUiIHN0cm9rZT0iY3VycmVudENvbG9yIiBzdHJva2Utd2lkdGg9IjEuNiI+PHJlY3QgeD0iMTAiIHk9IjYwIiB3aWR0aD0iMTcwIiBoZWlnaHQ9IjYwIiByeD0iOCIgLz48cmVjdCB4PSIyNzUiIHk9IjYwIiB3aWR0aD0iMTcwIiBoZWlnaHQ9IjYwIiByeD0iOCIgLz48cmVjdCB4PSI1NDAiIHk9IjYwIiB3aWR0aD0iMTcwIiBoZWlnaHQ9IjYwIiByeD0iOCIgLz48cmVjdCB4PSIxMCIgeT0iMjAwIiB3aWR0aD0iMTcwIiBoZWlnaHQ9IjYwIiByeD0iOCIgLz48cmVjdCB4PSIyNzUiIHk9IjMxMCIgd2lkdGg9IjE3MCIgaGVpZ2h0PSI2MCIgcng9IjgiIGZpbGw9ImN1cnJlbnRDb2xvciIgZmlsbC1vcGFjaXR5PSIwLjA4IiAvPjwvZz48ZyBmaWxsPSJjdXJyZW50Q29sb3IiIGZvbnQtZmFtaWx5PSJzYW5zLXNlcmlmIiBmb250LXNpemU9IjEzIj48dGV4dCB4PSIyMCIgeT0iNzgiIGZvbnQtZmFtaWx5PSJtb25vc3BhY2UiIGZvbnQtc2l6ZT0iMTIiIGZpbGwtb3BhY2l0eT0iMC43Ij4xPC90ZXh0Pjx0ZXh0IHg9Ijk1IiB5PSI5NiIgdGV4dC1hbmNob3I9Im1pZGRsZSI+UHJvZ3Jlc3Npb24tZnJlZTwvdGV4dD48dGV4dCB4PSIyODUiIHk9Ijc4IiBmb250LWZhbWlseT0ibW9ub3NwYWNlIiBmb250LXNpemU9IjEyIiBmaWxsLW9wYWNpdHk9IjAuNyI+MjwvdGV4dD48dGV4dCB4PSIzNjAiIHk9Ijg4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5Qcm9ncmVzc2VkPC90ZXh0Pjx0ZXh0IHg9IjM2MCIgeT0iMTA2IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj4obm8gc3dpdGNoKTwvdGV4dD48dGV4dCB4PSI1NTAiIHk9Ijc4IiBmb250LWZhbWlseT0ibW9ub3NwYWNlIiBmb250LXNpemU9IjEyIiBmaWxsLW9wYWNpdHk9IjAuNyI+MzwvdGV4dD48dGV4dCB4PSI2MjUiIHk9Ijk2IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5Td2l0Y2hlZDwvdGV4dD48dGV4dCB4PSIyMCIgeT0iMjE4IiBmb250LWZhbWlseT0ibW9ub3NwYWNlIiBmb250LXNpemU9IjEyIiBmaWxsLW9wYWNpdHk9IjAuNyI+NTwvdGV4dD48dGV4dCB4PSI5NSIgeT0iMjI4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5Qcm9ncmVzc2lvbi1mcmVlPC90ZXh0Pjx0ZXh0IHg9Ijk1IiB5PSIyNDYiIHRleHQtYW5jaG9yPSJtaWRkbGUiPmFmdGVyIGEgc3dpdGNoPC90ZXh0Pjx0ZXh0IHg9IjI4NSIgeT0iMzI4IiBmb250LWZhbWlseT0ibW9ub3NwYWNlIiBmb250LXNpemU9IjEyIiBmaWxsLW9wYWNpdHk9IjAuNyI+NDwvdGV4dD48dGV4dCB4PSIzNjAiIHk9IjM0NiIgdGV4dC1hbmNob3I9Im1pZGRsZSI+RGVhZDwvdGV4dD48L2c+PC9zdmc+)

The states of the model of one group and the arguments of
[`ms_model()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md)
that give the transitions between them. State 5 is used only when
`switch_pre` is positive.

PFS is the time to leave the progression-free states 1 and 5, and OS is
the time to enter state 4. Since both come from the same course of each
patient, the dependence between them is part of the model and needs no
separate assumption. A switch (also called treatment switching or
crossover) is a move to the experimental treatment or to another
treatment; it can be specified in either group:

- at progression, with probability `switch_prob`, or after progression,
  with hazard `switch_rate`;
- before progression, with hazard `switch_pre`; PFS is then followed up
  to progression or death, or censored at the switch (`pfs_switch` in
  [`ms_scenario()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md));
- from a calendar time on, for example after an interim analysis
  (`switch_start` in
  [`ms_scenario()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md)).

All intensities are piecewise constant, and the number of intervals is
not limited, so that a delayed or diminishing treatment effect is given
through the intervals. Before progression the intensities are functions
of the time since randomization. After progression they are functions of
the time since progression (`clock = "reset"`) or since randomization
(`"forward"`). The death hazard after a switch can be measured on the
same clock or from the switch (`switch_clock = "switch"`), or the switch
can have an accelerated failure time effect (`switch_af`). Dropout is
piecewise exponential by group, the same in all states, and
`dropout_pfs` adds a censoring of PFS only.

## What can be specified

| Feature | Options | Status |
|----|----|----|
| Accrual | piecewise uniform, by rates or by proportions and a duration | supported |
| Allocation | any ratio, for example 2:1 | supported |
| Dropout | piecewise exponential by group; a censoring of PFS only | supported |
| Hazards | piecewise constant, with delayed or diminishing effects through the intervals | supported |
| Switching at and after progression | probability at progression, hazard after progression, death hazard after the switch | supported |
| Effect of a switch | death hazard on the post-progression clock or from the switch, or an accelerated failure time effect | supported |
| Switching before progression | PFS followed up or censored at the switch | supported |
| Switching from a calendar time | OS tested alone | partly |
| Timing of the analyses | events of OS or PFS, calendar times, or a mix; final analysis at a given follow-up after the end of accrual | supported |
| Efficacy bounds | alpha spending: O’Brien-Fleming and Pocock types, Hwang-Shih-DeCani, power family, user-defined | supported |
| Futility bounds | fixed (Z or hazard ratio scale) or beta spending; binding or non-binding | supported |
| Two endpoints | PFS and OS by fixed sequence (two procedures) or Bonferroni, with recycling of the level after a rejection | supported |
| Familywise error rate | under the global null hypothesis and under partial null scenarios given by the user | supported |
| Solving for | number of events, sample size, accrual duration or follow-up | supported |
| Approximations | mean of the statistic (two choices) and covariance (three choices) | supported |
| Smooth hazards (Weibull), cure fractions, delays that vary between patients | piecewise constant approximation of the hazard of each group (one endpoint) | approximate |
| Not covered | weighted log-rank tests, restricted mean survival time, stratification, non-inferiority and equivalence, a fixed follow-up for each patient, repeated switches back and forth, dropout that depends on the state, switching started by interim results, more than two groups | not supported |

PFS and OS can be tested together with switching and with group
sequential designs; the exception is switching from a calendar time,
which covers the OS statistic alone. The other articles describe each
part in detail.

## Building a design

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

The following example goes through the four steps for a single endpoint
with a delayed treatment effect.

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
