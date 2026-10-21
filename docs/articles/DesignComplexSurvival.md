# Designing complex time-to-event trials without simulation

``` r

library(DesignComplexSurvival)
```

DesignComplexSurvival computes the operating characteristics of
two-group time-to-event trials analysed with the log-rank test, without
patient-level simulation. Each group follows a four-state model
(progression-free, progressed, switched, dead) with piecewise constant
intensities. This vignette shows a single endpoint with a delayed
treatment effect, then overall survival (OS) with treatment switching
after progression, and finally progression-free survival (PFS) and OS
tested together.

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
hazards (default) or with the Schoenfeld approximation, which uses the
event-weighted average of the log hazard ratio.

``` r

t_fin <- des$study_duration
rbind(direct = lr_moments(scn, 700, t_fin)$z_mean,
      schoenfeld = lr_moments(scn, 700, t_fin, mean = "schoenfeld")$z_mean)
#>                [,1]
#> direct     3.395909
#> schoenfeld 3.417165
```

## Overall survival with switching after progression

In the control group the hazard of progression is 0.10 per month in the
first 4 months and 0.08 afterwards. At progression, 30 percent of the
patients switch to the experimental treatment; the others may switch
later with a hazard of 0.04 per month. After progression the hazard of
death is 0.06 without switch and 0.035 after a switch. The treatment
group has a lower hazard of progression and no switching.

``` r

ctl_sw <- ms_model(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01,
                   switch_prob = 0.3, death_prog = 0.06, switch_rate = 0.04,
                   death_switch = 0.035, dropout = 0.003)
trt_sw <- ms_model(time = c(0, 4), prog = c(0.06, 0.05), death = 0.008,
                   death_prog = 0.055, dropout = 0.003)
ms_prob(ctl_sw, time = c(6, 12, 24))
#>   time progression_free progressed  switched      dead
#> 1    6        0.5379444  0.2102858 0.1417831 0.1099867
#> 2   12        0.3134862  0.2176410 0.2161661 0.2527068
#> 3   24        0.1064585  0.1329666 0.2539857 0.5065892
```

The post-progression intensities are measured from the time of
progression (`clock = "reset"`, the default). OS and progression-free
survival (PFS) of the control group at 12 months are 0.747 and 0.313.

``` r

scn_sw <- ms_scenario(ctl_sw, trt_sw, accrual = ms_accrual(rate = 25))
pw <- ms_power(scn_sw, n = 500, looks = ms_looks("os", fraction = c(0.6, 1)),
               test = ms_test(), events = 300)
pw
#> Group sequential log-rank design
#>   Tested endpoint: OS; looks driven by OS events
#>   Sample size: 500; accrual duration: 20; final analysis at: 43.89
#>   Power: 0.2206
#>   Type I error: 0.025 (non-binding), 0.025 (binding)
#>  look  time events info_frac z_upper  p_upper hr_upper z_lower hr_lower
#>     1 28.82    180       0.6   2.669 0.003808   0.6718    -Inf       NA
#>     2 43.89    300       1.0   1.981 0.023798   0.7955   1.981   0.7955
#>  prob_upper prob_lower cum_power
#>     0.08221     0.0000   0.08221
#>     0.13839     0.7794   0.22059
```

With 300 OS events the power is 0.221. Switching to the experimental
treatment in the control group moves the hazard ratio for OS towards 1.
For comparison, the same design without switching (all switching
intensities set to zero) has the power

``` r

ctl_ns <- ms_model(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01,
                   death_prog = 0.06, dropout = 0.003)
scn_ns <- ms_scenario(ctl_ns, trt_sw, accrual = ms_accrual(rate = 25))
pw_ns <- ms_power(scn_ns, n = 500, looks = ms_looks("os", fraction = c(0.6, 1)),
                  test = ms_test(), events = 300)
pw_ns$power
#> [1] 0.8426871
```

## PFS and OS tested together

Both groups follow the model of Fleischer et al. (2009): the times to
progression and to death are exponential and the hazard of death does
not change at progression, so that PFS and OS are exponential. The
medians of PFS and OS are 4 and 12 months in the control group and 6 and
16 months in the treatment group, and patients are enrolled at 25 per
month.

PFS is tested at 300 PFS events, where OS is also tested, and OS is
tested again at 220 OS events. With the fixed-sequence procedure, both
hypotheses are tested at level 0.025 and OS is tested only after PFS has
been rejected (procedure (b) of Glimm et al. 2010). Each endpoint has
its own efficacy bounds of O’Brien-Fleming type.

``` r

fl <- ms_scenario(ms_fleischer(median_pfs = 4, median_os = 12),
                  ms_fleischer(median_pfs = 6, median_os = 16),
                  accrual = ms_accrual(rate = 25))
lk2 <- ms_looks(c("pfs", "os"), fraction = c(1, 1),
                test = list(c("pfs", "os"), "os"))
fs <- ms_power(fl, n = 450, looks = lk2,
               test = ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence"),
               events = c(pfs = 300, os = 220))
fs
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

The power is 0.939 for PFS and 0.545 for OS, and both hypotheses are
rejected with probability 0.545. The final OS analysis takes place at
month 23.1.

The Bonferroni procedure tests PFS at level 0.02 and OS at level 0.005,
whatever happens to the other hypothesis.

``` r

bf <- ms_power(fl, n = 450, looks = lk2,
               test = ms_test(endpoints = c("pfs", "os"), mtp = "bonferroni",
                              weights = c(0.8, 0.2)),
               events = c(pfs = 300, os = 220))
rbind(fixed_sequence = fs$power, bonferroni = bf$power)
#>                      pfs        os      both       any
#> fixed_sequence 0.9387378 0.5451639 0.5451639 0.9387378
#> bonferroni     0.9265375 0.3237669 0.3173922 0.9329122
```

The familywise error rate is computed under the global null hypothesis,
where the correlation between the PFS and OS statistics is that of
Danzer et al. (2025). A partial null hypothesis, in which OS is the same
in both groups but PFS is not, can be given through `null_scenarios`:

``` r

os_null <- ms_scenario(ms_fleischer(median_pfs = 4, median_os = 12),
                       ms_fleischer(median_pfs = 6, median_os = 12),
                       accrual = ms_accrual(rate = 25))
ms_power(fl, n = 450, looks = lk2,
         test = ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence"),
         events = c(pfs = 300, os = 220), null_scenarios = list(os = os_null))$fwer
#>              scenario null_endpoints nonbinding    binding
#> 1              global        PFS, OS 0.02500000 0.02500000
#> 2 null_scenarios[[1]]             OS 0.02405394 0.02405394
```

## Choice of the covariance

The argument `variance` of
[`ms_test()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md)
(and of
[`lr_moments()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_moments.md)
and
[`lr_cov()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md))
chooses the covariance of the standardized statistics. By default
(`"null"`) the statistics of one endpoint have unit variances and the
correlations of independent increments. With `"alternative"`, the
variance and covariance of the score under the alternative hypothesis
(Luo et al. 2019) are used instead. With `"studentized"`, the covariance
is the asymptotic covariance of the statistic that is computed in
practice, standardized by the estimated variance and, for event-driven
analyses, taken at a random analysis time.

``` r

mo <- lapply(c(alternative = "alternative", studentized = "studentized"),
             function(v) lr_moments(scn, n = 700, calendar_time = c(24, t_fin),
                                    variance = v))
sapply(mo, function(m) m$z_var)
#>      alternative studentized
#> [1,]   0.9985616   1.0004021
#> [2,]   0.9530797   0.9900072
```

Here the largest difference of the variances from 1 is 0.047 with
`"alternative"` and 0.01 with `"studentized"`. The choice affects the
power and the futility bounds by beta spending, but not the efficacy
bounds or the non-binding type I error. In the simulations of
`inst/validation/` with PFS and OS tested together, `"studentized"`
agreed best with the simulated standard deviations and correlations of
the statistics, and the powers with `"null"` and `"studentized"` agreed
with the simulated powers within the Monte Carlo error.

## References

Danzer MF, Rufibach K, Beyersmann J, Schmidt R (2025). Exhausting the
type I error level in event-driven group-sequential designs with a
closed testing procedure for progression-free and overall survival.
arXiv:2512.08658.

Fleischer F, Gaschler-Markefski B, Bluhmki E (2009). A statistical model
for the dependence between progression-free survival and overall
survival. Statistics in Medicine 28:2669-2686.

Glimm E, Maurer W, Bretz F (2010). Hierarchical testing of multiple
endpoints in group-sequential trials. Statistics in Medicine 29:219-228.

Luo X, Mao X, Chen X, Qiu J, Bai S, Quan H (2019). Design and monitoring
of survival trials in complex scenarios. Statistics in Medicine
38:192-209.
