# PFS and OS tested together

``` r

library(DesignComplexSurvival)
```

Progression-free survival (PFS) and overall survival (OS) come from the
same multistate model, so that their log-rank statistics are correlated.
The package computes the joint distribution of the two statistics over
calendar time and tests the two hypotheses with a group sequential
fixed-sequence or Bonferroni procedure, the latter with recycling of the
level after a rejection.

## Scenario and analyses

Both groups follow the model of Fleischer et al. (2009): the times to
progression and to death are exponential and the hazard of death does
not change at progression, so that PFS and OS are exponential. The
medians of PFS and OS are 4 and 12 months in the control group and 6 and
16 months in the treatment group, and 450 patients are enrolled at 25
per month.

PFS is tested at 300 PFS events, where OS is also tested, and OS is
tested again at 220 OS events. Each look has its own trigger, and `test`
lists the endpoints tested at each look:

``` r

fl <- ms_scenario(ms_fleischer(median_pfs = 4, median_os = 12),
                  ms_fleischer(median_pfs = 6, median_os = 16),
                  accrual = ms_accrual(rate = 25))
lk2 <- ms_looks(c("pfs", "os"), fraction = c(1, 1),
                test = list(c("pfs", "os"), "os"))
ev2 <- c(pfs = 300, os = 220)
```

## Fixed sequence

With the fixed-sequence procedure, both hypotheses are tested at level
0.025 and OS is tested only after PFS has been rejected: at the analysis
at which PFS is rejected and at the later analyses (procedure (b) of
Glimm et al. 2010, `hierarchical = "overall"`). Each endpoint has its
own efficacy bounds of O’Brien-Fleming type.

``` r

fs <- ms_power(fl, n = 450, looks = lk2,
               test = ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence"),
               events = ev2)
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

With `hierarchical = "stagewise"`, OS is tested only at the analysis at
which PFS is rejected, where the trial stops (procedure (a) of Glimm et
al. 2010; Tamhane et al. 2010). Here PFS is tested at the first analysis
only, so OS can be rejected only there:

``` r

fs_a <- ms_power(fl, n = 450, looks = lk2,
                 test = ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence",
                                hierarchical = "stagewise"),
                 events = ev2)
rbind(overall = fs$power, stagewise = fs_a$power)
#>                 pfs        os      both       any
#> overall   0.9387378 0.5451639 0.5451639 0.9387378
#> stagewise 0.9387378 0.2965913 0.2965913 0.9387378
```

## Bonferroni and recycling of the level

The Bonferroni procedure with the weights 0.8 and 0.2 tests PFS at level
0.02 and OS at level 0.005, whatever happens to the other hypothesis.
With `transition` in
[`ms_test()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md)
the level of a rejected hypothesis is passed on to the other in the
graphical procedure of Maurer and Bretz (2013). With
`transition = c(1, 1)`, OS is tested at level 0.025 once PFS has been
rejected, and PFS at level 0.025 once OS has been rejected.

``` r

bf <- ms_power(fl, n = 450, looks = lk2,
               test = ms_test(endpoints = c("pfs", "os"), mtp = "bonferroni",
                              weights = c(0.8, 0.2)),
               events = ev2)
rc <- ms_power(fl, n = 450, looks = lk2,
               test = ms_test(endpoints = c("pfs", "os"), mtp = "bonferroni",
                              weights = c(0.8, 0.2), transition = c(1, 1)),
               events = ev2)
rbind(fixed_sequence = fs$power, bonferroni = bf$power, recycling = rc$power)
#>                      pfs        os      both       any
#> fixed_sequence 0.9387378 0.5451639 0.5451639 0.9387378
#> bonferroni     0.9265375 0.3237669 0.3173922 0.9329122
#> recycling      0.9266956 0.5479409 0.5417243 0.9329122
rc$looks[, c("endpoint", "look", "time", "z_upper", "z_upper_recycled",
             "p_upper", "p_upper_recycled")]
#>   endpoint look     time  z_upper z_upper_recycled     p_upper p_upper_recycled
#> 1      pfs    1 18.62269 2.053749         1.959964 0.020000000      0.025000000
#> 2       os    1 18.62269 3.059549         2.357185 0.001108353      0.009207032
#> 3       os    2 23.05363 2.600794         2.009474 0.004650414      0.022243431
```

The power for OS is 0.324 without and 0.548 with recycling. The columns
`z_upper_recycled` and `p_upper_recycled` give the bounds at the level
after the other hypothesis has been rejected. The procedure controls the
familywise error rate when the nominal levels of the spending functions
do not decrease with the level (Maurer and Bretz 2013, Section 3.1), and
[`ms_power()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)
checks at every look that the bounds after recycling do not exceed the
bounds at the initial level.

## Familywise error rate

The familywise error rate is computed under the global null hypothesis,
where the correlation between the PFS and OS statistics is that of
Danzer et al. (2025, Theorem 1). A partial null hypothesis, in which OS
is the same in both groups but PFS is not, cannot be derived from the
design scenario and is given through `null_scenarios`:

``` r

os_null <- ms_scenario(ms_fleischer(median_pfs = 4, median_os = 12),
                       ms_fleischer(median_pfs = 6, median_os = 12),
                       accrual = ms_accrual(rate = 25))
fw <- ms_power(fl, n = 450, looks = lk2,
               test = ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence"),
               events = ev2, null_scenarios = list(os = os_null))
fw$fwer
#>              scenario null_endpoints nonbinding    binding
#> 1              global        PFS, OS 0.02500000 0.02500000
#> 2 null_scenarios[[1]]             OS 0.02405394 0.02405394
```

## An interim analysis of PFS

The following design adds an interim analysis of PFS at 60 percent of
the PFS events, as in the example of
[`ms_looks()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_looks.md):

``` r

lk3 <- ms_looks(c("pfs", "pfs", "os"), fraction = c(0.6, 1, 1),
                test = list("pfs", c("pfs", "os"), "os"))
fs3 <- ms_power(fl, n = 450, looks = lk3,
                test = ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence"),
                events = ev2)
fs3$looks[, c("endpoint", "look", "trigger", "time", "events", "z_upper",
              "prob_upper")]
#>   endpoint look trigger     time   events  z_upper prob_upper
#> 1      pfs    1     pfs 13.17448 180.0000 2.668630  0.5180390
#> 2      pfs    2     pfs 18.62269 300.0000 1.980965  0.4190402
#> 3       os    2     pfs 18.62269 162.9636 2.357185  0.2964566
#> 4       os    3      os 23.05363 220.0000 2.009474  0.2482003
rbind(two_analyses = fs$power, three_analyses = fs3$power)
#>                      pfs        os      both       any
#> two_analyses   0.9387378 0.5451639 0.5451639 0.9387378
#> three_analyses 0.9370793 0.5446569 0.5446569 0.9370793
```

The interim analysis of PFS takes place at month 13.2. Futility bounds,
when given, are placed on the first endpoint of the test.

## Correlation between the statistics

[`lr_cov()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md)
gives the joint distribution of the statistics. At month 24 the
correlation between the PFS and OS statistics is

``` r

cv <- lr_cov(fl, n = 450, calendar_time = 24, endpoint = c("pfs", "os"),
             variance = "studentized")
cv0 <- lr_cov(ms_null(fl), n = 450, calendar_time = 24, endpoint = c("pfs", "os"))
c(design = stats::cov2cor(cv$sigma)[1, 2], global_null = stats::cov2cor(cv0$sigma)[1, 2])
#>      design global_null 
#>   0.4518932   0.4204810
```

[`ms_size()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md)
finds the sample size, number of events, accrual duration or follow-up
for a target power of the first endpoint, of each endpoint, of both or
of at least one (argument `target`).

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

Maurer W, Bretz F (2013). Multiple testing in group sequential trials
using graphical approaches. Statistics in Biopharmaceutical Research
5:311-320.

Tamhane AC, Mehta CR, Liu L (2010). Testing a primary and a secondary
endpoint in a group sequential design. Biometrics 66:1174-1184.
