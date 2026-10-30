# Multistate model and treatment switching

``` r

library(DesignComplexSurvival)
```

[`ms_model()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md)
defines the model of one treatment group. The patient moves between five
states with piecewise constant intensities:

1.  progression-free and alive,
2.  progressed without switch,
3.  switched,
4.  dead,
5.  progression-free after a switch before progression (used only when
    `switch_pre` is positive).

Progression-free survival (PFS) is the time to leave the
progression-free states 1 and 5, and overall survival (OS) is the time
to enter state 4. Before progression the intensities are functions of
the time since randomization; after progression they are functions of
the time since progression or since randomization (the clock, see
below). Dropout is independent of the event process. This article shows
how switching at, after and before progression enters the model, and how
the options change OS, PFS and the power of the log-rank test.

## Switching at and after progression

In the control group the hazard of progression is 0.10 per month in the
first 4 months and 0.08 afterwards. At progression, 30 percent of the
patients switch to the experimental treatment, and the others may switch
later with a hazard of 0.04 per month. After progression the hazard of
death is 0.06 without switch and 0.035 after a switch. The treatment
group has a lower hazard of progression and no switching.

``` r

ctl_sw <- ms_model(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01,
                   switch_prob = 0.3, death_prog = 0.06, switch_rate = 0.04,
                   death_switch = 0.035, dropout = 0.003)
trt_sw <- ms_model(time = c(0, 4), prog = c(0.06, 0.05), death = 0.008,
                   death_prog = 0.055, dropout = 0.003)
pr <- ms_prob(ctl_sw, time = c(6, 12, 24))
pr
#>   time progression_free progressed  switched      dead
#> 1    6        0.5379444  0.2102858 0.1417831 0.1099867
#> 2   12        0.3134862  0.2176410 0.2161661 0.2527068
#> 3   24        0.1064585  0.1329666 0.2539857 0.5065892
```

The state probabilities do not include dropout. At 24 months, 25 percent
of the control patients are alive after a switch and 51 percent have
died.

The control group without switching has the same intensities with
`switch_prob` and `switch_rate` set to zero. The OS curves of the three
models:

``` r

ctl_ns <- ms_model(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01,
                   death_prog = 0.06, dropout = 0.003)
models <- list(control = ctl_sw, control_no_switch = ctl_ns, treatment = trt_sw)
tt <- seq(0, 48, by = 1)
os <- sapply(models, function(m) ms_surv(m, tt, "os")$surv)
matplot(tt, os, type = "l", lty = 1:3, col = 1, ylim = c(0, 1),
        xlab = "Months since randomization", ylab = "Overall survival")
legend("topright", c("Control with switching", "Control without switching",
                     "Treatment"), lty = 1:3, bty = "n")
```

![Overall survival of the control group with and without switching and
of the treatment group over 48
months](switching_files/figure-html/os-curves-1.png)

``` r

median_os <- function(m) {
  stats::uniroot(function(t) ms_surv(m, t, "os")$surv - 0.5, c(0.1, 500))$root
}
med <- sapply(models, median_os)
round(med, 1)
#>           control control_no_switch         treatment 
#>              23.6              20.5              26.6
```

The median OS is 23.6 months in the control group with switching, 20.5
months without switching and 26.6 months in the treatment group. PFS
does not depend on the switching after progression.

The design has an interim analysis at 60 percent of the final 300 OS
events and 500 patients enrolled at 25 per month:

``` r

looks <- ms_looks("os", fraction = c(0.6, 1))
pw <- ms_power(ms_scenario(ctl_sw, trt_sw, accrual = ms_accrual(rate = 25)),
               n = 500, looks = looks, test = ms_test(), events = 300)
pw_ns <- ms_power(ms_scenario(ctl_ns, trt_sw, accrual = ms_accrual(rate = 25)),
                  n = 500, looks = looks, test = ms_test(), events = 300)
c(switching = pw$power, no_switching = pw_ns$power)
#>    switching no_switching 
#>    0.2205938    0.8426871
```

Switching to the experimental treatment in the control group moves the
hazard ratio for OS towards 1: the power is 0.221 with switching and
0.843 without.

## The clock of the post-progression intensities

The post-progression intensities are piecewise constant on intervals
that start at `time_post`. With `clock = "reset"` (the default) the
intervals refer to the time since progression (a semi-Markov model);
with `clock = "forward"` they refer to the time since randomization (a
Markov model). In the following model the hazard of death after
progression is 0.08 in the first 6 months and 0.04 afterwards:

``` r

m_reset <- ms_model(prog = 0.08, death = 0.01, time_post = c(0, 6),
                    death_prog = c(0.08, 0.04))
m_fwd <- ms_model(prog = 0.08, death = 0.01, time_post = c(0, 6),
                  death_prog = c(0.08, 0.04), clock = "forward")
tc <- c(12, 24, 36)
os_clock <- sapply(list(reset = m_reset, forward = m_fwd),
                   function(m) ms_surv(m, tc, "os")$surv)
rownames(os_clock) <- paste("month", tc)
round(os_clock, 3)
#>          reset forward
#> month 12 0.709   0.756
#> month 24 0.458   0.525
#> month 36 0.290   0.344
```

With the reset clock, a patient who progresses at month 10 has the
hazard 0.08 until month 16. With the forward clock the hazard is 0.04
from month 6 on, so this patient has the hazard 0.04 at once, and OS is
longer. PFS is the same under both clocks. When the post-progression
intensities are constant, the two clocks give the same model.

## The death hazard after a switch

By default (`switch_clock = "same"`) the hazard `death_switch` is
measured on the clock of the other post-progression intensities. With
`switch_clock = "switch"` it is measured from the time of the switch. In
the following model the patients switch after progression with a hazard
of 0.05 per month, and the death hazard after a switch is 0.03 in the
first interval and 0.05 afterwards:

``` r

m_same <- ms_model(prog = 0.08, death = 0.01, time_post = c(0, 6),
                   death_prog = c(0.08, 0.06), switch_rate = 0.05,
                   death_switch = c(0.03, 0.05))
m_since <- ms_model(prog = 0.08, death = 0.01, time_post = c(0, 6),
                    death_prog = c(0.08, 0.06), switch_rate = 0.05,
                    death_switch = c(0.03, 0.05), switch_clock = "switch")
os_switch <- sapply(list(same = m_same, since_switch = m_since),
                    function(m) ms_surv(m, tc, "os")$surv)
rownames(os_switch) <- paste("month", tc)
round(os_switch, 3)
#>           same since_switch
#> month 12 0.711        0.714
#> month 24 0.430        0.440
#> month 36 0.246        0.257
```

A patient who switches 4 months after progression has the hazard 0.03
for 2 months with `"same"` and for 6 months with `"switch"`. With
`clock = "reset"` the two choices differ only for switches after
progression (`switch_rate`); switches at progression (`switch_prob`)
start both clocks at the same time.

## An accelerated failure time effect of the switch

With `switch_af`, the switch has an accelerated failure time effect: the
remaining time after the switch is the remaining time without switch
divided by `switch_af`, as in the rank preserving structural failure
time model of Robins and Tsiatis (1991). A value below 1 prolongs
survival. When `death_prog` is constant, this is the same as the death
hazard `switch_af * death_prog` after the switch:

``` r

m_af1 <- ms_model(prog = 0.08, death = 0.01, death_prog = 0.06,
                  switch_rate = 0.05, switch_af = 0.7)
m_ph1 <- ms_model(prog = 0.08, death = 0.01, death_prog = 0.06,
                  switch_rate = 0.05, death_switch = 0.7 * 0.06)
d_af <- max(abs(ms_surv(m_af1, tt, "os")$surv - ms_surv(m_ph1, tt, "os")$surv))
d_af
#> [1] 5.551115e-17
```

The largest difference between the two OS curves over 48 months is
5.6e-17. When `death_prog` changes over time, the accelerated failure
time effect also changes the time scale on which the hazard after the
switch is read, and the two models differ:

``` r

m_af2 <- ms_model(prog = 0.08, death = 0.01, time_post = c(0, 6),
                  death_prog = c(0.08, 0.04), switch_rate = 0.05, switch_af = 0.7)
m_ph2 <- ms_model(prog = 0.08, death = 0.01, time_post = c(0, 6),
                  death_prog = c(0.08, 0.04), switch_rate = 0.05,
                  death_switch = 0.7 * c(0.08, 0.04))
os_af <- sapply(list(aft = m_af2, hazard_ratio = m_ph2),
                function(m) ms_surv(m, tc, "os")$surv)
rownames(os_af) <- paste("month", tc)
round(os_af, 3)
#>            aft hazard_ratio
#> month 12 0.716        0.718
#> month 24 0.479        0.483
#> month 36 0.325        0.328
```

## Switching before progression

A patient who is progression-free may switch with the hazard
`switch_pre` and enters state 5. There the hazards of progression and
death are `prog_switch` and `death_switch_pre`; at progression the
patient enters state 3 and dies with the hazard `death_switch`.

``` r

ctl_pre <- ms_model(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01,
                    death_prog = 0.06, death_switch = 0.035, switch_pre = 0.02,
                    prog_switch = 0.06, death_switch_pre = 0.007, dropout = 0.003)
ms_prob(ctl_pre, time = c(6, 12, 24))
#>   time progression_free progression_free_switched progressed   switched
#> 1    6       0.47711392                0.06723685  0.3236912 0.01369166
#> 2   12       0.24659696                0.07873972  0.3638703 0.03572292
#> 3   24       0.06587475                0.05593000  0.2637657 0.06368279
#>        dead
#> 1 0.1182664
#> 2 0.2750701
#> 3 0.5507468
```

The column `progression_free_switched` is the probability of state 5.
`pfs_switch` in
[`ms_scenario()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md)
sets whether PFS is followed up to progression or death after such a
switch (`"event"`, the default) or censored at the switch (`"censor"`,
for example when the start of a new treatment censors PFS). The
censoring appears among the PFS dropouts of
[`ms_events()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_events.md):

``` r

acc <- ms_accrual(rate = 25)
scn_event <- ms_scenario(ctl_pre, trt_sw, accrual = acc)
scn_censor <- ms_scenario(ctl_pre, trt_sw, accrual = acc, pfs_switch = "censor")
ev <- rbind(event = ms_events(scn_event, n = 500, calendar_time = 24, endpoint = "pfs"),
            censor = ms_events(scn_censor, n = 500, calendar_time = 24, endpoint = "pfs"))
ev[, c("events_control", "dropouts_control", "events_treatment", "dropouts_treatment")]
#>        events_control dropouts_control events_treatment dropouts_treatment
#> event        169.6267         5.294087         134.7651           6.520328
#> censor       155.9178        35.881942         134.7651           6.520328
pfs_power <- sapply(list(event = scn_event, censor = scn_censor), function(s)
  ms_power(s, n = 500, looks = ms_looks("pfs"), test = ms_test(endpoints = "pfs"),
           events = 300)$power)
pfs_power
#>     event    censor 
#> 0.9582631 0.9776734
```

With 300 PFS events at a single analysis, the power of the PFS test is
0.958 when the switch does not change the definition of PFS and 0.978
when it censors PFS.

## Censoring of PFS only

`dropout` censors both endpoints. `dropout_pfs` adds a censoring of PFS
only, for example the end of tumour assessments:

``` r

ctl_dp <- ms_model(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01,
                   death_prog = 0.06, dropout = 0.003, dropout_pfs = 0.02)
scn_dp <- ms_scenario(ctl_dp, trt_sw, accrual = acc)
rbind(pfs = ms_events(scn_dp, n = 500, calendar_time = 24, endpoint = "pfs"),
      os = ms_events(scn_dp, n = 500, calendar_time = 24, endpoint = "os"))[
  , c("events_control", "dropouts_control")]
#>     events_control dropouts_control
#> pfs      155.91779        35.881942
#> os        80.73701         8.383434
```

## References

Robins JM, Tsiatis AA (1991). Correcting for non-compliance in
randomized trials using rank preserving structural failure time models.
Communications in Statistics - Theory and Methods 20:2609-2631.
