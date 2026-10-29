# Two-group trial scenario

Combines the models of the control and treatment groups, the allocation
ratio and the accrual.

## Usage

``` r
ms_scenario(
  control,
  treatment,
  allocation = 1,
  accrual,
  pfs_switch = c("event", "censor"),
  switch_start = 0
)
```

## Arguments

- control, treatment:

  `ms_model` objects of the control and treatment groups.

- allocation:

  Allocation ratio, treatment to control (for example 2 for 2:1).

- accrual:

  An `ms_accrual` object.

- pfs_switch:

  Definition of PFS at a switch before progression: `"event"` (followed
  up to progression or death) or `"censor"` (censored at the switch). It
  applies to both groups.

- switch_start:

  Calendar time (time since the start of accrual) from which switching
  after progression is possible, nonnegative. The default 0 means that
  switching is possible from randomization, as in the models.

## Value

An object of class `ms_scenario`, a list with the two models,
`allocation`, the allocation probabilities `prob` (control, treatment),
`accrual`, `pfs_switch` and `switch_start`. The models carry
`pfs_switch` as well.

## Details

`pfs_switch` defines PFS for patients who switch before progression
(`switch_pre` in
[`ms_model`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md)).
With `"event"`, the switch does not change the definition: PFS is the
time to progression or death. With `"censor"`, PFS is censored at the
switch (for example when the start of a new treatment censors PFS); the
log-rank test of PFS then compares the hazard `prog + death` of the
patients who are progression-free and have not switched.

`switch_start` is a calendar time (time since the start of accrual)
before which nobody switches, for example when switching is allowed only
after an interim analysis. A patient enrolled at time e can then switch
from the follow-up time `max(switch_start - e, 0)` on: progressions
before that time do not lead to a switch at progression, and patients
who progressed before it can switch with `switch_rate` from that time
on. The distribution of OS then depends on the entry time, and the
computations for OS integrate over the entry time. They cover the OS
log-rank statistic alone with `mean = "direct"`; two endpoints together,
an OS test at looks driven by PFS events, and a PFS test at looks driven
by OS events with a variance other than `"null"` give an error. PFS is
not affected by the switching. Under the null hypothesis the score no
longer has independent increments when the switching changes the death
hazard: the correlations between the looks are computed from the
per-patient influence functions instead of the information (see
[`lr_cov`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md)).
`switch_start` cannot be combined with switching before progression
(`switch_pre`).

## See also

[`ms_model`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md),
[`ms_accrual`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_accrual.md),
[`ms_null`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_null.md)

## Examples

``` r
ctl <- ms_model(death = log(2) / 12)
trt <- ms_model(death = 0.7 * log(2) / 12)
ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#> Two-group scenario, allocation 1:1 (treatment:control)
#> 
#> Control group. Four-state model for one group (clock after progression: "reset")
#> Before progression (time since randomization):
#>  time prog      death switch_prob dropout
#>     0    0 0.05776227           0       0
#> 
#> Treatment group. Four-state model for one group (clock after progression: "reset")
#> Before progression (time since randomization):
#>  time prog      death switch_prob dropout
#>     0    0 0.04043359           0       0
#> 
#> Piecewise uniform accrual by rates (patients per unit time):
#>  time rate
#>     0   20

# Switching after progression in the control group, possible only from
# month 15 (for example after an interim analysis)
ctl_sw <- ms_model(prog = 0.1, death = 0.01, death_prog = 0.05,
                   switch_prob = 0.3, switch_rate = 0.04, death_switch = 0.03)
trt_sw <- ms_model(prog = 0.06, death = 0.008, death_prog = 0.045)
ms_scenario(ctl_sw, trt_sw, accrual = ms_accrual(rate = 20), switch_start = 15)
#> Two-group scenario, allocation 1:1 (treatment:control)
#> 
#> Control group. Four-state model for one group (clock after progression: "reset")
#> Before progression (time since randomization):
#>  time prog death switch_prob dropout
#>     0  0.1  0.01         0.3       0
#> After progression (time since progression):
#>  time_post death_prog switch_rate death_switch
#>          0       0.05        0.04         0.03
#> 
#> Treatment group. Four-state model for one group (clock after progression: "reset")
#> Before progression (time since randomization):
#>  time prog death switch_prob dropout
#>     0 0.06 0.008           0       0
#> After progression (time since progression):
#>  time_post death_prog switch_rate death_switch
#>          0      0.045           0        0.045
#> 
#> Switching possible from calendar time 15 (time since the start of accrual)
#> 
#> Piecewise uniform accrual by rates (patients per unit time):
#>  time rate
#>     0   20
```
