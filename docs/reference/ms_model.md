# Multistate model for one treatment group

Defines the four-state model of one group: progression-free and alive
(state 1), progressed without switch (state 2), switched (state 3) and
dead (state 4). Progression-free survival (PFS) is the time to leave
state 1 and overall survival (OS) is the time to enter state 4.

## Usage

``` r
ms_model(
  time = 0,
  prog = 0,
  death,
  death_prog = NULL,
  switch_prob = 0,
  switch_rate = 0,
  death_switch = NULL,
  dropout = 0,
  time_post = 0,
  clock = c("reset", "forward")
)
```

## Arguments

- time:

  Interval start times for the pre-progression intensities and the
  dropout hazard. Must start at 0 and be strictly increasing.

- prog:

  Progression hazard in each interval (length 1 or `length(time)`).

- death:

  Hazard of death before progression in each interval.

- death_prog:

  Hazard of death after progression without switch, in each interval of
  `time_post`. Required when `prog` is positive in some interval.

- switch_prob:

  Probability of switching immediately at progression, in each interval
  of `time`.

- switch_rate:

  Hazard of switching after progression, in each interval of
  `time_post`.

- death_switch:

  Hazard of death after a switch, in each interval of `time_post`.
  Defaults to `death_prog`, that is, switching has no effect on
  survival.

- dropout:

  Dropout hazard in each interval of `time`.

- time_post:

  Interval start times for the post-progression intensities. Must start
  at 0 and be strictly increasing.

- clock:

  Clock of the post-progression intensities: `"reset"` (time since
  progression) or `"forward"` (time since randomization).

## Value

An object of class `ms_model`: a list with the intensity vectors
expanded to the lengths of `time` and `time_post`, and `clock`.

## Details

Before progression, the intensities are piecewise constant in the time
since randomization, with intervals starting at `time`: progression
(`prog`), death before progression (`death`), and the probability
`switch_prob` of switching immediately at progression (evaluated at the
time of progression). After progression without switch, the patient dies
with intensity `death_prog` or switches with intensity `switch_rate`;
after a switch, the patient dies with intensity `death_switch`. The
post-progression intensities are piecewise constant with intervals
starting at `time_post`.

With `clock = "reset"` (default) the post-progression intensities are
functions of the time since progression (a semi-Markov model). With
`clock = "forward"` they are functions of the time since randomization
(a Markov model). The two clocks give the same model when the
post-progression intensities are constant.

Dropout is piecewise exponential with intervals starting at `time`,
independent of the event process, and censors both endpoints.

A model without progression (`prog = 0`) is a piecewise exponential
model for a single endpoint with hazard `death`. With
`death_prog = death`, `switch_prob = 0` and `switch_rate = 0`, the model
is that of Fleischer et al. (2009); see
[`ms_fleischer`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_fleischer.md).

## References

Fleischer F, Gaschler-Markefski B, Bluhmki E (2009). A statistical model
for the dependence between progression-free survival and overall
survival. Statistics in Medicine 28:2669-2686.

## See also

[`ms_fleischer`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_fleischer.md),
[`ms_scenario`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md),
[`ms_surv`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_surv.md),
[`ms_prob`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_prob.md)

## Examples

``` r
# Single endpoint: exponential survival with a median of 12 months
ms_model(death = log(2) / 12)
#> Four-state model for one group (clock after progression: "reset")
#> Before progression (time since randomization):
#>  time prog      death switch_prob dropout
#>     0    0 0.05776227           0       0

# Progression, with 30 percent of the patients switching at progression
ms_model(time = c(0, 6), prog = c(0.10, 0.08), death = 0.01,
         death_prog = 0.05, switch_prob = 0.3, death_switch = 0.03)
#> Four-state model for one group (clock after progression: "reset")
#> Before progression (time since randomization):
#>  time prog death switch_prob dropout
#>     0 0.10  0.01         0.3       0
#>     6 0.08  0.01         0.3       0
#> After progression (time since progression):
#>  time_post death_prog switch_rate death_switch
#>          0       0.05           0         0.03
```
