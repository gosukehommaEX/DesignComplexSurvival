# Multistate model for one treatment group

Defines the multistate model of one group: progression-free and alive
(state 1), progressed without switch (state 2), switched (state 3) and
dead (state 4), and, when `switch_pre` is positive, progression-free
after a switch before progression (state 5). Progression-free survival
(PFS) is the time to leave the progression-free states and overall
survival (OS) is the time to enter state 4.

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
  clock = c("reset", "forward"),
  switch_clock = c("same", "switch"),
  switch_af = NULL,
  switch_pre = 0,
  prog_switch = NULL,
  death_switch_pre = NULL,
  dropout_pfs = 0
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

- switch_clock:

  Clock of `death_switch`: `"same"` (the clock of the other
  post-progression intensities) or `"switch"` (time since the switch).

- switch_af:

  Acceleration factor of an accelerated failure time effect of the
  switch (positive), or `NULL` (default) for the hazard `death_switch`.
  When given, `death_switch` must not be given and `switch_clock` must
  be `"same"`.

- switch_pre:

  Hazard of switching before progression, in each interval of `time`.

- prog_switch:

  Hazard of progression after a switch before progression, in each
  interval of `time`. Defaults to `prog`.

- death_switch_pre:

  Hazard of death before progression after a switch, in each interval of
  `time`. Defaults to `death`.

- dropout_pfs:

  Hazard of a censoring of PFS only, in each interval of `time`.

## Value

An object of class `ms_model`: a list with the intensity vectors
expanded to the lengths of `time` and `time_post`, `clock`,
`switch_clock`, `switch_af` and `pfs_switch` (`"event"`; set by
[`ms_scenario`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md)).
With `switch_af`, the element `death_switch` holds
`switch_af * death_prog`, the death hazard after a switch when
`death_prog` is constant.

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

The death hazard after a switch at the post-progression clock value
\\\sigma\\ is, at the clock value \\\tau\\, `death_switch` at \\\tau\\
with `switch_clock = "same"` (default), and `death_switch` at \\\tau -
\sigma\\ (the time since the switch) with `switch_clock = "switch"`.
With `clock = "reset"` the two differ only for switches after
progression (`switch_rate`); with `clock = "forward"` they differ also
for switches at progression. With `switch_af` given, the effect of a
switch is an accelerated failure time effect: the death hazard after the
switch is \\\phi c(\sigma + \phi (\tau - \sigma))\\, where \\\phi\\ is
`switch_af` and \\c\\ is `death_prog`, so that the remaining time is
that without switch divided by \\\phi\\ (the rank preserving structural
failure time model of Robins and Tsiatis 1991; \\\phi \< 1\\ prolongs
survival).

A switch before progression moves a patient from state 1 to state 5 with
hazard `switch_pre`. In state 5 the patient progresses with hazard
`prog_switch` and dies with hazard `death_switch_pre`; at progression
the patient enters state 3 and then dies with the hazard `death_switch`
on the post-progression clock, whatever `switch_clock` (with
`switch_af`, the hazard `switch_af * death_prog`). Whether a switch
before progression censors PFS is set in
[`ms_scenario`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md)
(`pfs_switch`).

Dropout is piecewise exponential with intervals starting at `time`,
independent of the event process, and censors both endpoints. The hazard
`dropout_pfs` adds a censoring that acts on PFS only (for example the
end of tumor assessments), independent of everything else.

A model without progression (`prog = 0`) is a piecewise exponential
model for a single endpoint with hazard `death`. With
`death_prog = death`, `switch_prob = 0` and `switch_rate = 0`, the model
is that of Fleischer et al. (2009); see
[`ms_fleischer`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_fleischer.md).

## References

Fleischer F, Gaschler-Markefski B, Bluhmki E (2009). A statistical model
for the dependence between progression-free survival and overall
survival. Statistics in Medicine 28:2669-2686.

Robins JM, Tsiatis AA (1991). Correcting for non-compliance in
randomized trials using rank preserving structural failure time models.
Communications in Statistics - Theory and Methods 20:2609-2631.

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

# Switching after progression with an accelerated failure time effect
ms_model(prog = 0.1, death = 0.01, death_prog = 0.05, switch_rate = 0.04,
         switch_af = 0.7)
#> Four-state model for one group (clock after progression: "reset")
#> Before progression (time since randomization):
#>  time prog death switch_prob dropout
#>     0  0.1  0.01           0       0
#> After progression (time since progression):
#>  time_post death_prog switch_rate death_switch
#>          0       0.05        0.04        0.035
#> Death after a switch: accelerated failure time effect with factor 0.7 ('death_switch' is 0.7 * death_prog).

# Switching before progression, and censoring of PFS only
ms_model(prog = 0.1, death = 0.01, death_prog = 0.05, death_switch = 0.03,
         switch_pre = 0.02, dropout_pfs = 0.01)
#> Five-state model for one group (clock after progression: "reset")
#> Before progression (time since randomization):
#>  time prog death switch_prob dropout dropout_pfs
#>     0  0.1  0.01           0       0        0.01
#> Switch before progression (time since randomization):
#>  time switch_pre prog_switch death_switch_pre
#>     0       0.02         0.1             0.01
#> After progression (time since progression):
#>  time_post death_prog switch_rate death_switch
#>          0       0.05           0         0.03
```
