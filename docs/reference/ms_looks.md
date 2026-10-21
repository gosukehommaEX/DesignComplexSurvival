# Timing of the analyses

Defines when the interim and final analyses take place: when the
expected number of events of an endpoint reaches given fractions of its
final number of events, or at given calendar times, and which endpoints
are tested at each analysis.

## Usage

``` r
ms_looks(
  trigger = c("os", "pfs", "calendar"),
  fraction = NULL,
  time = NULL,
  test = NULL
)
```

## Arguments

- trigger:

  `"os"` or `"pfs"` for event-driven looks, or `"calendar"`; one value
  for all looks or one value per look.

- fraction:

  Event fractions of the event-driven looks: with one trigger, strictly
  increasing and ending at 1 (defaults to 1, a single analysis); with
  one trigger per look, one value per look in (0, 1\] (`NA` at
  calendar-driven looks).

- time:

  Strictly increasing calendar times of the calendar-driven looks (with
  one trigger per look, one value per look, `NA` at event-driven looks).

- test:

  `NULL` (all endpoints at all looks), or a list with one element per
  look giving the endpoints tested at the look (`"os"`, `"pfs"` or
  both). A character vector with one endpoint per look is also accepted.

## Value

An object of class `ms_looks`.

## Details

With one value of `trigger`, all looks are driven by the events of the
same endpoint (`"os"` or `"pfs"`, with `fraction`) or take place at
calendar times (`"calendar"`, with `time`). The final number of events
is given in
[`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)
through `events` or `followup`, and the look k takes place when the
expected number of events of the trigger endpoint reaches `fraction[k]`
times the final number.

With one value of `trigger` per look, each look has its own trigger.
`fraction` then has one value per look, `NA` at the calendar-driven
looks, and `time` one value per look, `NA` at the event-driven looks.
The fractions of the looks driven by the events of an endpoint refer to
the final number of events of that endpoint (see
[`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md))
and must be strictly increasing; an event-driven final look must have
the fraction 1. The calendar times of the looks, which depend on the
scenario, must be strictly increasing.

By default every look tests all endpoints of the test
([`ms_test`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md)).
`test` gives the endpoints tested at each look, for example
`list("pfs", c("pfs", "os"), "os")` for an interim analysis of PFS, the
final analysis of PFS with an interim analysis of OS, and the final
analysis of OS.

## See also

[`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md),
[`ms_size`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md)

## Examples

``` r
ms_looks("os", fraction = c(0.5, 0.75, 1))
#> 3 look(s) driven by OS events at fractions: 0.50, 0.75, 1.00
ms_looks("calendar", time = c(24, 36))
#> 2 look(s) at calendar times: 24, 36
# PFS interim and final analyses driven by PFS events, an OS interim
# analysis at the final PFS analysis and the final OS analysis driven by
# OS events
ms_looks(c("pfs", "pfs", "os"), fraction = c(0.6, 1, 1),
         test = list("pfs", c("pfs", "os"), "os"))
#> 3 look(s):
#>   1: driven by PFS events at fraction 0.6; tests PFS
#>   2: driven by PFS events at fraction 1; tests PFS and OS
#>   3: driven by OS events at fraction 1; tests OS
```
