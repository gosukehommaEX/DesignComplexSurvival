# Number of events, sample size, accrual duration or follow-up for a target power

Solves for one design parameter so that the design of
[`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)
has the target power, with the other parameters fixed.

## Usage

``` r
ms_size(
  scenario,
  looks,
  test,
  power = 0.9,
  solve_for = c("events", "n", "accrual_duration", "followup"),
  n = NULL,
  followup = NULL,
  interval = NULL,
  round_up = TRUE,
  target = c("first", "each", "both", "any"),
  null_scenarios = NULL
)
```

## Arguments

- scenario:

  An `ms_scenario` object.

- looks:

  An `ms_looks` object.

- test:

  An `ms_test` object.

- power:

  Target power.

- solve_for:

  `"events"`, `"n"`, `"accrual_duration"` or `"followup"`.

- n:

  Total sample size (fixed for `"events"`, `"followup"` and
  `"accrual_duration"`).

- followup:

  Time from the end of accrual to the final analysis (fixed for `"n"`
  and `"accrual_duration"` with event-driven looks).

- interval:

  Optional search interval for the unknown (for `"events"` and
  `"followup"`, an interval of the calendar time of the final analysis).

- round_up:

  Whether to round the number of events or the sample size up to an
  integer.

- target:

  For two endpoints, `"first"`, `"each"`, `"both"` or `"any"` (see
  Details).

- null_scenarios:

  For two endpoints, null scenarios for the familywise error rate of the
  solved design (see
  [`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)).

## Value

An object of class `c("ms_size", "ms_power")`: the result of
[`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)
for the solved design, with `solve_for`, the unrounded `solution` (named
by endpoint for the final numbers of events of two trigger endpoints),
the `target_power` and, for two endpoints, the `target`.

## Details

The unknown is one of

- `"events"`: the final number of events of the trigger endpoint, with
  `n` and the accrual fixed (event-driven looks);

- `"followup"`: the time from the end of accrual to the final analysis,
  with `n` and the accrual fixed (event-driven looks);

- `"n"`: the sample size, with the accrual rates (or, for an accrual
  given by proportions, the accrual duration) fixed and the final
  analysis at `followup` after the end of accrual (event-driven looks)
  or at the calendar times of the looks;

- `"accrual_duration"`: the accrual duration of an accrual given by
  proportions, with `n` fixed and the final analysis at `followup` after
  the end of accrual (event-driven looks) or at the calendar times of
  the looks.

For `"events"` and `"followup"` the search is over the calendar time of
the final analysis, so that the two give the same design (before
rounding). The search doubles or halves the unknown until the power
brackets the target, evaluates the power on a grid within the bracket
(or within `interval`) and refines the first grid interval in which the
target is reached, so that the smallest solution is returned when the
power is not monotone. Beta spending uses one minus `power` as the total
type II error.

When the looks are driven by the events of two endpoints, the final
numbers of events of both are their expected numbers at the final
analysis (see
[`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)),
so that the looks driven by the endpoint other than that of the final
look need fractions below 1. With `"events"` the solution then has the
final numbers of events of both trigger endpoints, each rounded up with
`round_up = TRUE`.

With two endpoints, `target` chooses the quantity that reaches `power`:
the power of the first endpoint (`"first"`), the smaller of the powers
of the two endpoints (`"each"`), the probability of rejecting both
hypotheses (`"both"`) or that of rejecting at least one (`"any"`), all
under the multiple testing procedure of `test`.

## See also

[`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)

## Examples

``` r
ctl <- ms_model(death = log(2) / 12)
trt <- ms_model(death = 0.7 * log(2) / 12)
scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
ms_size(scn, looks = ms_looks("os", fraction = c(0.5, 1)),
        test = ms_test(), power = 0.9, solve_for = "events", n = 500)
#> Solved for the final number of events: 332.317 (target power 0.9)
#> 
#> Group sequential log-rank design
#>   Tested endpoint: OS; looks driven by OS events
#>   Sample size: 500; accrual duration: 25; final analysis at: 36.49
#>   Power: 0.9006
#>   Type I error: 0.025 (non-binding), 0.025 (binding)
#>  look  time events info_frac z_upper  p_upper hr_upper z_lower hr_lower
#>     1 21.75  166.5       0.5   2.963 0.001525   0.6318    -Inf       NA
#>     2 36.49  333.0       1.0   1.969 0.024500   0.8059   1.969   0.8059
#>  prob_upper prob_lower cum_power
#>      0.2516    0.00000    0.2516
#>      0.6489    0.09941    0.9006

# \donttest{
# Sample size for a power of 80 per cent for both PFS and OS, tested in
# this order at one analysis 12 months after the end of accrual
fl <- ms_scenario(ms_fleischer(4, 12), ms_fleischer(6, 16),
                  accrual = ms_accrual(rate = 25))
ms_size(fl, looks = ms_looks("os"),
        test = ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence"),
        power = 0.8, solve_for = "n", followup = 12, target = "each")
#> Solved for the sample size: 567.43 (target power of PFS and of OS 0.8)
#> 
#> Group sequential log-rank design with two endpoints
#>   Endpoints: PFS, then OS (fixed sequence, overall); looks driven by OS events
#>   Sample size: 568; accrual duration: 22.72; final analysis at: 34.72
#>   Power: PFS 0.9965, OS 0.8005; both 0.8005; at least one 0.9965
#>   FWER under global null: 0.025 (non-binding), 0.025 (binding)
#>  endpoint look  time events info_frac z_upper p_upper hr_upper z_lower hr_lower
#>       PFS    1 34.72  534.1         1    1.96   0.025   0.8440    1.96   0.8440
#>        OS    1 34.72  381.5         1    1.96   0.025   0.8182    1.96   0.8182
#>  prob_upper prob_lower cum_power
#>      0.9965          0    0.9965
#>      0.8005          0    0.8005
# }
```
