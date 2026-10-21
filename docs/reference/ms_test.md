# Test and stopping rules

Defines the one-sided log-rank test, the efficacy and futility bounds
and the approximation of the distribution of the statistics.

## Usage

``` r
ms_test(
  alpha = 0.025,
  endpoints = "os",
  sf_eff = "OF",
  param_eff = NULL,
  futility = c("none", "fixed", "beta"),
  fut_bound = NULL,
  fut_scale = c("z", "hr"),
  sf_fut = "OF",
  param_fut = NULL,
  beta = 0.1,
  binding = FALSE,
  mean = c("direct", "schoenfeld"),
  variance = c("null", "alternative", "studentized"),
  spending_time = c("events", "information"),
  mtp = c("none", "bonferroni", "fixed_sequence"),
  weights = NULL,
  hierarchical = c("overall", "stagewise")
)
```

## Arguments

- alpha:

  One-sided type I error (with two endpoints, the familywise error
  rate).

- endpoints:

  Tested endpoint, `"os"` or `"pfs"`, or both in the order of testing
  (for example `c("pfs", "os")`). With two endpoints the futility bounds
  are placed on the first, and with `mtp = "fixed_sequence"` the first
  is tested first.

- sf_eff, param_eff:

  Alpha spending function and its parameter (see
  [`gs_spend`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_spend.md)).
  For `"user"`, `param_eff` gives the cumulative proportions of `alpha`
  at the looks. With two endpoints, `sf_eff` has one value or one value
  per endpoint, and `param_eff` is a list with one element per endpoint,
  or a number or one number per endpoint.

- futility:

  `"none"`, `"fixed"` or `"beta"` (for the first endpoint).

- fut_bound:

  Fixed futility bounds at the interim looks (for `futility = "fixed"`).

- fut_scale:

  Scale of `fut_bound`: `"z"` or `"hr"`.

- sf_fut, param_fut:

  Beta spending function and its parameter.

- beta:

  Total type II error for beta spending in
  [`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md);
  [`ms_size`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md)
  uses one minus the target power instead.

- binding:

  Whether the futility bounds are binding.

- mean:

  Approximation of the mean of the statistic: `"direct"` or
  `"schoenfeld"` (see
  [`lr_moments`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_moments.md)).

- variance:

  Approximation of the covariance of the statistics: `"null"`,
  `"alternative"` or `"studentized"` (see
  [`lr_cov`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md)).

- spending_time:

  `"events"` or `"information"`.

- mtp:

  Multiple testing procedure: `"none"` (one endpoint), `"bonferroni"` or
  `"fixed_sequence"` (two endpoints).

- weights:

  Bonferroni weights of the two endpoints (positive, adding up to 1;
  defaults to equal weights).

- hierarchical:

  `"overall"` or `"stagewise"` (for `mtp = "fixed_sequence"`).

## Value

An object of class `ms_test`.

## Details

The efficacy bounds are obtained by alpha spending under the null
hypothesis. The correlation of the statistics under the null hypothesis
is taken as `sqrt(d_j / d_k)`, where d are the expected numbers of
events of the tested endpoint at the looks; it is exact under the null
hypothesis when the looks are driven by the events of the tested
endpoint. The spending time is the event fraction
(`spending_time = "events"`) or the information fraction under the
design scenario (`"information"`).

Futility bounds may be fixed, on the Z scale (`fut_scale = "z"`) or on
the hazard ratio scale (`fut_scale = "hr"`, converted with
`z = -log(hr) sqrt(p0 p1 d)`), or obtained by beta spending under the
design scenario. With `binding = FALSE` the efficacy bounds ignore the
futility bounds.

With two endpoints, each endpoint has its own efficacy bounds, computed
as above from the looks at which it is tested (see
[`ms_looks`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_looks.md)).
With `mtp = "bonferroni"`, the endpoint k is tested at level
`weights[k] * alpha`; the levels are not recycled after a rejection.
With `mtp = "fixed_sequence"`, both endpoints are tested at level
`alpha` and the second endpoint is tested only after the first has been
rejected: with `hierarchical = "overall"` at the analysis at which the
first is rejected and at its later looks until it is rejected (procedure
(b) of Glimm et al. 2010), and with `"stagewise"` only at the analysis
at which the first is rejected, where the trial stops (procedure (a) of
Glimm et al. 2010, Tamhane et al. 2010). The bounds of the second
endpoint do not depend on the correlation between the endpoints; with
them the familywise error rate is at most `alpha` (Glimm et al. 2010,
Tamhane et al. 2010). Futility bounds are placed on the first endpoint
only: the trial stops when its statistic falls below the futility bound
at an interim look at which it has not been rejected. See
[`gs_prob_mtp`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob_mtp.md)
for the stopping rules.

The efficacy bounds and the non-binding type I error (with two
endpoints, the familywise error rate under the global null hypothesis)
do not depend on `variance`; the power, the futility bounds by beta
spending and the familywise error rate under `null_scenarios` of
[`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)
do. In the simulations of `inst/validation/` with PFS and OS tested
together, the powers with `"null"` and `"studentized"` agreed with the
simulated powers within the Monte Carlo error, `"alternative"` gave
slightly larger powers, and `"studentized"` agreed best with the
simulated standard deviations and correlations of the statistics.

## References

Glimm E, Maurer W, Bretz F (2010). Hierarchical testing of multiple
endpoints in group-sequential trials. Statistics in Medicine 29:219-228.

Tamhane AC, Mehta CR, Liu L (2010). Testing a primary and a secondary
endpoint in a group sequential design. Biometrics 66:1174-1184.

## See also

[`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md),
[`gs_bounds`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_bounds.md),
[`gs_prob_mtp`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob_mtp.md)

## Examples

``` r
ms_test()
#> One-sided log-rank test of OS at level 0.025
#>   Efficacy: alpha spending "OF"
#>   Futility: none
#>   Mean: "direct"; variance: "null"; spending time: "events"
ms_test(futility = "beta", sf_fut = "HSD", param_fut = -2)
#> One-sided log-rank test of OS at level 0.025
#>   Efficacy: alpha spending "OF"
#>   Futility: beta spending "HSD" with beta = 0.1 (non-binding)
#>   Mean: "direct"; variance: "null"; spending time: "events"
ms_test(futility = "fixed", fut_bound = 1, fut_scale = "hr")
#> One-sided log-rank test of OS at level 0.025
#>   Efficacy: alpha spending "OF"
#>   Futility: fixed bounds on the hazard ratio scale (1) (non-binding)
#>   Mean: "direct"; variance: "null"; spending time: "events"
# PFS first, then OS, with Pocock-type spending for OS
ms_test(endpoints = c("pfs", "os"), sf_eff = c("OF", "Pocock"),
        mtp = "fixed_sequence")
#> One-sided log-rank tests of PFS and OS at familywise level 0.025
#>   Fixed sequence: PFS, then OS (overall)
#>   Efficacy: alpha spending "OF" (PFS), "Pocock" (OS)
#>   Futility (PFS): none
#>   Mean: "direct"; variance: "null"; spending time: "events"
ms_test(endpoints = c("pfs", "os"), mtp = "bonferroni", weights = c(0.8, 0.2))
#> One-sided log-rank tests of PFS and OS at familywise level 0.025
#>   Bonferroni: PFS at level 0.02, OS at level 0.005
#>   Efficacy: alpha spending "OF" (PFS), "OF" (OS)
#>   Futility (PFS): none
#>   Mean: "direct"; variance: "null"; spending time: "events"
```
