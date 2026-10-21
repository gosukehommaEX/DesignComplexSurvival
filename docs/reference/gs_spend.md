# Error spending functions

Cumulative error spent at information fractions `t` for a total
(one-sided) error `alpha`.

## Usage

``` r
gs_spend(
  t,
  alpha = 0.025,
  type = c("OF", "Pocock", "HSD", "power", "user"),
  param = NULL
)
```

## Arguments

- t:

  Information fractions between 0 and 1.

- alpha:

  Total error to spend.

- type:

  Type of spending function: `"OF"`, `"Pocock"`, `"HSD"`, `"power"` or
  `"user"`.

- param:

  Parameter of the spending function (required for `"HSD"`, `"power"`
  and `"user"`).

## Value

Numeric vector of cumulative error spent.

## Details

The available types are

- `"OF"`: Lan-DeMets spending function of O'Brien-Fleming type,
  `2 - 2 pnorm(qnorm(1 - alpha / 2) / sqrt(t))`.

- `"Pocock"`: Lan-DeMets spending function of Pocock type,
  `alpha log(1 + (exp(1) - 1) t)`.

- `"HSD"`: Hwang-Shih-DeCani family with parameter `param` (gamma),
  `alpha (1 - exp(-gamma t)) / (1 - exp(-gamma))`, and `alpha t` for
  gamma = 0.

- `"power"`: Kim-DeMets power family with parameter `param` (rho),
  `alpha t^rho`.

- `"user"`: cumulative proportions of `alpha` given in `param`, one
  value per element of `t`.

The same functions are used for beta spending, with `alpha` replaced by
the total type II error.

## References

Lan KKG, DeMets DL (1983). Discrete sequential boundaries for clinical
trials. Biometrika 70:659-663.

## See also

[`gs_bounds`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_bounds.md)

## Examples

``` r
gs_spend(c(0.5, 0.75, 1), alpha = 0.025, type = "OF")
#> [1] 0.001525323 0.009649325 0.025000000
gs_spend(c(0.5, 0.75, 1), alpha = 0.1, type = "HSD", param = -2)
#> [1] 0.02689414 0.05449458 0.10000000
```
