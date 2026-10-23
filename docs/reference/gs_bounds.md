# Bounds of a group sequential design

Computes one-sided efficacy bounds by alpha spending and, optionally,
futility bounds that are fixed or obtained by beta spending, for a
vector of standardized statistics with a multivariate normal
distribution.

## Usage

``` r
gs_bounds(
  info_frac,
  alpha = 0.025,
  sf_eff = "OF",
  param_eff = NULL,
  futility = c("none", "fixed", "beta"),
  fut_bound = NULL,
  sf_fut = "OF",
  param_fut = NULL,
  beta = 0.1,
  binding = FALSE,
  mean = NULL,
  corr = NULL,
  sigma_alt = NULL
)
```

## Arguments

- info_frac:

  Information fractions used for spending (strictly increasing and
  ending at 1).

- alpha:

  One-sided type I error.

- sf_eff, param_eff:

  Alpha spending function and its parameter (see
  [`gs_spend`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_spend.md)).

- futility:

  `"none"`, `"fixed"` or `"beta"`.

- fut_bound:

  Fixed futility bounds on the Z scale at looks 1 to K - 1 (for
  `futility = "fixed"`).

- sf_fut, param_fut:

  Beta spending function and its parameter.

- beta:

  Total type II error for beta spending.

- binding:

  Whether the futility bounds are binding.

- mean:

  Mean vector of the statistics under the alternative (required for
  `futility = "beta"`).

- corr:

  Correlation matrix of the statistics under the null hypothesis.
  Defaults to `sqrt(t_j / t_k)` for j \<= k with `t = info_frac`.

- sigma_alt:

  Covariance matrix of the statistics under the alternative. Defaults to
  `corr`.

## Value

A data frame with one row per look: `look`, `info_frac`, cumulative
`alpha_spent` and `beta_spent`, the efficacy bound `upper`, the futility
bound `lower` (`-Inf` when there is none) and the one-sided nominal
p-value of the efficacy bound `nominal_p`.

## Details

The efficacy bound at look k solves, under the null hypothesis (mean 0,
correlation `corr`), the equation that the probability of continuing to
look k and then crossing the bound equals the alpha spent at look k.
With `binding = TRUE` the continuation region at the earlier looks
includes the futility bounds; with `binding = FALSE` the futility bounds
are ignored for the efficacy bounds.

With `futility = "beta"`, the futility bound at look k solves, under the
alternative (mean `mean`, covariance `sigma_alt`), the equation that the
probability of continuing to look k and then falling below the bound
equals the beta spent at look k. The mean may take any value at each
look, so that non-proportional hazards are allowed. A futility bound
above the efficacy bound is set to the efficacy bound. At the last look
the futility bound equals the efficacy bound. With `binding = TRUE`, a
futility bound set to the efficacy bound at an interim look stops every
trial at that look, so that no alpha can be spent at the later looks;
the function then stops with an error that names that look.

The probabilities are computed as in
[`gs_prob`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob.md):
with the Miwa algorithm for rectangles of at most three coordinates and
for covariance matrices of Markov sequences, such as the default `corr`,
except for rectangles of nine or more coordinates with two finite limits
for some coordinate (nine or more looks with futility bounds at the
earlier looks), and with the Genz-Bretz algorithm otherwise (for example
for `sigma_alt` of `variance = "alternative"` or `"studentized"` with
four or more looks).

## References

Lan KKG, DeMets DL (1983). Discrete sequential boundaries for clinical
trials. Biometrika 70:659-663.

## See also

[`gs_spend`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_spend.md),
[`gs_prob`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob.md),
[`ms_power`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_power.md)

## Examples

``` r
# O'Brien-Fleming type efficacy bounds at 50, 75 and 100 percent
gs_bounds(c(0.5, 0.75, 1))
#>   look info_frac alpha_spent beta_spent    upper    lower   nominal_p
#> 1    1      0.50 0.001525323         NA 2.962588     -Inf 0.001525323
#> 2    2      0.75 0.009649325         NA 2.359018     -Inf 0.009161689
#> 3    3      1.00 0.025000000         NA 2.014084 2.014084 0.022000375

# Non-binding futility by beta spending with a delayed effect
gs_bounds(c(0.5, 0.75, 1), futility = "beta", sf_fut = "HSD",
          param_fut = -2, mean = c(0.8, 1.9, 2.8))
#>   look info_frac alpha_spent beta_spent    upper      lower   nominal_p
#> 1    1      0.50 0.001525323 0.02689414 2.962588 -1.1285377 0.001525323
#> 2    2      0.75 0.009649325 0.05449458 2.359018  0.1920799 0.009161689
#> 3    3      1.00 0.025000000 0.10000000 2.014084  2.0140837 0.022000375
```
