# State occupation probabilities

Probabilities of being in each state of the four-state model at given
times since randomization, without dropout.

## Usage

``` r
ms_prob(model, time)
```

## Arguments

- model:

  An `ms_model` object.

- time:

  Times since randomization (nonnegative).

## Value

A data frame with `time` and the probabilities `progression_free`,
`progressed` (without switch), `switched` and `dead`.

## Details

The probability of being progression-free and alive is `exp(-H(t))`,
where H is the cumulative hazard of progression and death. The
probabilities of the progressed and switched states are integrals over
the time of progression of the progression density times the transition
probabilities after progression, which have closed forms for piecewise
constant intensities. The integrals are computed with Gauss-Legendre
rules.

## See also

[`ms_model`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_model.md),
[`ms_surv`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_surv.md)

## Examples

``` r
m <- ms_model(time = c(0, 6), prog = c(0.10, 0.08), death = 0.01,
              death_prog = 0.05, switch_prob = 0.3, death_switch = 0.03)
ms_prob(m, time = c(6, 12, 24))
#>   time progression_free progressed  switched      dead
#> 1    6        0.5168513  0.2612947 0.1194071 0.1024469
#> 2   12        0.3011942  0.3079500 0.1519437 0.2389121
#> 3   24        0.1022842  0.2572271 0.1491483 0.4913404
```
