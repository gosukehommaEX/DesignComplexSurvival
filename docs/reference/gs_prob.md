# Stopping probabilities of a group sequential design

For statistics with a multivariate normal distribution, computes at each
look the probability of continuing to the look and then crossing the
efficacy bound (`prob_upper`) or falling below the futility bound
(`prob_lower`).

## Usage

``` r
gs_prob(lower, upper, mean = 0, sigma)
```

## Arguments

- lower:

  Futility bounds (`-Inf` for none). At the last look the bound is
  usually equal to the efficacy bound.

- upper:

  Efficacy bounds (`Inf` for none).

- mean:

  Mean vector (recycled if of length 1).

- sigma:

  Covariance matrix.

## Value

A data frame with `look`, `prob_upper`, `prob_lower` and their
cumulative sums `cum_upper` and `cum_lower`.

## Details

Rectangle probabilities are computed with
[`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html).
The Miwa algorithm, which is deterministic, is used for rectangles of at
most three coordinates and when the covariance matrix is that of a
Markov sequence (its inverse is tridiagonal), as for statistics with
independent increments. For four or more coordinates with other
covariance matrices, such as those of `variance = "alternative"` and
`"studentized"` in
[`ms_test`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_test.md),
the Miwa algorithm was found to be inaccurate, and for nine or more
coordinates with two finite limits for some coordinate (futility bounds
at eight or more earlier looks), for which mvtnorm leaves out terms of
the inclusion-exclusion of the Miwa algorithm, the Genz-Bretz algorithm
is used, with a fixed seed, so that the results are reproducible and the
random number generator of the session is not changed, and with an
absolute error tolerance of 1e-5 for each rectangle, which can be
changed with `options(DesignComplexSurvival.genz_abseps = ...)`.

## See also

[`gs_bounds`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_bounds.md)

## Examples

``` r
b <- gs_bounds(c(0.5, 1))
corr <- matrix(c(1, sqrt(0.5), sqrt(0.5), 1), 2)
gs_prob(b$lower, b$upper, mean = 0, sigma = corr)
#>   look  prob_upper prob_lower   cum_upper cum_lower
#> 1    1 0.001525323      0.000 0.001525323     0.000
#> 2    2 0.023474677      0.975 0.025000000     0.975
gs_prob(b$lower, b$upper, mean = c(2, 2.9), sigma = corr)
#>   look prob_upper prob_lower cum_upper cum_lower
#> 1    1  0.1678772  0.0000000 0.1678772 0.0000000
#> 2    2  0.6570068  0.1751161 0.8248839 0.1751161
```
