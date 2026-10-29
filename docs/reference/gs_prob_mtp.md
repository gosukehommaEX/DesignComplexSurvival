# Rejection probabilities of a group sequential design with two endpoints

For the standardized statistics of two endpoints with a joint
multivariate normal distribution, computes under a multiple testing
procedure the probability of rejecting each null hypothesis at each
look, the probability of rejecting both and the probability of rejecting
at least one.

## Usage

``` r
gs_prob_mtp(
  upper,
  lower = NULL,
  mean = 0,
  sigma,
  analysis = NULL,
  mtp = c("fixed_sequence", "bonferroni"),
  hierarchical = c("overall", "stagewise"),
  upper_recycled = NULL
)
```

## Arguments

- upper:

  List of two numeric vectors: the efficacy bounds of the first and of
  the second endpoint at their looks (`Inf` for no efficacy test at a
  look).

- lower:

  Futility bounds of the first endpoint at its looks (`-Inf` for none),
  or `NULL` for no futility bounds.

- mean:

  Mean vector of the statistics, first the looks of the first endpoint
  and then those of the second (recycled if of length 1).

- sigma:

  Covariance matrix of the statistics in the same order, for example
  `sigma` of
  [`lr_cov`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md)
  with two endpoints.

- analysis:

  List of two integer vectors: the analysis at which each look of the
  first and of the second endpoint takes place (strictly increasing
  within an endpoint; equal numbers mean the same analysis). `NULL`
  means that both endpoints are tested at analyses 1 to K.

- mtp:

  `"fixed_sequence"` or `"bonferroni"`.

- hierarchical:

  `"overall"` or `"stagewise"` (used with `mtp = "fixed_sequence"`).

- upper_recycled:

  For `mtp = "bonferroni"`, an optional list of two numeric vectors: the
  efficacy bounds of the first and of the second endpoint at their looks
  at the level after the other hypothesis has been rejected (not above
  `upper`). `NULL` means no recycling.

## Value

A list with

- `looks`:

  A data frame with one row per look: `endpoint` (1 or 2), `look`,
  `analysis`, the probability of rejecting the hypothesis of the
  endpoint at the look (`prob_upper`), its cumulative sum (`cum_upper`)
  and the probability of stopping for futility at the look
  (`prob_lower`, 0 for the second endpoint and for the last look of the
  first).

- `reject`:

  The probabilities of rejecting the first hypothesis (`first`), the
  second (`second`), both (`both`) and at least one (`any`).

## Details

Each endpoint is tested at its own looks, and `analysis` places the
looks of both endpoints on one sequence of analyses; an analysis may
include looks of one endpoint or of both. A hypothesis is rejected at
the first look at which its statistic is at or above its efficacy bound
and is not tested again. The trial stops for futility at an interim look
of the first endpoint at which its hypothesis has not been rejected and
its statistic is below the futility bound. The futility bound at the
last look of the first endpoint is not used.

With `mtp = "bonferroni"`, each endpoint is tested at all its looks
until the trial stops, with bounds computed at levels that add up to
alpha. The trial stops when both hypotheses are rejected or for
futility.

With `mtp = "bonferroni"` and `upper_recycled`, the level of a rejected
hypothesis is passed on to the other hypothesis (the graphical procedure
of Maurer and Bretz 2013, Algorithm 1, for two hypotheses). Each
hypothesis is tested with its bounds in `upper` until the other
hypothesis is rejected, and with its bounds in `upper_recycled` from the
analysis at which the other is rejected on: a hypothesis rejected with
its bound in `upper` at an analysis makes the other be tested with its
bound in `upper_recycled` at the same analysis. The bounds in
`upper_recycled` are those of the level after recycling, for example
`(w1 + g21 w2) alpha` for the first hypothesis with the weights w and
the fraction g21 of the level of the second hypothesis that is passed on
to the first. They must not exceed the bounds in `upper`, which holds
when the nominal levels of the spending function do not decrease with
the level (Maurer and Bretz 2013, Section 3.1). The futility stop of the
first endpoint also applies after the second hypothesis has been
rejected.

With `mtp = "fixed_sequence"`, the second hypothesis is tested only
after the first has been rejected. With `hierarchical = "overall"` it is
tested at the analysis at which the first hypothesis is rejected and at
the later looks of the second endpoint, until it is rejected (strategy
(b) of Glimm et al. 2010). With `hierarchical = "stagewise"` the trial
stops at the analysis at which the first hypothesis is rejected, and the
second hypothesis is tested only at that analysis (strategy (a) of Glimm
et al. 2010, the procedure of Tamhane et al. 2010).

The probabilities are sums of rectangle probabilities of the
multivariate normal distribution, one for each combination of the
outcome of the first endpoint (rejection or stop for futility at a look,
or no rejection) and the look at which the second hypothesis is
rejected. With `upper_recycled`, the probabilities are summed over the
outcomes at each analysis: at an analysis at which both hypotheses are
tested and neither has been rejected, with the bounds A1, A2 in `upper`
and B1, B2 in `upper_recycled`, both hypotheses are rejected when Z1 \>=
A1 and Z2 \>= B2 or when B1 \<= Z1 \< A1 and Z2 \>= A2, only the first
when Z1 \>= A1 and Z2 \< B2, only the second when Z1 \< B1 and Z2 \>=
A2, and neither when Z1 \< A1 and Z2 \< A2; at an analysis at which one
hypothesis is tested, it is rejected or not. The futility stop splits
the outcomes in which the first hypothesis is tested and not rejected at
an interim look. The rectangle probabilities are computed with
[`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html) (at
most 20 looks in total): with the Genz-Bretz algorithm for the
rectangles that involve both endpoints and have four or more coordinates
and, as in
[`gs_prob`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob.md),
for the other rectangles of four or more coordinates whose covariance
matrix is not that of a Markov sequence and of nine or more coordinates
with two finite limits for some coordinate, for which the Miwa algorithm
can be inaccurate, and with the Miwa algorithm otherwise. The Genz-Bretz
algorithm uses a fixed seed, so that the results are reproducible and
the random number generator of the session is not changed, and an
absolute error tolerance of 1e-5 for each rectangle, which can be
changed with `options(DesignComplexSurvival.genz_abseps = ...)`.

## References

Glimm E, Maurer W, Bretz F (2010). Hierarchical testing of multiple
endpoints in group-sequential trials. Statistics in Medicine 29:219-228.

Maurer W, Bretz F (2013). Multiple testing in group sequential trials
using graphical approaches. Statistics in Biopharmaceutical Research
5:311-320.

Tamhane AC, Mehta CR, Liu L (2010). Testing a primary and a secondary
endpoint in a group sequential design. Biometrics 66:1174-1184.

## See also

[`gs_prob`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_prob.md),
[`gs_bounds`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/gs_bounds.md),
[`lr_cov`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md)

## Examples

``` r
# Two looks at information fractions 0.5 and 1 for both endpoints and
# correlation 0.5 between the endpoints
tt <- c(0.5, 1)
r <- sqrt(outer(tt, tt, pmin) / outer(tt, tt, pmax))
sigma <- rbind(cbind(r, 0.5 * r), cbind(0.5 * r, r))
b1 <- gs_bounds(tt, sf_eff = "OF")$upper
b2 <- gs_bounds(tt, sf_eff = "Pocock")$upper
gs_prob_mtp(list(b1, b2), mean = c(2, 2.9, 1.5, 2.1), sigma = sigma)$reject
#>     first    second      both       any 
#> 0.8248839 0.4412122 0.4412122 0.8248839 
# Probability of rejecting the second hypothesis when it is true
gs_prob_mtp(list(b1, b2), mean = c(2, 2.9, 0, 0), sigma = sigma)$reject
#>      first     second       both        any 
#> 0.82488393 0.02068571 0.02068571 0.82488393 
# Bonferroni with the level alpha / 2 for each endpoint and the level
# alpha after the other hypothesis has been rejected (full recycling)
h1 <- gs_bounds(tt, alpha = 0.0125, sf_eff = "OF")$upper
h2 <- gs_bounds(tt, alpha = 0.0125, sf_eff = "Pocock")$upper
gs_prob_mtp(list(h1, h2), mean = c(2, 2.9, 1.5, 2.1), sigma = sigma,
            mtp = "bonferroni", upper_recycled = list(b1, b2))$reject
#>     first    second      both       any 
#> 0.7620396 0.4654847 0.4417205 0.7858038 
```
