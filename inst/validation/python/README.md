# Independent Python implementation

The scripts in this folder compute the expected values used in
`tests/testthat`. They are written independently of the R code and use
different numerical methods:

| Quantity | R code | Python reference |
|---|---|---|
| State probabilities, forward clock | closed-form transitions and Gauss-Legendre rules | products of matrix exponentials (`scipy.linalg.expm`) |
| State probabilities, reset clock | closed-form transitions and Gauss-Legendre rules | adaptive quadrature (`scipy.integrate.quad`) over the progression time, transitions by `expm` |
| Expected events, log-rank mean and information | accrual handled through cumulative integrals and Gauss-Legendre rules | adaptive quadrature with the enrolled fraction inside the integrand |
| Score covariance under the alternative | per-patient influence functions | formulas of Luo et al. (2019, Section 3.1) as printed, by nested quadrature |
| Bounds and stopping probabilities | `mvtnorm::pmvnorm()` (Miwa algorithm) | recursive conditioning with adaptive quadrature |

The influence-function derivation used in the R code was checked against the
formulas of Luo et al. (2019) in Python (function `influence_cov()` in
`ref_logrank.py`): for scenario B, the covariance at months 25 and 15 and the
variance at month 25 agree to about 1e-14 (`results/alt_check.csv`).

## Files

* `ref_model.py`: multistate model of one group.
* `ref_logrank.py`: accrual, expected events, log-rank moments, Luo et al.
  (2019) variance and covariance, and the influence-function check.
* `ref_gsd.py`: spending functions, bounds and stopping probabilities.
* `expected_values.py`: computes all expected values and writes
  `results/*.csv`.
* `simulate_ms.py`: patient-level simulation of the validation scenarios
  V6-V8, written to `../results/sim_python.csv` (see `../README.md`).

## Running

```
python3 expected_values.py            # all parts (about 30 minutes)
python3 expected_values.py gs power   # selected parts
python3 simulate_ms.py                # validation scenarios V6-V8 (about 15 minutes)
```

The results in `results/` were computed with Python 3.13, numpy 2.5.3 and
scipy 1.18.1.
