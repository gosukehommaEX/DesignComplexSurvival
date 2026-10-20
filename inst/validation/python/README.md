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
| Moments of the score, the observed variance and the event indicator (studentized covariance, stage 2) | cumulative functions on Gauss-Legendre rules, entry by Gauss-Legendre rules | cumulative functions by an ODE solver (`solve_ivp`), entry by adaptive quadrature |
| Products of the PFS and OS functionals | four terms of counting-process integrals over the PFS and OS times (covariance note, equations (6)-(9)) | Fleischer model: conditioning on the death time; scenarios with switching: Monte Carlo with 10 million patients |
| Time derivatives of the score mean, information and events | integrals against the accrual density | central differences of the reference moments |
| Rejection probabilities of two endpoints (`gs_prob_mtp()`) | sums of rectangle probabilities over the outcomes of the first endpoint (Miwa algorithm) | numerical integration over the analyses of the joint density, with the decisions applied on the sample path (canonical joint distribution); for two looks, the formulas of Glimm et al. (2010) as one-dimensional integrals; Monte Carlo with 2 million trials |
| Designs with two endpoints or with looks driven by the events of different endpoints (`ms_power()`, stage 2) | `.look_times()`, `.design_core2()` and `.fwer_core()` with the functions above | analysis times by root finding on the reference events, moments and covariance from the references above, bounds by recursive conditioning, rejection probabilities as sums over the sample-path outcomes of rectangle probabilities by recursive conditioning with adaptive quadrature (any covariance matrix); Monte Carlo with 4 million trials |

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
* `ref_zeta.py`: moments of the score, the observed variance and the event
  indicator of one endpoint, products of the PFS and OS functionals in the
  Fleischer model, and the covariance of the standardized statistics
  (Planning/covariance_note in the project folder).
* `expected_zeta.py`: expected values for the studentized and PFS-OS
  covariance, written to `results/zeta_*.csv` and read by the tests through
  `tests/testthat/helper-expected.R`.
* `ref_mtp.py`: rejection probabilities of two endpoints under the
  Bonferroni and fixed-sequence procedures (recursive integration, Monte
  Carlo simulation and the formulas of Glimm et al. 2010 and Tamhane et
  al. 2010; the last are also used by tools/reproduce_python.py in the
  project folder).
* `expected_mtp.py`: expected values for `gs_prob_mtp()`, written to
  `results/gs_prob_mtp_configs.csv`, `results/gs_prob_mtp.csv` and the Monte
  Carlo check `results/gs_prob_mtp_mc.csv`.
* `expected_design.py`: expected values for `ms_power()` with two endpoints
  and with looks driven by the events of PFS and OS (designs D1 to D4 of
  `tests/testthat/helper-designs.R`), written to `results/design2.csv`, and
  the Monte Carlo check of the rejection probabilities
  `results/design2_mc.csv` (`ref_mtp.simulate()` with the mean vector and
  covariance matrix of each design).
* `simulate_ms.py`: patient-level simulation of the validation scenarios
  V6-V8, written to `../results/sim_python.csv` (see `../README.md`).

## Running

```
python3 expected_values.py            # all parts (about 30 minutes)
python3 expected_values.py gs power   # selected parts
python3 expected_zeta.py              # stage 2 moments (about 10 minutes)
python3 expected_mtp.py               # rejection probabilities of two endpoints (about 2 minutes)
python3 expected_design.py            # designs with two endpoints (about 30 minutes)
python3 simulate_ms.py                # validation scenarios V6-V8 (about 15 minutes)
```

The results in `results/` were computed with Python 3.13, numpy 2.5.3 and
scipy 1.18.1.
