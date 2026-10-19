# Validation against independent routes

The scripts in this folder compare the analytic results of
DesignComplexSurvival with patient-level simulations and with the package
PWEALL. The expected values of the unit tests are computed separately by the
Python scripts in `python/` (see `python/README.md`).

## How to run

With the working directory at the package root:

```r
devtools::load_all()
source("inst/validation/run_all.R")
```

`python/simulate_ms.py` (scenarios V6-V8, about 15 minutes) is run in Python
beforehand; its result `results/sim_python.csv` is part of the folder and was
computed with Python 3.13, numpy 2.5.3 and scipy 1.18.1:

```
python3 inst/validation/python/simulate_ms.py
```

The FastSurvival simulation (V1-V5, about 30 minutes) is reused by later runs
as long as `helpers.R` and `validate_fastsurvival.R` are unchanged; set
`rerun_fastsurvival <- TRUE` before sourcing `run_all.R` to force a new
simulation.

## Scenarios

| Scenario | Simulated by | Content | n | Calendar looks | OS event looks |
|---|---|---|---|---|---|
| V1 | FastSurvival | test scenario A (exponential, HR 0.7) | 400 | 24, 36 | 200, 300 |
| V2 | FastSurvival | test scenario B (piecewise exponential, 2:1, piecewise accrual) | 300 | 15, 25 | 75, 150 |
| V3 | FastSurvival | innovaTV 301, exponential PFS and OS (Fleischer et al. 2009) | 482 | 27, 35 | 252, 336 |
| V4 | FastSurvival | test scenario C without switching after progression | 500 | 20, 30 | 150, 250 |
| V5 | FastSurvival | Barthel et al. (2006), Table V, row 9 (crossover in both groups) | 292 | 3, 4 years | none |
| V6 | Python | test scenario C (switching at and after progression), clock reset | 500 | 20, 30 | 150, 250 |
| V7 | Python | test scenario C, clock forward | 500 | 20, 30 | 150, 250 |
| V8 | Python | Luo et al. (2019), Example 1, delayed effect with dilution | 4600 | 30, 40 | 361, 722 |

Each scenario is simulated 20000 times under the alternative and under the
null hypothesis (`ms_null()`). Only summaries are kept. The Python simulation
uses simple randomization; FastSurvival allocates fixed group sizes and a
fixed number of patients to each accrual interval. FastSurvival switches
treatment only after an observed progression, so its latent quantities come
from a separate simulation without dropout.

## Judgements

| Quantity | Criterion |
|---|---|
| Group-wise OS and PFS events and dropouts at calendar times | `4se`: within 4 Monte Carlo standard errors |
| State occupation probabilities, OS and PFS survival without dropout | `4se_p`: within 4 standard errors computed from the package value |
| O'Brien-Fleming bounds of the event-driven looks (Python, `ref_gsd.py`) | `abs:1e-6` |
| PWEALL: expected events, mean, information and variance of the score | `rel:0.002` |
| PWEALL: power (variance under the alternative and under the null) | `abs:0.002` |
| Log-rank Z (mean, standard deviation, correlation), event-driven look times and events, power, type I error | `info` (asymptotic approximations; difference and standard error reported) |

PASS, FAIL and INFO as in `inst/reproduce/`; a FAIL whose id is registered in
`explained.csv` is EXPLAINED.

## Checks made by run_all.R

In the perturbation check the package values of the judged rows are
multiplied by 1.1. Rows for which the 10% change exceeds the tolerance plus
the observed difference (sensitive rows) must then FAIL; the number of the
other rows is reported.
The run stops with an error, after writing the outputs, when the Python
results are missing, when the perturbation check fails, or when an id of
`explained.csv` is not judged EXPLAINED.

## Files

| File | Content |
|---|---|
| `helpers.R` | Scenario definitions V1-V8, judgement and summary functions |
| `validate_fastsurvival.R` | FastSurvival simulation of V1-V5, writes `results/sim_fastsurvival.csv` |
| `python/simulate_ms.py` | Python simulation of V6-V8, writes `results/sim_python.csv` |
| `compare_simulation.R` | Simulation summaries versus the package |
| `compare_pweall.R` | PWEALL versus the package (Luo et al. 2019, Examples 1-3; Barthel et al. 2006, Table IV, rows 5 and 12) |
| `run_all.R` | Runs the comparisons and writes the summaries |
| `explained.csv` | Registered explanations (id, reason) |
| `results/` | Outputs: `compare_simulation.csv`, `compare_pweall.csv`, `summary.csv`, `summary.md`, `perturb.csv` (generated; do not edit) |
