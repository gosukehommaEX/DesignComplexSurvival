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

FastSurvival and PWEALL are used only by these scripts, so they are not
declared in DESCRIPTION. `run_all.R` records their versions in
`results/summary.md`, together with the version of lrstat (in Suggests),
which is used for the comparison with `lrstat::fseqbon()`.

`python/simulate_ms.py` (scenarios V6-V10, about 17 minutes) is run in Python
beforehand; its result `results/sim_python.csv` is part of the folder and was
computed with Python 3.13, numpy 2.5.3 and scipy 1.18.1:

```
python3 inst/validation/python/simulate_ms.py
```

The comparison of PFS and OS together uses `python/simulate_two.py` (V3, V4,
V6, V7, V11 and V12, about 15 minutes) and `python/danzer_check.py` (about 2
minutes), also run beforehand; their results `results/sim_python_two.csv` and
`results/danzer_check.csv` are part of the folder:

```
python3 inst/validation/python/simulate_two.py
python3 inst/validation/python/danzer_check.py
```

The comparison of switching from a calendar time uses
`python/simulate_calendar.py` (V13-V15, about 3 minutes), run beforehand in the
same way; its result `results/sim_python_calendar.csv` is part of the folder:

```
python3 inst/validation/python/simulate_calendar.py
```

The larger simulations of V2 and V7 under the null hypothesis use
`python/simulate_type1.py` (about 11 minutes); its result
`results/sim_python_type1.csv` is part of the folder and is not read by
`run_all.R` (see "Type I error at the looks driven by OS events" below):

```
python3 inst/validation/python/simulate_type1.py
```

The FastSurvival simulation (V1-V5, about 3 minutes) is reused by later runs
as long as `helpers.R` and `validate_fastsurvival.R` are unchanged; set
`rerun_fastsurvival <- TRUE` before sourcing `run_all.R` to force a new
simulation. The comparison with `lrstat::fseqbon()` (`compare_graph_lrstat.R`,
1.2 million simulated trials) takes about 20 minutes.

At each look driven by OS deaths, `validate_fastsurvival.R` checks that every
trial has the target number of OS deaths, and in trials 1 to 20 of the first
batch it compares the log-rank Z of `FastSurvival::analysis_fast()` with that
of `survival::survdiff()`; the run stops when a check fails, and the results
are written to `results/sim_fastsurvival_info.txt`. Up to version 0.3.2 the
death at the cut was lost in 10% to 31% of the trials, because the follow-up
time cut - entry could be smaller than its time by rounding. The correction
changed only the numbers of events at these looks, not the log-rank Z.

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
| V9 | Python | test scenario C, clock forward, death hazard after a switch measured from the switch (`switch_clock = "switch"`) | 500 | 20, 30 | 150, 250 |
| V10 | Python | test scenario C, clock reset, accelerated failure time effect of the switch (`switch_af = 0.6`) | 500 | 20, 30 | 150, 250 |

Each scenario is simulated 20000 times under the alternative and under the
null hypothesis (`ms_null()`). Only summaries are kept. The Python simulation
uses simple randomization; FastSurvival allocates fixed group sizes and a
fixed number of patients to each accrual interval. FastSurvival switches
treatment only after an observed progression, so its latent quantities come
from a separate simulation without dropout. V9 and V10 (stage 3) are defined
in `compare_simulation.R` (not in `helpers.R`, whose MD5 decides whether the
FastSurvival simulation is reused) and, independently, in
`python/simulate_ms.py`. V11-V15 are described in the sections below.

## PFS and OS together (stage 2)

`python/simulate_two.py` simulates V3, V4, V6 and V7 again in Python, and
V11 and V12 (stage 3), 20000 times under the alternative and under the global
null hypothesis, and analyses PFS and OS together at the calendar looks and
at the looks of a design with two endpoints (defined in
`compare_two_endpoints.R` and, independently, in `simulate_two.py`). V11 and
V12 are test scenario C (clock reset) with switching before progression in
the control group (`switch_pre` 0.02, `prog_switch` 0.07 and 0.06,
`death_switch_pre` 0.008); in V12 PFS is censored at that switch
(`pfs_switch = "censor"`) and both groups have a censoring of PFS only
(`dropout_pfs` 0.01). Their patients are generated by `ref_model5.py`.

| Scenario | Looks | Tested | Order and spending |
|---|---|---|---|
| V3 | 252 and 336 OS events | PFS and OS at both looks | OS (O'Brien-Fleming type), then PFS (Pocock type) |
| V4, V6, V7 | 330 PFS events, then 250 OS events | PFS and OS at look 1, OS at look 2 | PFS, then OS (both O'Brien-Fleming type) |
| V11, V12 | calendar months 20 and 30 | PFS and OS at look 1, OS at look 2 | PFS, then OS (both O'Brien-Fleming type) |

Each design is analysed with the fixed-sequence procedure (`hierarchical =
"overall"` and `"stagewise"`, each endpoint at level 0.025), with the
Bonferroni procedure (weights 0.8 and 0.2) and with the same Bonferroni
procedure with recycling of the level (`transition = c(1, 1)`; stage 3). The
script also estimates, from 10 million simulated patients, the per-patient
covariances of PFS and OS at the calendar looks under the alternative (V3,
V4, V6 and V7), and for V3 computes them deterministically with
`ref_zeta.py`. `python/danzer_check.py` applies the
covariance estimator of Danzer et al. (2025, Section 2) to trials of 1,000,000
patients under the global null hypothesis (V3 and V6, 10 trials each). V5 and
V8 are not included, as their progression is crossover or discontinuation of
treatment.

## Switching from a calendar time (stage 3)

`python/simulate_calendar.py` simulates test scenario C with switching
possible from calendar month 15 (`switch_start = 15` in `ms_scenario()`):

| Scenario | Content | n | Calendar looks | OS event looks |
|---|---|---|---|---|
| V13 | clock forward | 500 | 24, 36 | 150, 250 |
| V14 | clock reset | 500 | 24, 36 | 150, 250 |
| V15 | clock forward, accelerated failure time effect of the switch (`switch_af = 0.6`) | 500 | 24, 36 | 150, 250 |

The per-patient OS events, score mean and information at the calendar looks
are estimated from 100 trials of 20000 patients; Z, its correlation between
the looks, the look times and the rejection rates of the O'Brien-Fleming
design at the event looks from 20000 trials of 500 patients, under the
alternative and under the null hypothesis. `compare_calendar.R` defines the
same scenarios and compares the results with the package.

## Graphical procedure with recycling (stage 3)

`compare_graph_lrstat.R` checks the rejection probabilities of the
Bonferroni procedure with recycling of the level (`ms_test(transition =
...)`, `gs_prob_mtp(upper_recycled = ...)`) against the decisions of
`lrstat::fseqbon()` (lrstat 0.3.4 or later, with `lookback = FALSE`). For
each design, 100000 sets of statistics are drawn from the multivariate
normal distribution with the mean and covariance of `ms_power()`, and with
the mean of both endpoints or of one of them set to 0 (global and partial
null hypotheses). `fseqbon()` has no futility bounds; in L3 the
non-binding futility stop of the package is applied to its decisions
afterwards. The decisions are also compared trial by trial with the
procedure applied directly with the bounds of the package, and the bounds
with `lrstat::getBound()`.

| Design | Scenario | Looks | Order and spending | Weights | `transition` |
|---|---|---|---|---|---|
| L1 | V3 | 252 and 336 OS events, both endpoints at both looks | OS (O'Brien-Fleming type), then PFS (Pocock type) | 0.8, 0.2 | c(1, 1) |
| L2 | V6 | 330 PFS events (PFS and OS), then 250 OS events (OS) | PFS, then OS (both O'Brien-Fleming type) | 0.8, 0.2 | c(1, 1) |
| L3 | V3 | as L1 | OS (O'Brien-Fleming type, non-binding futility bound by O'Brien-Fleming type beta spending, beta 0.2), then PFS (Hwang-Shih-DeCani, gamma -4) | 0.6, 0.4 | c(0.5, 1) |

## Type I error at the looks driven by OS events

Under the null hypothesis, the simulated rejection rates of the
O'Brien-Fleming design at the looks driven by OS events (`info` rows of
`results/summary.csv`) were 0.0287 in V2 and 0.0285 in V7, about 3 standard
errors above the level 0.025 computed by the package, while those of the other
10 scenarios with such looks were within 2 standard errors of it.
`python/simulate_type1.py` repeats the two simulations with new seeds and more
trials (standard errors at most 0.0004):

| Simulation | Trials | Rejection in the direction of benefit | Rejection in the opposite direction |
|---|---|---|---|
| V7 | 200000 | 0.025 | 0.025 |
| Scenario B with the group sizes of V2 (100 control, 200 treatment) | 400000 | 0.028 | 0.023 |
| Scenario B with 200 control and 100 treatment | 400000 | 0.023 | 0.028 |
| Scenario B with 150 in each group | 400000 | 0.025 | 0.025 |

The rate of V7 in `results/summary.csv` is thus a chance deviation. The rate
of V2 is a property of the log-rank test with unequal allocation in small
trials, which the normal approximation of the package does not describe: at
each death the indicator of the treatment group has the third cumulant
q (1 - q) (1 - 2 q), with q the proportion of the treatment group in the risk
set, so that after d deaths the skewness of Z is about
(2 q - 1) / sqrt(q (1 - q) d), 0.08 and 0.06 for q = 2/3 and 75 and 150
deaths (simulated: 0.08 and 0.06). With the treatment group the larger group
this raises the probability of rejecting in the direction of benefit, and with
the control group the larger group it lowers it; the effect decreases as the
number of deaths increases.

## Judgements

| Quantity | Criterion |
|---|---|
| Group-wise OS and PFS events and dropouts at calendar times | `4se`: within 4 Monte Carlo standard errors |
| State occupation probabilities, OS and PFS survival without dropout | `4se_p`: within 4 standard errors computed from the package value |
| O'Brien-Fleming bounds of the event-driven looks (Python, `ref_gsd.py`) | `abs:1e-6` |
| PWEALL: expected events, mean, information and variance of the score | `rel:0.002` |
| PWEALL: power (variance under the alternative and under the null) | `abs:0.002` |
| Log-rank Z (mean, standard deviation, correlation), event-driven look times and events, power, type I error | `info` (asymptotic approximations; difference and standard error reported) |
| PFS and OS: per-patient covariances at the calendar looks (10 million simulated patients) versus `lr_cov()` with each `variance` | `4se` |
| PFS and OS in V3: the same covariances computed by `ref_zeta.py` | `rel:1e-6` |
| PFS and OS: efficacy bounds of the designs (`ref_gsd.py`) | `abs:1e-6` |
| PFS and OS: correlations of Z at the calendar looks, look times, events, Z, rejection rates of the procedures (power and familywise error rate), estimate of Danzer et al. (2025) | `info` |
| Switching from a calendar time: per-patient OS events by group and in total, score mean and information at the calendar looks versus `ms_events()` and `lr_moments()` | `4se` |
| Switching from a calendar time: O'Brien-Fleming bounds of the event looks (`ref_gsd.py`) | `abs:1e-6` |
| Switching from a calendar time: Z (mean, standard deviation, correlation), look times and rejection rates | `info` |
| Recycling: probabilities of `gs_prob_mtp()` and `ms_power()` of rejecting each hypothesis at each look, the first, the second, both and at least one, and of stopping for futility, versus the proportions of the trials rejected by `lrstat::fseqbon()` | `4se_p` |
| Recycling: efficacy bounds at the initial and at the recycled levels versus `lrstat::getBound()` | `abs:1e-5` |
| Recycling: proportion of the trials in which the decisions of `fseqbon()` differ from the procedure applied directly with the bounds of the package | `abs:1e-4` |

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
| `validate_fastsurvival.R` | FastSurvival simulation of V1-V5, writes `results/sim_fastsurvival.csv` and `results/sim_fastsurvival_info.txt` (versions, times and checks) |
| `python/simulate_ms.py` | Python simulation of V6-V10, writes `results/sim_python.csv` |
| `compare_simulation.R` | Simulation summaries versus the package (definitions of V9 and V10) |
| `python/simulate_two.py` | Python simulation of PFS and OS together (V3, V4, V6, V7, V11, V12), writes `results/sim_python_two.csv` |
| `python/danzer_check.py` | Estimator of Danzer et al. (2025) under the global null hypothesis, writes `results/danzer_check.csv` |
| `compare_two_endpoints.R` | Designs with two endpoints (definitions of V11 and V12); PFS and OS results versus the package |
| `python/simulate_calendar.py` | Python simulation of switching from a calendar time (V13-V15), writes `results/sim_python_calendar.csv` |
| `compare_calendar.R` | Definitions of V13-V15; switching from a calendar time versus the package |
| `python/simulate_type1.py` | Larger simulations of V2 and V7 under the null hypothesis, writes `results/sim_python_type1.csv` (not read by `run_all.R`) |
| `compare_pweall.R` | PWEALL versus the package (Luo et al. 2019, Examples 1-3; Barthel et al. 2006, Table IV, rows 5 and 12) |
| `compare_graph_lrstat.R` | Graphical procedure with recycling of the level versus `lrstat::fseqbon()` and `lrstat::getBound()` |
| `run_all.R` | Runs the comparisons and writes the summaries |
| `explained.csv` | Registered explanations (id, reason) |
| `results/` | Outputs: `compare_simulation.csv`, `compare_two_endpoints.csv`, `compare_calendar.csv`, `compare_pweall.csv`, `compare_graph_lrstat.csv`, `summary.csv`, `summary.md`, `perturb.csv` (generated; do not edit) |
