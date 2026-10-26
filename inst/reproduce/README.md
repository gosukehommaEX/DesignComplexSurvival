# Reproduction of published results

The scripts in this folder recompute published sample sizes, numbers of
events, analysis times, bounds and powers with DesignComplexSurvival and
compare them with the values printed in the articles.

## How to run

With the working directory at the package root:

```r
devtools::load_all()
source("inst/reproduce/run_all.R")
```

A single script can also be run, for example
`source("inst/reproduce/reproduce_Barthel2006.R")`; it writes its own
results file.

## Files

| File | Content |
|---|---|
| `published_values.csv` | Published values: id, source, file and page, location, quantity, printed value, search pattern, rounding, criterion, explained_id, note |
| `explained.csv` | Registered explanations: id, source, reason, evidence and the judgement expected in the perturbation check |
| `helpers.R` | Judgement, output, common model parts and the SIZE algorithm (discrete-time Markov chain of Lakatos 1988) |
| `reproduce_<Source>.R` | One script per source |
| `run_all.R` | Runs all scripts and writes the summaries |
| `results/` | Outputs (generated; do not edit) |

| Script | Source |
|---|---|
| `reproduce_Barthel2006.R` | Barthel et al. (2006), Statistics in Medicine 25:2521-2542, Tables II-V (ART) |
| `reproduce_Erdmann2025.R` | Erdmann et al. (2025), Biometrical Journal 67:e70017, Tables 2 and 3 (illness-death model; settings of the authors' R code) |
| `reproduce_Glimm2010.R` | Glimm et al. (2010), Statistics in Medicine 29:219-228, Sections 2-5 (strategies (a), (b) and (d) by `gs_prob_mtp()`) |
| `reproduce_Lakatos1988.R` | Lakatos (1988), Biometrics 44:229-241, Section 3 example and Tables 2-3 |
| `reproduce_Lakatos2002.R` | Lakatos (2002), Statistics in Medicine 21:1969-1989, Tables VIII-IX |
| `reproduce_Lu2021.R` | Lu (2021), Pharmaceutical Statistics, doi:10.1002/pst.2069, Tables 1-3 (lrstat) |
| `reproduce_Luo2019.R` | Luo et al. (2019), Statistics in Medicine 38:192-209, Examples 1-3 (PWEALL) |
| `reproduce_MaurerBretz2013.R` | Maurer and Bretz (2013), Statistics in Biopharmaceutical Research 5:311-320, Tables 1 and 2 and Section 4 (nominal levels and repeated p-values, ADDPLAN 6.0.1) |
| `reproduce_Porcher2002.R` | Porcher et al. (2002), Controlled Clinical Trials 23:650-661, Tables 1-2 and the MAG90 example |
| `reproduce_Shih1995.R` | Shih (1995), Controlled Clinical Trials 16:395-407, examples computed with SIZE |
| `reproduce_Sugimoto2017.R` | Sugimoto et al. (2017), Statistics in Medicine 36:1363-1382, Tables III and VI (PFS and OS, rows with rho = 0) |
| `reproduce_Tamhane2010.R` | Tamhane et al. (2010), Biometrics 66:1174-1184, Table 1 (by `gs_prob_mtp()`), Sections 3.2 and 5 |
| `reproduce_Tang2022.R` | Tang (2022), Statistics in Medicine, doi:10.1002/sim.9256, Example 4 |
| `reproduce_TorresSaavedra2026.R` | Torres-Saavedra et al. (2026), Clinical Trials 23:219-224, Simulations and Results sections (crossover at progression in the rank preserving structural failure time model; probability of progression and proportions of censored observations) |
| `reproduce_Vergote2024.R` | innovaTV 301 statistical analysis plan, version 2.0 (Vergote et al. 2024, doi:10.1056/NEJMoa2313811), Sections 6.2 and 6.7 (EAST) |
| `reproduce_Xia2016.R` | Xia et al. (2016), Statistics in Biopharmaceutical Research 8:12-21, CALGB 30607 example |

## Judgements

Each published value is judged by the criterion in `published_values.csv`
(P is the published value and x the computed value):

| Criterion | Pass when |
|---|---|
| `ceil` | P - 1.05 < x <= P + 0.05 (values rounded up) |
| `even_ceil` | P - 2.05 < x <= P + 0.05 (values rounded up to even numbers) |
| `ceil_to:step` | P - 1.05 step < x <= P + 0.05 step (values rounded up to a multiple of step) |
| `abs:tol` | abs(x - P) <= tol |
| `power_at_time:target:tol` | the power at the published time P is within tol of target |
| `round_n_from_d` | x rounded to the nearest integer equals P |
| `info` | not judged |

- PASS: the criterion is met.
- EXPLAINED: the criterion is not met, the row has an `explained_id`, and the
  evidence registered in `explained.csv`, computed by the script, meets the
  criterion (and any additional condition stated there).
- FAIL: the criterion is not met and there is no registered explanation, or
  the evidence does not meet the criterion.
- INFO: the published value was computed by a method that the package does
  not reproduce (for example simulation or a discrete-time approximation);
  the difference is reported without a judgement.

Where the article prints an input or a heading that appears to be a misprint
(Tang2022-Ex4-lag, Luo2019-Ex2-labels), the rows are computed as printed and
the corrected reading is the evidence. The same applies to the values of
Table 1 of Tamhane et al. (2010) that appear to be misprints
(Tamhane2010-T1-OF1OF2-rho0.4-d and -delta,
Tamhane2010-T1-OF1AH2-rho0.4-delta).

## Checks made by run_all.R

The run stops with an error, after writing the outputs, when a published
value does not have exactly one result row, when an `explained_id` is not
registered, when a registered explanation never gives EXPLAINED, when a
diagnostic is outside its tolerance, or when the perturbation check fails.

In the perturbation check, each judged published value is multiplied by 1.1,
one row at a time, and judged again; for the criterion `power_at_time` the
power is recomputed at 1.1 times the published time. A value printed with one
significant digit, whose change by 10 per cent is within the tolerance of its
`abs` criterion (MaurerBretz2013-T1-a1-0.0125, -T1-a2-0.01875 and
-S4-a1-0.025), is moved by three times the tolerance instead. Every row is expected
to FAIL, except Barthel2006-T2-r2 and Tamhane2010-T1-OF1OF2-rho0.4-d, whose
evidence is compared with a value printed in another row (the same design in
Barthel et al.; the Delta_1 of the same cell in Tamhane et al.) and which are
therefore expected to stay EXPLAINED.

## Outputs (results/)

| File | Content |
|---|---|
| `<Source>.csv` | All rows of a source: published and computed values, difference, evidence, judgement, perturbed values, note |
| `summary.csv` | All rows of all sources |
| `summary.md` | Counts by source and judgement, FAIL, EXPLAINED (with reasons) and INFO rows, perturbation check, diagnostics, run times |
| `perturb.csv` | Judgements of the perturbation check |
| `diagnostics.csv` | Checks of the approximation of simultaneous entry by uniform accrual over 1e-6 times the study length |
| `helpers_<Source>.csv` | Values of the article-specific algorithms (SIZE chain, the K-interval formula of Porcher et al., formula (2) of Sugimoto et al.) for an independent check |
