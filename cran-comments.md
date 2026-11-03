## Test environments

* Local: Windows, R 4.6.1
* win-builder: R-devel (2026-11-02 r90827 ucrt) and R 4.6.2 (2026-10-29 ucrt)
* GitHub Actions: macOS (R release), Windows (R release), Ubuntu (R-devel,
  R release and R oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

* Possibly misspelled words in DESCRIPTION: Bonferroni, Bretz, Luo, Maurer,
  et, al and multistate. Bonferroni, Bretz, Luo and Maurer are names of
  authors and of a procedure, "et al." is part of a reference, and
  "multistate" is the usual term for multistate models in survival analysis.

* On macOS (GitHub Actions, R 4.6.2) the check reports the installed size as
  INFO: 5.7 Mb, with 3.5 Mb in inst/validation. This directory holds the
  scripts and stored results that compare the package with published values,
  other R packages and simulations, so that users can rerun the validation.
