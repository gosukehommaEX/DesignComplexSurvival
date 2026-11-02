# Expected values computed by the independent Python implementation in
# inst/validation/python (see the README there); read directly from the CSV
# files so that the values are not copied by hand.
expected_csv <- function(name) {
  path <- system.file("validation", "python", "results", name,
                      package = "DesignComplexSurvival")
  if (!nzchar(path)) {
    stop("Expected values not found: ", name)
  }
  read.csv(path, stringsAsFactors = FALSE)
}

# Covariance matrix of calendar_zeta.csv (switching from a calendar time) for a design ("alt" or
# "null"), a trigger ("none" or "os") and a variance
expected_zeta_cal <- function(design, trigger, variance) {
  ex <- expected_csv("calendar_zeta.csv")
  e <- ex[ex$design == design & ex$trigger == trigger & ex$variance == variance, ]
  matrix(e$value[order(e$j, e$i)], max(e$i), max(e$j))
}

# Expected values of version 0.4.0 (expected_stage4.csv, expected_stage4.py) of a design (config),
# a quantity and a hypothesis ("design" or "null"): the expected values (by_look = FALSE) or the
# values of the looks under the null scenario in the order of the looks (by_look = TRUE)
stage4_expected <- function(config, quantity, hypothesis, by_look = FALSE) {
  ex <- expected_csv("expected_stage4.csv")
  e <- ex[ex$config == config & ex$quantity == quantity & ex$hypothesis == hypothesis &
            is.na(ex$look) != by_look, ]
  if (nrow(e) == 0L) {
    stop("No expected values for ", config, ", ", quantity, ", ", hypothesis)
  }
  e$value[order(e$look)]
}

# Conditional power of version 0.4.0 (cond_power_stage4.csv, expected_stage4.py) of a case and a
# quantity: a matrix with the rows i and columns j (quantities with both indices), the values in
# the order of i, or a single value
stage4_cp <- function(case, quantity) {
  ex <- expected_csv("cond_power_stage4.csv")
  e <- ex[ex$case == case & ex$quantity == quantity, ]
  if (nrow(e) == 0L) {
    stop("No conditional power values for ", case, ", ", quantity)
  }
  if (!anyNA(e$j)) {
    return(matrix(e$value[order(e$j, e$i)], max(e$i), max(e$j)))
  }
  e$value[order(e$i)]
}
