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
