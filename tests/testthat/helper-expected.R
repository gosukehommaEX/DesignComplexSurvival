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
