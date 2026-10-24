# Expected values: inst/validation/python/results/model5_from5.csv (ref_model5.from5:
# adaptive quadrature over the time of progression)

test_that(".from_state5 matches the reference values", {
  ex <- expected_csv("model5_from5.csv")
  for (name in unique(ex$model)) {
    for (clock in c("reset", "forward")) {
      e <- ex[ex$model == name & ex$clock == clock, ]
      f <- .from_state5(models_s3(name, clock), e$w, e$u)
      expect_lt(max(abs(f$surv - e$surv)), 1e-12)
      expect_lt(max(abs(f$dens - e$dens)), 1e-12)
    }
  }
})

test_that(".from_state5 is exponential without progression after the switch", {
  m <- models_s3("P", prog_switch = 0)
  f <- .from_state5(m, c(0, 2), c(5, 30))
  expect_equal(f$surv, exp(-0.007 * c(5, 28)), tolerance = 1e-14)
  expect_equal(f$dens, 0.007 * f$surv, tolerance = 1e-14)
})
