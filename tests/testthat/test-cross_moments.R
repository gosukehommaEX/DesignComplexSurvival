# Expected values: inst/validation/python/results/zeta_cross.csv (scenario F,
# ref_zeta.cross_moments_fleischer: conditioning on the death time, ODE-based
# cumulative functions and nested adaptive quadrature) and zeta_cross_mc.csv
# (Monte Carlo with 10 million patients for scenarios C and Cf)

test_that(".cross_moments matches the reference in a Fleischer model", {
  scn <- scn_f()
  acc <- .acc_setup(scn$accrual, 400)
  ex <- expected_csv("zeta_cross.csv")
  keys <- c("psi", "phi", "dlt")
  for (tp in c(12, 20)) {
    for (to in c(12, 20)) {
      cm <- .cross_moments(scn, acc, tp, to)
      e <- ex[ex$t_p == tp & ex$t_o == to, ]
      expect_equal(cm[cbind(match(e$key_p, keys), match(e$key_o, keys))], e$value,
                   tolerance = 1e-9)
    }
  }
})

test_that(".cross_moments satisfies the identity of the event indicators", {
  # at the same time an observed death is an observed PFS event
  scn <- scn_c()
  acc <- .acc_setup(scn$accrual, 500)
  cm <- .cross_moments(scn, acc, 20, 20)
  ev <- .events_core(scn, acc, 20, "os")
  expect_equal(cm["dlt", "dlt"], (ev$ev0 + ev$ev1) / 500, tolerance = 1e-10)
  expect_equal(unname(.cross_moments(scn, acc, 0, 20)), matrix(0, 3, 3))
})

test_that(".cross_moments agrees with Monte Carlo with switching", {
  ex <- expected_csv("zeta_cross_mc.csv")
  keys <- c("psi", "phi", "dlt")
  for (clock in c("reset", "forward")) {
    scn <- scn_c(clock)
    acc <- .acc_setup(scn$accrual, 500)
    cm <- .cross_moments(scn, acc, 12, 20)
    e <- ex[ex$scenario == if (clock == "reset") "C" else "Cf", ]
    z <- (cm[cbind(match(e$key_p, keys), match(e$key_o, keys))] - e$estimate) / e$se
    expect_lt(max(abs(z)), 4.5)
  }
})
