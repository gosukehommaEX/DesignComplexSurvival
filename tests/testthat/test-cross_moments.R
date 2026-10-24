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

# Expected values of version 0.3.0: inst/validation/python/results/model5_cross.csv
# (ref_cross5.py: Markov models, transition matrices and nested adaptive quadrature) and
# model5_cross_mc.csv (Monte Carlo with 10 million patients)

test_that(".cross_moments matches the reference with a switch before progression", {
  ex <- expected_csv("model5_cross.csv")
  keys <- c("psi", "phi", "dlt")
  for (name in unique(ex$scenario)) {
    scn <- scn_s3(name)
    acc <- .acc_setup(scn$accrual, 500)
    for (tp in c(12, 20)) {
      e <- ex[ex$scenario == name & ex$t_p == tp, ]
      cm <- .cross_moments(scn, acc, tp, e$t_o[1L])
      expect_lt(max(abs(cm[cbind(match(e$key_p, keys), match(e$key_o, keys))] - e$value)), 1e-11)
    }
  }
})

test_that(".cross_moments agrees with Monte Carlo with the options of version 0.3.0", {
  ex <- expected_csv("model5_cross_mc.csv")
  keys <- c("psi", "phi", "dlt")
  for (name in unique(ex$scenario)) {
    scn <- scn_s3(name)
    acc <- .acc_setup(scn$accrual, 500)
    for (tp in c(12, 20)) {
      e <- ex[ex$scenario == name & ex$t_p == tp, ]
      cm <- .cross_moments(scn, acc, tp, e$t_o[1L])
      z <- (cm[cbind(match(e$key_p, keys), match(e$key_o, keys))] - e$estimate) / e$se
      expect_lt(max(abs(z)), 4.5)
    }
  }
})

test_that(".cross_moments keeps the identity of the event indicators and the reductions", {
  # PFS followed up after a switch: an observed death is an observed PFS event
  scn <- scn_s3("PR")
  acc <- .acc_setup(scn$accrual, 500)
  ev <- .events_core(scn, acc, 20, "os")
  expect_equal(.cross_moments(scn, acc, 20, 20)["dlt", "dlt"], (ev$ev0 + ev$ev1) / 500,
               tolerance = 1e-10)
  # a switch before progression without any effect gives the model without it
  ctl <- list(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01, time_post = c(0, 3),
              death_prog = c(0.05, 0.06), dropout = 0.005)
  acc_c <- ms_accrual(time = c(0, 6), prop = c(0.3, 0.7), duration = 18)
  s0 <- ms_scenario(do.call(ms_model, ctl), models_s3("T0"), accrual = acc_c)
  s1 <- ms_scenario(do.call(ms_model, c(ctl, switch_pre = 0.03)), models_s3("T0"),
                    accrual = acc_c)
  expect_equal(.cross_moments(s1, acc, 12, 20), .cross_moments(s0, acc, 12, 20),
               tolerance = 1e-12)
})
