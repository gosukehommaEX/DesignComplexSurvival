# Expected values: inst/validation/python/results/zeta_same.csv
# (ref_zeta.same_moments: cumulative functions by an ODE solver and the
# expectation over the entry time by adaptive quadrature)

check_same <- function(scn, name, n, tt, ep) {
  acc <- .acc_setup(scn$accrual, n)
  sm <- .same_moments(scn, acc, tt, ep)
  ex <- expected_csv("zeta_same.csv")
  ex <- ex[ex$scenario == name & ex$endpoint == ep, ]
  em <- ex[ex$quantity == "mean", ]
  expect_equal(sm$mean[cbind(match(em$time1, tt), match(em$key1, colnames(sm$mean)))],
               em$value, tolerance = 1e-9)
  ep2 <- ex[ex$quantity == "prod", ]
  keys <- c("psi", "phi", "dlt")
  got <- sm$prod[cbind(match(ep2$time1, tt), match(ep2$time2, tt),
                       match(ep2$key1, keys), match(ep2$key2, keys))]
  expect_equal(got, ep2$value, tolerance = 1e-9)
  sm
}

test_that(".same_moments agrees with the first-stage moments", {
  scn <- scn_b()
  acc <- .acc_setup(scn$accrual, 300)
  tt <- c(15, 25)
  sm <- .same_moments(scn, acc, tt, "os")
  core <- .lr_core(scn, acc, tt, "os")
  expect_equal(sm$mean[, "psi"], core$score_mean / 300, tolerance = 1e-12)
  expect_equal(sm$mean[, "phi"], core$info / 300, tolerance = 1e-12)
  expect_equal(sm$mean[, "dlt"], core$events / 300, tolerance = 1e-12)
  cv <- sm$prod[, , "psi", "psi"] - outer(sm$mean[, "psi"], sm$mean[, "psi"])
  expect_equal(cv, .lr_alt_cov(scn, acc, tt, "os", core$score_mean / 300),
               tolerance = 1e-12)
})

test_that(".same_moments matches the reference without progression", {
  check_same(scn_b(), "B", 300, c(15, 25), "os")
})

test_that(".same_moments matches the reference with progression and switching", {
  check_same(scn_c(), "C", 500, c(10, 20), "pfs")
  check_same(scn_c(), "C", 500, c(20, 30), "os")
})

test_that(".same_moments matches the reference with the forward clock", {
  check_same(scn_c("forward"), "Cf", 500, c(20, 30), "os")
})
