test_that(".cal_dt gives the time derivatives of the events and of the risk set", {
  # central differences of the per-patient events and of the integral of p G Y over the
  # follow-up time (step 1e-3: error of order 1e-7 relative)
  for (scn in list(scn_cal(), scn_cal(clock = "reset", tc = 13.7))) {
    acc <- .acc_setup(scn$accrual, 500)
    h <- 1e-3
    at <- function(t) {
      lk <- .cal_look(scn, acc, t)
      y <- sum(lk$g$w * (scn$prob[[1L]] * lk$G[[1L]] * lk$Y[[1L]] +
                           scn$prob[[2L]] * lk$G[[2L]] * lk$Y[[2L]]))
      c(sum(lk$ev), y)
    }
    lk <- .cal_look(scn, acc, 27.3)
    dt <- .cal_dt(scn, acc, lk)
    got <- c(sum(lk$g$w * (dt$dd[[1L]] + dt$dd[[2L]])), sum(lk$g$w * (dt$dy[[1L]] + dt$dy[[2L]])))
    expect_equal(got, (at(27.3 + h) - at(27.3 - h)) / (2 * h), tolerance = 1e-6)
  }
})
