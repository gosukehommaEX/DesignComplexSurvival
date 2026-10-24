# Expected values: inst/validation/python/results/model5_states.csv

test_that(".ep_eval gives the PFS of state 1 when a switch censors PFS", {
  ex <- expected_csv("model5_states.csv")
  ex <- ex[ex$pfs_switch == "censor", ]
  for (name in unique(ex$model)) {
    for (clock in c("reset", "forward")) {
      e <- ex[ex$model == name & ex$clock == clock, ]
      m <- models_s3(name, clock)
      m$pfs_switch <- "censor"
      ep <- .ep_eval(m, e$time, "pfs")
      expect_lt(max(abs(ep$surv - e$pfs_surv)), 1e-10)
      expect_lt(max(abs(ep$dens - e$pfs_dens)), 1e-10)
      expect_equal(ep$surv, e$p1, tolerance = 1e-10)
    }
  }
})

test_that(".ep_eval without a switch before progression does not depend on pfs_switch", {
  m <- models_c()$control
  m2 <- m
  m2$pfs_switch <- "censor"
  expect_identical(.ep_eval(m, c(1, 5, 20), "pfs"), .ep_eval(m2, c(1, 5, 20), "pfs"))
  expect_identical(.ep_eval(m, c(1, 5, 20), "os"), .ep_eval(m2, c(1, 5, 20), "os"))
})
