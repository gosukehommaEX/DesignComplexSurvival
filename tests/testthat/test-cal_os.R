test_that(".cal_os uses .os_avail for a group that switches and .ep_eval otherwise", {
  m <- models_c()
  s <- c(3, 12, 25)
  o <- c(5, 0, 10)
  expect_identical(.cal_os(m$control, s, o), .os_avail(m$control, s, o))
  ep <- .ep_eval(m$treatment, s, "os")
  expect_identical(.cal_os(m$treatment, s, o), list(surv = ep$surv, dens = ep$dens))
})
