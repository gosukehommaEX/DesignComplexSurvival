test_that(".grp_eval adds the censoring of PFS only to the dropout of PFS", {
  s <- c(0.5, 3, 7, 20)
  m <- models_s3("P", dropout_pfs = c(0.02, 0.03))
  m0 <- models_s3("P")
  m1 <- models_s3("P", dropout = c(0.025, 0.035))
  a <- .grp_eval(m, s, "pfs")
  b <- .grp_eval(m1, s, "pfs")
  for (k in c("S", "f", "R", "e", "c")) {
    expect_equal(a[[k]], b[[k]], tolerance = 1e-14)
  }
  expect_identical(.grp_eval(m, s, "os"), .grp_eval(m0, s, "os"))
})

test_that(".grp_eval treats a censoring switch as a censoring of PFS", {
  # with pfs_switch = "censor", the switch acts on PFS as a dropout with hazard switch_pre
  s <- c(0.5, 3, 7, 20)
  m <- scn_s3("PRC")$control
  m1 <- models_s3("P", switch_pre = 0, dropout = 0.005 + c(0.03, 0.02))
  a <- .grp_eval(m, s, "pfs")
  b <- .grp_eval(m1, s, "pfs")
  for (k in c("R", "e", "c")) {
    expect_equal(a[[k]], b[[k]], tolerance = 1e-14)
  }
  expect_equal(a$f / a$S, b$f / b$S, tolerance = 1e-14)
  expect_equal(.grp_eval(scn_s3("PR")$control, s, "os"), .grp_eval(m, s, "os"))
})
