test_that(".check_null_scenarios checks the names, allocation and accrual", {
  scn <- scn_f()
  eps <- c("pfs", "os")
  expect_true(.check_null_scenarios(NULL, scn, eps))
  expect_true(.check_null_scenarios(list(os = scn_f_null("os"), pfs = scn_f_null("pfs")),
                                    scn, eps))
  expect_error(.check_null_scenarios(scn_f_null("os"), scn, eps), "list")
  expect_error(.check_null_scenarios(list(scn_f_null("os")), scn, eps), "named")
  expect_error(.check_null_scenarios(list(os = 1), scn, eps), "list")
  s2 <- scn_f_null("os")
  s2$prob <- c(control = 0.5, treatment = 0.5)
  expect_error(.check_null_scenarios(list(os = s2), scn, eps), "allocation")
})

test_that(".check_null_scenarios checks the definition of PFS at a switch", {
  scn <- scn_s3("PRC")
  eps <- c("pfs", "os")
  expect_true(.check_null_scenarios(list(os = ms_null(scn)), scn, eps))
  expect_error(.check_null_scenarios(list(os = ms_null(scn_s3("PR"))), scn, eps), "pfs_switch")
})
