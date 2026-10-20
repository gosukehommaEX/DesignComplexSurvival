test_that(".fwer_core gives the familywise error rate of a single analysis", {
  scn <- scn_f()
  acc <- .acc_setup(scn$accrual, 400)
  lo <- ms_looks("os")
  te <- ms_test(endpoints = c("pfs", "os"), mtp = "bonferroni", weights = c(0.6, 0.4))
  bnd <- list(gs_bounds(1, alpha = 0.015), gs_bounds(1, alpha = 0.01))
  fw <- .fwer_core(scn, acc, .looks_table(lo), 220, .tested_looks(lo, te$endpoints), te, bnd)
  s0 <- ms_null(scn)
  t0 <- .event_times(s0, acc, 220, "os")
  r <- .zeta_cov(s0, acc, c("pfs", "os"), c(t0, t0), NULL, "null")$sigma[1, 2]
  # P(Z1 < a, Z2 < b) by conditioning on Z1
  a <- bnd[[1]]$upper
  b <- bnd[[2]]$upper
  f <- function(x) stats::dnorm(x) * stats::pnorm((b - r * x) / sqrt(1 - r^2))
  p_none <- stats::integrate(f, -Inf, a, rel.tol = 1e-12)$value
  expect_lt(abs(fw$type1[["nonbinding"]] - (1 - p_none)), 1e-8)
  expect_equal(fw$type1[["binding"]], fw$type1[["nonbinding"]])
  expect_identical(fw$fwer$scenario, "global")
  expect_equal(fw$cum_alpha, c(0.015, 0.01), tolerance = 1e-8)
  fq <- .fwer_core(scn, acc, .looks_table(lo), 220, .tested_looks(lo, te$endpoints), te, bnd,
                   list(os = scn_f_null("os")))
  expect_identical(fq$fwer$null_endpoints, c("PFS, OS", "OS"))
  expect_lt(abs(fq$fwer$nonbinding[2] - 0.01), 1e-6)
})
