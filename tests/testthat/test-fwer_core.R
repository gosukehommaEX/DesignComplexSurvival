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

test_that(".fwer_core applies the bounds after recycling", {
  scn <- scn_f()
  acc <- .acc_setup(scn$accrual, 400)
  lo <- ms_looks("os")
  te <- ms_test(endpoints = c("pfs", "os"), mtp = "bonferroni", weights = c(0.6, 0.4),
                transition = c(1, 1))
  tl <- .tested_looks(lo, te$endpoints)
  bnd <- list(gs_bounds(1, alpha = 0.015), gs_bounds(1, alpha = 0.01))
  ur <- list(gs_bounds(1, alpha = 0.025)$upper, gs_bounds(1, alpha = 0.025)$upper)
  fw0 <- .fwer_core(scn, acc, .looks_table(lo), 220, tl, te, bnd)
  fw <- .fwer_core(scn, acc, .looks_table(lo), 220, tl, te, bnd, upper_recycled = ur)
  # the probability of rejecting at least one hypothesis does not change
  expect_equal(fw$type1, fw0$type1, tolerance = 1e-8)
  # P(reject PFS) = P(Z1 >= a1) + P(b <= Z1 < a1, Z2 >= a2), by conditioning
  # on Z1
  s0 <- ms_null(scn)
  t0 <- .event_times(s0, acc, 220, "os")
  r <- .zeta_cov(s0, acc, c("pfs", "os"), c(t0, t0), NULL, "null")$sigma[1, 2]
  a1 <- bnd[[1]]$upper
  a2 <- bnd[[2]]$upper
  f <- function(x) stats::dnorm(x) * stats::pnorm((r * x - a2) / sqrt(1 - r^2))
  extra <- stats::integrate(f, ur[[1]], a1, rel.tol = 1e-12)$value
  expect_lt(abs(fw$cum_alpha[1] - (0.015 + extra)), 1e-8)
  expect_gt(extra, 1e-4)
})
