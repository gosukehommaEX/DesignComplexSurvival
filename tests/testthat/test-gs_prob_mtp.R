# Expected values: inst/validation/python/results/gs_prob_mtp.csv, computed by
# ref_mtp.recursion() (numerical integration over the analyses of the joint
# density, with the decisions of the procedure applied on the sample path)
# and, for two looks, by the formulas of Glimm et al. (2010, equations (1)
# and (2)) written as one-dimensional integrals (ref_mtp.glimm_two_look()).
# Configurations: gs_prob_mtp_configs.csv.

mtp_check_config <- function(id, method, mtp = NULL) {
  cf <- expected_csv("gs_prob_mtp_configs.csv")
  ex <- expected_csv("gs_prob_mtp.csv")
  cf <- cf[cf$config == id, ]
  ex <- ex[ex$config == id & ex$method == method, ]
  if (!is.null(mtp)) {
    ex <- ex[paste(ex$mtp, ex$hierarchical) == mtp, ]
  }
  runs <- unique(ex[, c("mtp", "hierarchical", "hypothesis")])
  a1 <- mtp_num(cf$a1)
  a2 <- mtp_num(cf$a2)
  n <- 0L
  for (r in seq_len(nrow(runs))) {
    e <- ex[ex$mtp == runs$mtp[r] & ex$hierarchical == runs$hierarchical[r] &
              ex$hypothesis == runs$hypothesis[r], ]
    cn <- mtp_canonical(mtp_num(cf$tau1), mtp_num(cf$tau2), cf$rho, e$theta1[1],
                        e$theta2[1], a1, a2)
    res <- gs_prob_mtp(list(mtp_num(cf$u1), mtp_num(cf$u2)), lower = mtp_num(cf$l1),
                       mean = cn$mean, sigma = cn$sigma, analysis = list(a1, a2),
                       mtp = runs$mtp[r], hierarchical = runs$hierarchical[r])
    q <- mtp_quantities(res)
    expect_true(all(e$quantity %in% names(q)))
    expect_lt(max(abs(q[e$quantity] - e$value)), 1e-6)
    n <- n + nrow(e)
  }
  n
}

# three looks of each endpoint and futility bounds: 3 hypotheses, 13 quantities
test_that("gs_prob_mtp matches the recursion with three looks (overall)", {
  expect_identical(mtp_check_config("three_looks", "recursion", "fixed_sequence overall"), 39L)
})

test_that("gs_prob_mtp matches the recursion with three looks (stagewise)", {
  expect_identical(mtp_check_config("three_looks", "recursion", "fixed_sequence stagewise"), 39L)
})

test_that("gs_prob_mtp matches the recursion with three looks (Bonferroni)", {
  expect_identical(mtp_check_config("three_looks", "recursion", "bonferroni overall"), 39L)
})

test_that("gs_prob_mtp matches the recursion with different analyses", {
  expect_identical(mtp_check_config("unequal", "recursion"), 99L)
})

test_that("gs_prob_mtp matches the recursion when the second endpoint starts later", {
  expect_identical(mtp_check_config("late_second", "recursion"), 108L)
})

test_that("gs_prob_mtp matches the formulas of Glimm et al. for two looks", {
  # strategies (a) and (b): 5 and 3 quantities for each correlation
  n <- mtp_check_config("glimm_rho0", "glimm") + mtp_check_config("glimm_rho08", "glimm")
  expect_identical(n, 16L)
  n <- mtp_check_config("glimm_rho0", "recursion") + mtp_check_config("glimm_rho08", "recursion")
  expect_identical(n, 60L)
})

test_that("gs_prob_mtp agrees with gs_prob for the marginal probabilities", {
  cn <- mtp_canonical(c(0.4, 0.7, 1), c(0.4, 0.7, 1), 0.6, 3, 2.5, 1:3, 1:3)
  u1 <- c(3, 2.4, 2)
  u2 <- c(2.6, 2.3, 2.1)
  l1 <- c(0, 1, 2)
  m1 <- gs_prob(l1, u1, cn$mean[1:3], cn$sigma[1:3, 1:3])
  m2 <- gs_prob(rep(-Inf, 3), u2, cn$mean[4:6], cn$sigma[4:6, 4:6])
  for (h in c("overall", "stagewise")) {
    r <- gs_prob_mtp(list(u1, u2), l1, cn$mean, cn$sigma, hierarchical = h)
    expect_equal(r$looks$prob_upper[1:3], m1$prob_upper, tolerance = 1e-8)
    expect_equal(r$looks$prob_lower[1:2], m1$prob_lower[1:2], tolerance = 1e-8)
    expect_equal(unname(r$reject["any"]),
                 unname(r$reject["first"] + r$reject["second"] - r$reject["both"]))
  }
  # Bonferroni without futility: the second endpoint is tested at all looks
  rb <- gs_prob_mtp(list(u1, u2), NULL, cn$mean, cn$sigma, mtp = "bonferroni")
  expect_equal(rb$looks$prob_upper[4:6], m2$prob_upper, tolerance = 1e-7)
  # fixed sequence with the first hypothesis always rejected at the first look
  ro <- gs_prob_mtp(list(c(-Inf, u1[2:3]), u2), NULL, cn$mean, cn$sigma)
  expect_equal(ro$looks$prob_upper[4:6], m2$prob_upper, tolerance = 1e-7)
  expect_equal(unname(ro$reject["both"]), unname(ro$reject["second"]), tolerance = 1e-7)
  # rejection of the second hypothesis: stagewise <= overall <= Bonferroni
  # (Glimm et al. 2010, inequality (5))
  rs <- gs_prob_mtp(list(u1, u2), NULL, cn$mean, cn$sigma, hierarchical = "stagewise")
  rv <- gs_prob_mtp(list(u1, u2), NULL, cn$mean, cn$sigma, hierarchical = "overall")
  expect_lt(rs$reject[["second"]], rv$reject[["second"]])
  expect_lt(rv$reject[["second"]], rb$reject[["second"]])
  # an infinite bound means no test at that look
  ri <- gs_prob_mtp(list(u1, c(Inf, u2[2:3])), NULL, cn$mean, cn$sigma)
  rj <- gs_prob_mtp(list(u1, u2[2:3]), NULL, cn$mean[-4], cn$sigma[-4, -4],
                    analysis = list(1:3, 2:3))
  expect_equal(ri$reject, rj$reject, tolerance = 1e-8)
  expect_identical(ri$looks$prob_upper[4], 0)
})

test_that("gs_prob_mtp checks its arguments", {
  s <- diag(4)
  expect_error(gs_prob_mtp(c(2, 2), sigma = s), "list of two")
  expect_error(gs_prob_mtp(list(c(2, 2), 2), sigma = diag(3)), "analysis")
  expect_error(gs_prob_mtp(list(c(2, 2), c(2, 2)), lower = c(3, 0), sigma = s), "exceed")
  expect_error(gs_prob_mtp(list(c(2, 2), c(2, 2)), lower = 0, sigma = s), "lower")
  expect_error(gs_prob_mtp(list(c(2, 2), c(2, 2)), mean = 1:3, sigma = s), "mean")
  expect_error(gs_prob_mtp(list(c(2, 2), c(2, 2)), sigma = diag(3)), "sigma")
  expect_error(gs_prob_mtp(list(c(2, 2), c(2, 2)), sigma = s, analysis = list(c(2, 1), 1:2)),
               "analysis")
  expect_error(gs_prob_mtp(list(rep(2, 11), rep(2, 10)), sigma = diag(21),
                           analysis = list(1:11, 1:10)), "20 looks")
  expect_error(gs_prob_mtp(list(c(2, 2), c(2, 2)), sigma = s, mtp = "holm"))
})
