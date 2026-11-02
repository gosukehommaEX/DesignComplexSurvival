# Expected values: inst/validation/python/results/cond_power_stage4.csv (cases cm_*; conditional
# moments from the precision matrix in expected_stage4.py)

test_that(".cond_moments gives the conditional mean and covariance", {
  mean <- c(0.5, 1, 1.5)
  sigma <- matrix(c(1, 0.6, 0.4, 0.6, 1.2, 0.7, 0.4, 0.7, 1.5), 3)
  cases <- list(list("cm_1", 1, 2:3, 0.2), list("cm_12", 1:2, 3, c(0.2, 1.4)),
                list("cm_13", c(1, 3), 2, c(0.2, 2)))
  for (cs in cases) {
    cm <- .cond_moments(mean, sigma, cs[[2]], cs[[3]], cs[[4]])
    expect_equal(cm$mean, stage4_cp(cs[[1]], "mean"), tolerance = 1e-12, label = cs[[1]])
    expect_equal(cm$sigma, stage4_cp(cs[[1]], "sigma"), tolerance = 1e-12, label = cs[[1]])
    expect_identical(cm$sigma, t(cm$sigma))
  }
})

test_that(".cond_moments of statistics with independent increments", {
  # Z_k = W_k / sqrt(t_k), where W has independent increments: given Z_1 = z (t_1 = 1), the mean
  # of Z_k is (z + m_k sqrt(t_k) - m_1) / sqrt(t_k) and the covariances are
  # (min(t_k, t_l) - 1) / sqrt(t_k t_l)
  tt <- c(1, 2, 4)
  m <- c(0.3, 0.6, 1.1)
  sigma <- sqrt(outer(tt, tt, pmin) / outer(tt, tt, pmax))
  cm <- .cond_moments(m, sigma, 1, 2:3, 0.7)
  t2 <- tt[2:3]
  expect_equal(cm$mean, (0.7 + m[2:3] * sqrt(t2) - m[1]) / sqrt(t2))
  expect_equal(cm$sigma, (outer(t2, t2, pmin) - 1) / sqrt(outer(t2, t2)))
})
