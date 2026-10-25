test_that(".recycled_upper gives the bounds of gs_bounds at the level after recycling", {
  tf <- c(0.4, 0.7, 1)
  corr <- sqrt(outer(tf, tf, pmin) / outer(tf, tf, pmax))
  a0 <- gs_bounds(tf, alpha = 0.0125, sf_eff = "OF")$upper
  b <- .recycled_upper(tf, 0.025, "OF", NULL, corr, a0, NULL, "pfs")
  expect_identical(b, gs_bounds(tf, alpha = 0.025, sf_eff = "OF")$upper)
  expect_true(all(b < a0))
  # binding futility bounds of the initial level
  lo <- c(0.2, 0.9, 2)
  a1 <- gs_bounds(tf, alpha = 0.0125, sf_eff = "Pocock", futility = "fixed",
                  fut_bound = lo[1:2], binding = TRUE)$upper
  b1 <- .recycled_upper(tf, 0.02, "Pocock", NULL, corr, a1, lo, "pfs")
  expect_identical(b1, gs_bounds(tf, alpha = 0.02, sf_eff = "Pocock", futility = "fixed",
                                 fut_bound = lo[1:2], binding = TRUE)$upper)
  nb <- gs_bounds(tf, alpha = 0.02, sf_eff = "Pocock")$upper
  expect_identical(b1[1], nb[1])
  expect_true(all(b1[2:3] < nb[2:3]))
  # user-defined spending
  pu <- c(0.2, 0.5, 1)
  au <- gs_bounds(tf, alpha = 0.01, sf_eff = "user", param_eff = pu)$upper
  bu <- .recycled_upper(tf, 0.025, "user", pu, corr, au, NULL, "os")
  expect_identical(bu, gs_bounds(tf, alpha = 0.025, sf_eff = "user", param_eff = pu)$upper)
})

test_that(".recycled_upper stops when the bounds after recycling are not valid", {
  tf <- c(0.5, 1)
  corr <- sqrt(outer(tf, tf, pmin) / outer(tf, tf, pmax))
  a0 <- gs_bounds(tf, alpha = 0.0125)$upper
  expect_error(.recycled_upper(tf, 0.025, "OF", NULL, corr, a0 - c(0, 0.5), NULL, "os"),
               "efficacy bound of OS at its look 2 at the level 0.025 after recycling exceeds")
  # binding futility bound at look 1 (3.8) above the bound after recycling
  # (3.71), so that no trial continues to look 2 at the level after recycling
  t3 <- c(1, 2, 3) / 3
  c3 <- sqrt(outer(t3, t3, pmin) / outer(t3, t3, pmax))
  a3 <- gs_bounds(t3, alpha = 0.0125)$upper
  expect_error(.recycled_upper(t3, 0.025, "OF", NULL, c3, a3, c(3.8, 0, 2), "pfs"),
               "Cannot compute the efficacy bounds of PFS at the level 0.025 after recycling")
})
