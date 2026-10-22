# Expected values: inst/validation/python/results/pmvn_genz.csv, computed
# independently by expected_genz.py. D1 (four coordinates): recursive
# conditioning with adaptive quadrature (ref_gsd.mvn_box()). F6 (six
# coordinates, three looks of each endpoint): quasi-Monte Carlo with the
# separation of variables of Genz and 24 scramblings of 2^20 Sobol points
# (three standard errors 3e-9). With 128 grid points, the Miwa algorithm is
# off by 9.6e-5 (D1) and 4.5e-4 (F6).

genz_d1 <- list(
  m = c(2.141636271, 2.748018643, 2.034474618, 2.609699421),
  S = matrix(c(1.0669695172, 0.8174879172, 0.3779545106, 0.2913956653,
               0.8174879172, 1.0538823953, 0.4843409699, 0.3728319907,
               0.3779545106, 0.4843409699, 1.0880711119, 0.8376117154,
               0.2913956653, 0.3728319907, 0.8376117154, 1.0767925675), 4),
  lo = c(-Inf, 1.980965029, -Inf, 1.981224804),
  hi = c(2.668630143, Inf, 2.664375831, Inf),
  exact = 0.22710652488131736)

genz_f6 <- list(
  m = c(2.3266834358, 2.91543213425, 3.21678539962, 1.66187457073, 2.19259585036,
        2.60969942072),
  S = matrix(c(1.0629926776, 0.832669469667, 0.748564986909, 0.502349318637,
               0.378467094341, 0.316332203304,
               0.832669469667, 1.04791190754, 0.934245939148, 0.395728002013,
               0.473693925547, 0.395540252994,
               0.748564986909, 0.934245939148, 1.02402167357, 0.355455009638,
               0.42309965197, 0.432352726809,
               0.502349318637, 0.395728002013, 0.355455009638, 1.08808430447,
               0.818750399697, 0.683601996344,
               0.378467094341, 0.473693925547, 0.42309965197, 0.818750399697,
               1.08404799311, 0.903700658043,
               0.316332203304, 0.395540252994, 0.432352726809, 0.683601996344,
               0.903700658043, 1.07679256746), 6),
  lo = c(-Inf, -Inf, 2.02976766455, -Inf, -Inf, 2.00053906808),
  hi = c(2.93857575243, 2.25655974097, Inf, 3.35686935609, 2.44454228117, Inf),
  exact = 0.05845214941591478)

genz_value <- function(x) .pmvn_genz(x$lo, x$hi, x$m, x$S)

test_that(".pmvn_genz matches the exact rectangle probabilities", {
  for (x in list(genz_d1, genz_f6)) {
    expect_lt(abs(with_abseps(genz_value(x), 1e-6) - x$exact), 2e-6)
    # default tolerance 1e-5
    expect_lt(abs(genz_value(x) - x$exact), 2e-5)
  }
})

test_that(".pmvn_genz is reproducible and keeps the random number state", {
  set.seed(123)
  r1 <- stats::runif(2)
  set.seed(123)
  v1 <- genz_value(genz_f6)
  r2 <- stats::runif(2)
  v2 <- genz_value(genz_f6)
  expect_identical(v1, v2)
  expect_identical(r1, r2)
})

test_that("gs_prob_mtp uses .pmvn_genz for rectangles over both endpoints", {
  # first endpoint at analyses 1 and 2, second at analyses 2 and 3 (fixed
  # sequence, overall): the second hypothesis is rejected at its look 2 after
  # the first is rejected at look 1 (three coordinates, exact
  # 0.12708345860555498 by ref_gsd.mvn_box()) or at look 2 (the four
  # coordinates of genz_d1)
  u1 <- c(genz_d1$hi[1], genz_d1$lo[2])
  u2 <- c(genz_d1$hi[3], genz_d1$lo[4])
  r <- with_abseps(gs_prob_mtp(list(u1, u2), mean = genz_d1$m, sigma = genz_d1$S,
                               analysis = list(1:2, 2:3)), 1e-6)
  expect_lt(abs(r$looks$prob_upper[4] - (0.12708345860555498 + genz_d1$exact)), 2e-6)
})
