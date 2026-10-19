test_that(".acc_integral_deriv is the derivative of .acc_integral", {
  scn <- scn_b()
  acc <- .acc_setup(scn$accrual, 300)
  brk <- .scn_breaks(scn)
  hfun <- function(s) .lr_integrands(scn, s, "os")[, c("ev0", "u", "v"), drop = FALSE]
  tt <- c(2, 8, 15)
  d <- .acc_integral_deriv(acc, tt, hfun, brk, hmax = 1)
  h <- 1e-5
  fd <- (.acc_integral(acc, tt + h, hfun, brk, hmax = 1) -
           .acc_integral(acc, tt - h, hfun, brk, hmax = 1)) / (2 * h)
  expect_equal(d, fd, tolerance = 1e-7)
  # with h = 1 the integral of Abar(t - s) has derivative Abar(t)
  one <- function(s) cbind(one = rep(1, length(s)))
  expect_equal(.acc_integral_deriv(acc, tt, one, brk, hmax = 1)[, 1L], .acc_frac(acc, tt),
               tolerance = 1e-12)
})
