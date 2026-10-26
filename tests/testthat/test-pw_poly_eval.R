test_that(".pw_poly_eval evaluates the polynomials and their integrals", {
  g <- .gl_pieces(0, 10, c(2.5, 7), hmax = 3, n = 6L)
  f <- function(x) cbind(x^5 - 3 * x^2 + 1, 2 - x)
  fi <- function(x) cbind(x^6 / 6 - x^3 + x, 2 * x - x^2 / 2)
  pp <- .pw_poly(g, f(g$x), 6L)
  x <- c(0, 0.3, 2.5, 4.1, 7, 9.99, 10)
  expect_equal(.pw_poly_eval(pp, x), f(x), tolerance = 1e-12)
  expect_equal(.pw_poly_eval(pp, x, integral = TRUE), fi(x), tolerance = 1e-12)
  # points outside the pieces are moved to the nearest end
  expect_equal(.pw_poly_eval(pp, c(-1, 12), integral = TRUE), fi(c(0, 10)), tolerance = 1e-12)
  expect_equal(nrow(.pw_poly_eval(pp, numeric(0))), 0L)
})
