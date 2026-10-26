test_that(".pw_poly reproduces polynomials of degree below the number of nodes", {
  g <- .gl_pieces(0, 10, c(2.5, 7), hmax = 3, n = 6L)
  f <- function(x) cbind(x^5 - 3 * x^2 + 1, 2 - x)
  pp <- .pw_poly(g, f(g$x), 6L)
  expect_equal(dim(pp$coef), c(6L, length(g$lo), 2L))
  # integrals from 0 to the start of each piece
  fi <- function(x) cbind(x^6 / 6 - x^3 + x, 2 * x - x^2 / 2)
  expect_equal(pp$cum_lo, fi(g$lo), tolerance = 1e-12)
})
