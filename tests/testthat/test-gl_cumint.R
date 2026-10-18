test_that(".gl_cumint gives the cumulative integral at the nodes", {
  q <- .gl_cumint(8L)
  x <- .gl_rule(8L)$x
  # integral of 4 x^3 + 1 from -1 to x is x^4 + x
  expect_equal(as.vector(q %*% (4 * x^3 + 1)), x^4 + x, tolerance = 1e-12)
  # degree 7 is exact
  expect_equal(as.vector(q %*% (8 * x^7)), x^8 - 1, tolerance = 1e-12)
})
