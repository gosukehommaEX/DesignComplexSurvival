test_that(".gl_rule integrates polynomials of degree 2n - 1 exactly", {
  r <- .gl_rule(6L)
  expect_equal(sum(r$w), 2, tolerance = 1e-14)
  expect_equal(sum(r$w * r$x^10), 2 / 11, tolerance = 1e-13)
  expect_equal(sum(r$w * r$x^11), 0, tolerance = 1e-13)
  expect_true(all(diff(r$x) > 0))
})
