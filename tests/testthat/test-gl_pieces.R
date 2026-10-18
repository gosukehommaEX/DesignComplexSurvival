test_that(".gl_pieces splits at breaks and by the maximum length", {
  g <- .gl_pieces(0, 10, breaks = c(3, 12), hmax = 4, n = 5L)
  expect_equal(g$lo, c(0, 3, 6.5))
  expect_equal(g$hi[length(g$hi)], 10)
  expect_equal(sum(g$w), 10, tolerance = 1e-12)
  expect_equal(sum(g$w * exp(-0.3 * g$x)), (1 - exp(-3)) / 0.3, tolerance = 1e-10)
  expect_length(.gl_pieces(1, 1)$x, 0)
})
