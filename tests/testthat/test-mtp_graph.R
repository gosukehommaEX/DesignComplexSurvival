# With independent statistics (identity covariance) the rectangle
# probabilities are products of normal probabilities, and the rejection
# probabilities of the graphical procedure have closed forms over the
# outcomes of Algorithm 1 of Maurer and Bretz (2013).
graph_box_indep <- function(mean) {
  function(idx, lo, hi) {
    if (any(lo >= hi)) {
      return(0)
    }
    prod(stats::pnorm(hi, mean[idx]) - stats::pnorm(lo, mean[idx]))
  }
}

test_that(".mtp_graph gives the closed forms for one look of each endpoint", {
  m <- c(1.8, 1.5)
  up <- function(x, j) stats::pnorm(x, m[j], lower.tail = FALSE)
  a1 <- 2.2
  b1 <- 1.9
  a2 <- 2.4
  b2 <- 2.0
  # both endpoints at analysis 1
  r <- .mtp_graph(list(a1, a2), list(b1, b2), -Inf, list(1, 1), graph_box_indep(m))
  both <- up(a1, 1) * up(b2, 2) + (up(b1, 1) - up(a1, 1)) * up(a2, 2)
  expect_equal(r$both, both, tolerance = 1e-12)
  expect_equal(r$p1, up(a1, 1) + (up(b1, 1) - up(a1, 1)) * up(a2, 2), tolerance = 1e-12)
  expect_equal(r$p2, up(a2, 2) + up(a1, 1) * (up(b2, 2) - up(a2, 2)), tolerance = 1e-12)
  expect_identical(r$pf, 0)
  # the second endpoint at analysis 1 and the first at analysis 2
  r <- .mtp_graph(list(a1, a2), list(b1, b2), -Inf, list(2, 1), graph_box_indep(m))
  expect_equal(r$p2, up(a2, 2), tolerance = 1e-12)
  expect_equal(r$p1, up(a2, 2) * up(b1, 1) + (1 - up(a2, 2)) * up(a1, 1), tolerance = 1e-12)
  expect_equal(r$both, up(a2, 2) * up(b1, 1), tolerance = 1e-12)
})

test_that(".mtp_graph applies the futility bound before and after a rejection", {
  # first endpoint at analyses 1 and 2 with a futility bound at look 1,
  # second endpoint at analysis 2; coordinates X1, X2 (first) and Y (second)
  m <- c(1.2, 1.9, 1.6)
  px <- function(lo, hi, j) stats::pnorm(hi, m[j]) - stats::pnorm(lo, m[j])
  r <- .mtp_graph(list(c(2.5, 2.0), 2.1), list(c(2.2, 1.8), 1.9), c(0.3, 1), list(1:2, 2),
                  graph_box_indep(m))
  q1 <- px(2.5, Inf, 1)
  c1 <- px(0.3, 2.5, 1)
  y_b <- px(1.9, Inf, 3)
  y_a <- px(2.1, Inf, 3)
  expect_equal(r$p1, c(q1, c1 * (px(2, Inf, 2) + px(1.8, 2, 2) * y_a)), tolerance = 1e-12)
  expect_equal(r$p2, q1 * y_b + c1 * (px(2, Inf, 2) * y_b + px(-Inf, 2, 2) * y_a),
               tolerance = 1e-12)
  expect_equal(r$both, q1 * y_b + c1 * (px(2, Inf, 2) * y_b + px(1.8, 2, 2) * y_a),
               tolerance = 1e-12)
  expect_equal(r$pf, c(px(-Inf, 0.3, 1), 0), tolerance = 1e-12)
  # futility bound (1.5) above the bound after recycling (1.0) at look 1 of
  # the first endpoint, with the second endpoint at analysis 1
  m <- c(1.0, 1.8, 2.0)
  px <- function(lo, hi, j) stats::pnorm(hi, m[j]) - stats::pnorm(lo, m[j])
  r <- .mtp_graph(list(c(2.6, 2.0), 2.2), list(c(1.0, 1.7), 1.9), c(1.5, 2), list(1:2, 1),
                  graph_box_indep(m))
  expect_equal(r$p1, c(px(2.6, Inf, 1) + px(1.0, 2.6, 1) * px(2.2, Inf, 3),
                       px(1.5, 2.6, 1) * px(-Inf, 2.2, 3) * px(2.0, Inf, 2)), tolerance = 1e-12)
  expect_equal(r$p2, px(2.6, Inf, 1) * px(1.9, Inf, 3) + px(-Inf, 2.6, 1) * px(2.2, Inf, 3),
               tolerance = 1e-12)
  expect_equal(r$both, px(2.6, Inf, 1) * px(1.9, Inf, 3) + px(1.0, 2.6, 1) * px(2.2, Inf, 3),
               tolerance = 1e-12)
  expect_equal(r$pf, c(px(-Inf, 1.0, 1) * px(2.2, Inf, 3) + px(-Inf, 1.5, 1) * px(-Inf, 2.2, 3),
                       0), tolerance = 1e-12)
})

test_that(".mtp_graph uses five rectangles at an analysis with both endpoints", {
  # three common analyses without futility bounds: 4 rectangles at analysis
  # 1, then 1 per analysis after one rejection and 4 + 1 + 1 + 4 from the
  # outcome with no rejection (18 in total); none after an infinite bound
  n <- 0L
  count_box <- function(idx, lo, hi) {
    n <<- n + 1L
    0.5
  }
  .mtp_graph(list(c(3, 2.5, 2), c(3, 2.5, 2)), list(c(2.8, 2.2, 1.9), c(2.7, 2.3, 1.8)),
             rep(-Inf, 3), list(1:3, 1:3), count_box)
  expect_identical(n, 18L)
  n <- 0L
  .mtp_graph(list(c(Inf, 2.5, 2), c(3, 2.5, 2)), list(c(Inf, 2.2, 1.9), c(2.7, 2.3, 1.8)),
             rep(-Inf, 3), list(1:3, 1:3), count_box)
  expect_identical(n, 16L)
})
