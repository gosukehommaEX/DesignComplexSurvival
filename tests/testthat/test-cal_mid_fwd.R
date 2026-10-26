# Reference: the integral over the entry time of .os_avail() (Gauss-Legendre rules with pieces
# at most 0.25 long ending at the accrual breaks and the lines of .cal_lines()), the
# computation that .cal_group() uses for the reset clock

mid_by_entry <- function(model, acc, s, lo, hi, tc) {
  ln <- .cal_lines(model)
  a <- acc$alpha
  out <- matrix(0, length(s), 2L)
  for (k in seq_along(s)) {
    brk <- .cal_breaks(c(a, tc - (ln[, 1L] + ln[, 2L] * s[k])), lo[k], hi[k])
    g <- .gl_pieces(lo[k], hi[k], brk, hmax = 0.25, n = 10L)
    w <- g$w * .pw_value(g$x, a[-length(a)], acc$rho)
    av <- .os_avail(model, rep(s[k], length(g$x)), tc - g$x)
    out[k, ] <- c(sum(w * av$surv), sum(w * av$dens))
  }
  out
}

test_that(".cal_mid_fwd agrees with the integral over the entry time of .os_avail", {
  acc <- .acc_setup(scn_cal()$accrual, 500)
  for (cs in list(list(13.7, 33.1), list(4.2, 30))) {
    tc <- cs[[1L]]
    t <- cs[[2L]]
    s <- seq(0.5, t - 0.5, length.out = 9)
    lo <- pmax(tc - s, 0)
    hi <- pmin(t - s, acc$duration, tc)
    k <- hi > lo
    for (nm in c("forward_same", "forward_switch", "forward_aft0.6", "forward_aft1.5")) {
      m <- model_post_cal(nm)
      got <- .cal_mid_fwd(m, acc, s[k], lo[k], hi[k], tc)
      ref <- mid_by_entry(m, acc, s[k], lo[k], hi[k], tc)
      expect_lt(max(abs(got$surv - ref[, 1L])), 1e-13)
      expect_lt(max(abs(got$dens - ref[, 2L])), 1e-14)
    }
  }
})

test_that(".cal_mid_fwd handles the follow-up times in chunks and empty input", {
  acc <- .acc_setup(scn_cal()$accrual, 500)
  m <- model_post_cal("forward_switch")
  s <- seq(1, 20, length.out = 7)
  lo <- pmax(15 - s, 0)
  hi <- pmin(36 - s, 18, 15)
  expect_equal(.cal_mid_fwd(m, acc, s, lo, hi, 15, chunk = 3L), .cal_mid_fwd(m, acc, s, lo, hi, 15),
               tolerance = 1e-14)
  expect_equal(.cal_mid_fwd(m, acc, numeric(0), numeric(0), numeric(0), 15),
               list(surv = numeric(0), dens = numeric(0)))
})
