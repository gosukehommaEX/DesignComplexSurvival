#' Products of per-patient functionals at two calendar times when switching
#' starts at a calendar time
#'
#' Sum over the groups, weighted by the allocation probabilities, of
#' E[r(t1) r'(t2)] for the functionals psi (score), phi (observed variance)
#' and dlt (event indicator) of the OS log-rank statistic at the looks
#' \code{l1} and \code{l2} (\code{t1 <= t2}). Given the entry time e, the
#' expectation is equation (5) of the covariance note (W = t1 - e,
#' W' = t2 - e); integrated over e it splits into
#' \itemize{
#'   \item the integral over the follow-up time s < W, which over e gives the
#'     group functions \code{D} and \code{Y} of look \code{l1};
#'   \item r0(W) times the integral from W to W' of the functionals of
#'     \code{l2} against the OS density and survival given the entry, and
#'   \item r0(W) r0'(W') times the survival and the dropout survival at W',
#' }
#' where the last two are integrated over W (Gauss-Legendre rules with 12
#' nodes, pieces ending at the pieces of both looks, at \code{t1} minus the
#' accrual breaks, where the availability \code{o} is 0, and where (W, o) or
#' (W', o) lies on a line of \code{.cal_lines()}) with the accrual density at
#' \code{t1 - W}. The inner integral over s uses rules with 10 nodes on
#' pieces ending at the pieces of \code{l2}, the model breaks and the points
#' where (s, o) lies on a line of \code{.cal_lines()}. For a group without
#' switching the inner integral is a difference of a cumulative integral.
#' The functionals at other points than the nodes of their look are the
#' piecewise polynomials of \code{.pw_poly()}.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param acc Output of \code{.acc_setup()}.
#' @param l1,l2 Outputs of \code{.cal_look()} with \code{l1$t <= l2$t}.
#' @param f1,f2 Functionals of the looks by group: lists (control,
#'   treatment) of lists with \code{rho}, \code{beta} (node values, columns
#'   psi, phi, dlt), \code{pr} and \code{pb} (their \code{.pw_poly()}).
#' @return A 3 x 3 matrix (rows: functionals at \code{t1}, columns: at
#'   \code{t2}).
#' @noRd
.cal_pair <- function(scenario, acc, l1, l2, f1, f2) {
  keys <- c("psi", "phi", "dlt")
  t1 <- l1$t
  dlt_t <- l2$t - t1
  tc <- scenario$switch_start
  a <- acc$alpha
  models <- list(scenario$control, scenario$treatment)
  p <- scenario$prob
  same <- identical(models[[1L]], models[[2L]])
  s <- l1$g$x
  w <- l1$g$w
  out <- matrix(0, 3L, 3L, dimnames = list(keys, keys))
  # rule over W = t1 - e
  # lines o = c0 + c1 s; with o = tc - t1 + W they are met at s = W and s = W + dlt_t
  ln <- rbind(.cal_lines(models[[1L]]), .cal_lines(models[[2L]]))
  nh <- abs(ln[, 2L] - 1) > 1e-12
  c0 <- ln[nh, 1L]
  c1 <- ln[nh, 2L]
  wl <- c((c0 - tc + t1) / (1 - c1), (c0 + c1 * dlt_t - tc + t1) / (1 - c1))
  w0 <- max(t1 - acc$duration, 0)
  gw <- .gl_pieces(w0, t1,
                   .cal_breaks(c(l1$g$lo, l1$g$hi, l2$g$lo - dlt_t, l2$g$hi - dlt_t, t1 - a,
                                 t1 - tc, wl), w0, t1),
                   hmax = l1$hmax, n = 12L)
  sl <- abs(ln[, 2L]) > 1e-12
  wv <- gw$x
  ew <- t1 - wv
  ow <- pmax(tc - ew, 0)
  aw <- gw$w * .pw_value(ew, a[-length(a)], acc$rho)
  if (dlt_t > 0) {
    brk_s <- c(l2$g$lo, l2$g$hi, .scn_breaks(scenario))
    inner <- lapply(seq_along(wv), function(k) {
      g <- .gl_pieces(wv[k], wv[k] + dlt_t,
                      .cal_breaks(c(brk_s, (ow[k] - ln[sl, 1L]) / ln[sl, 2L]), wv[k],
                                  wv[k] + dlt_t),
                      hmax = l2$hmax, n = 10L)
      list(x = g$x, w = g$w)
    })
    len <- vapply(inner, function(z) length(z$x), integer(1))
    iw <- rep.int(seq_along(wv), len)
    xs <- unlist(lapply(inner, `[[`, "x"), use.names = FALSE)
    ws <- unlist(lapply(inner, `[[`, "w"), use.names = FALSE)
  }
  os_w <- NULL
  os_in <- NULL
  for (j in 1:2) {
    m <- models[[j]]
    fa <- f1[[j]]
    fb <- f2[[j]]
    r0a <- .pw_poly_eval(fa$pb, s, integral = TRUE)
    r1a <- fa$rho + r0a
    if (dlt_t > 0) {
      r0b <- .pw_poly_eval(fb$pb, s, integral = TRUE)
      r1b <- .pw_poly_eval(fb$pr, s) + r0b
    } else {
      r0b <- r0a
      r1b <- r1a
    }
    gs <- l1$G[[j]]
    term <- crossprod(r1a * (w * gs * l1$D[[j]]), r1b) +
      crossprod(r0a * (w * gs * l1$kap[[j]] * l1$Y[[j]]), r0b)
    # survival at W' given the entry
    if (j == 1L || !same) {
      os_w <- .cal_os(m, wv + dlt_t, ow)
    }
    gw2 <- exp(-.pw_cumhaz(wv + dlt_t, m$time, m$dropout))
    r0a_w <- .pw_poly_eval(fa$pb, wv, integral = TRUE)
    r0b_w2 <- .pw_poly_eval(fb$pb, wv + dlt_t, integral = TRUE)
    term <- term + crossprod(r0a_w * (aw * os_w$surv * gw2), r0b_w2)
    if (dlt_t > 0) {
      if (.has_switch(m)) {
        if (j == 1L || !same) {
          os_in <- .os_avail(m, xs, ow[iw])
        }
        g_in <- exp(-.pw_cumhaz(xs, m$time, m$dropout))
        k_in <- .pw_value(pmax(xs, 0), m$time, m$dropout)
        r0b_in <- .pw_poly_eval(fb$pb, xs, integral = TRUE)
        r1b_in <- .pw_poly_eval(fb$pr, xs) + r0b_in
        hv <- .sum_by(ws * g_in * (r1b_in * os_in$dens + r0b_in * k_in * os_in$surv),
                      iw, length(wv))
      } else {
        # no dependence on the entry: difference of a cumulative integral over the rule of l2
        s2 <- l2$g$x
        ep2 <- .ep_eval(m, s2, "os")
        g2 <- l2$G[[j]]
        hf <- g2 * ((fb$rho + .pw_poly_eval(fb$pb, s2, integral = TRUE)) * ep2$dens +
                      .pw_poly_eval(fb$pb, s2, integral = TRUE) * l2$kap[[j]] * ep2$surv)
        ph <- .pw_poly(l2$g, hf, 12L)
        hv <- .pw_poly_eval(ph, wv + dlt_t, integral = TRUE) -
          .pw_poly_eval(ph, wv, integral = TRUE)
      }
      term <- term + crossprod(r0a_w * aw, hv)
    }
    out <- out + p[[j]] * term
  }
  out
}
