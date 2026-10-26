#' Functions of both groups at one calendar time when switching starts at a
#' calendar time
#'
#' Gauss-Legendre rule with 12 nodes per piece over the follow-up time from 0
#' to \code{t}, with the break points of \code{.scn_breaks()}, \code{t} minus
#' the accrual breaks and \code{switch_start}, \code{switch_start} minus the
#' accrual breaks, \code{switch_start}, and the follow-up times of
#' \code{.cal_cross()} where a line of \code{.cal_lines()} of the groups meets
#' an edge of the range of the availability (near duplicates removed by
#' \code{.cal_breaks()}); pieces at most \code{.scn_hmax()} long. At the
#' nodes, the group functions of \code{.cal_group()}, the dropout survival
#' and hazard, the proportion of the treatment group in the risk set
#' \code{q1}, the pooled hazard \code{lbar}, and the expected events,
#' dropouts, score mean and information per randomized patient. Identical
#' models are evaluated once.
#'
#' @param scenario An \code{ms_scenario} object with \code{switch_start}
#'   positive.
#' @param acc Output of \code{.acc_setup()}.
#' @param t Calendar time (positive).
#' @return A list with \code{t}, the rule \code{g}, \code{hmax}, the lists by
#'   group (control, treatment) \code{Y}, \code{D}, \code{G} and \code{kap},
#'   the node values \code{q1} and \code{lbar}, the events \code{ev} and
#'   dropouts \code{do} by group per randomized patient (including the
#'   allocation probability), the score mean \code{mu} and the information
#'   \code{v}.
#' @noRd
.cal_look <- function(scenario, acc, t) {
  tc <- scenario$switch_start
  a <- acc$alpha
  hmax <- .scn_hmax(scenario, max(t, 1))
  ln <- rbind(.cal_lines(scenario$control), .cal_lines(scenario$treatment))
  brk <- .cal_breaks(c(.scn_breaks(scenario), t - a, t - tc, tc - a, tc,
                       .cal_cross(ln, t, tc, acc)), 0, t)
  g <- .gl_pieces(0, t, brk, hmax = hmax, n = 12L)
  s <- g$x
  models <- list(scenario$control, scenario$treatment)
  p <- scenario$prob
  grp <- vector("list", 2L)
  grp[[1L]] <- .cal_group(models[[1L]], acc, t, s, tc, hmax)
  grp[[2L]] <- if (identical(models[[1L]], models[[2L]])) {
    grp[[1L]]
  } else {
    .cal_group(models[[2L]], acc, t, s, tc, hmax)
  }
  gg <- lapply(models, function(m) exp(-.pw_cumhaz(s, m$time, m$dropout)))
  kap <- lapply(models, function(m) .pw_value(pmax(s, 0), m$time, m$dropout))
  y <- lapply(1:2, function(j) p[[j]] * gg[[j]] * grp[[j]]$Y)
  d <- lapply(1:2, function(j) p[[j]] * gg[[j]] * grp[[j]]$D)
  ys <- y[[1L]] + y[[2L]]
  ok <- ys > 0
  ns <- length(s)
  q1 <- rep(p[[2L]], ns)
  lbar <- numeric(ns)
  u <- numeric(ns)
  q1[ok] <- y[[2L]][ok] / ys[ok]
  lbar[ok] <- (d[[1L]][ok] + d[[2L]][ok]) / ys[ok]
  u[ok] <- (y[[1L]][ok] * d[[2L]][ok] - y[[2L]][ok] * d[[1L]][ok]) / ys[ok]
  v <- q1 * (1 - q1) * (d[[1L]] + d[[2L]])
  list(t = t, g = g, hmax = hmax, Y = lapply(grp, `[[`, "Y"), D = lapply(grp, `[[`, "D"),
       G = gg, kap = kap, q1 = q1, lbar = lbar,
       ev = vapply(d, function(z) sum(g$w * z), numeric(1)),
       do = vapply(1:2, function(j) sum(g$w * y[[j]] * kap[[j]]), numeric(1)),
       mu = sum(g$w * u), v = sum(g$w * v))
}
