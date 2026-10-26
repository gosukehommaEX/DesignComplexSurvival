#' Follow-up times where a line of the (s, o) plane meets an edge of the
#' range of the availability
#'
#' At a look at calendar time t, the follow-up times in (0, t) at which a
#' line of \code{.cal_lines()} meets an edge of the range of the availability
#' o, or two edges meet, inside that range. The edges are o = s
#' (entry time \code{tc - s}), o = tc - t + s (entry time \code{t - s}),
#' o = tc - a for the accrual breaks a, and o = 0; the range of o at the
#' follow-up time s is from max(tc - t + s, tc - duration, 0) to min(s, tc).
#' The group functions of \code{.cal_group()} are not smooth in s at these
#' times. Points where two lines meet inside the range are not included:
#' there the group functions are smoother, and with them the rules of
#' \code{.cal_look()} would get many short pieces.
#'
#' @param lines Matrix with columns \code{c0} and \code{c1} (lines
#'   o = c0 + c1 s).
#' @param t Calendar time of the look.
#' @param tc Calendar time from which switching is possible.
#' @param acc Output of \code{.acc_setup()}.
#' @return Numeric vector of follow-up times.
#' @noRd
.cal_cross <- function(lines, t, tc, acc) {
  a <- acc$alpha
  ed <- rbind(cbind(c(tc - a, 0), 0), c(0, 1), c(tc - t, 1))
  al <- rbind(unname(lines), ed)
  d1 <- outer(al[, 2L], al[, 2L], "-")
  d0 <- outer(al[, 1L], al[, 1L], function(x, y) y - x)
  # pairs i < j with j an edge
  ok <- upper.tri(d1) & abs(d1) > 1e-12 & col(d1) > nrow(lines)
  sx <- d0[ok] / d1[ok]
  i <- row(d1)[ok]
  ox <- al[i, 1L] + al[i, 2L] * sx
  lo <- pmax(tc - t + sx, tc - acc$duration, 0)
  hi <- pmin(sx, tc)
  sx[sx > 0 & sx < t & ox >= lo - 1e-9 & ox <= hi + 1e-9]
}
