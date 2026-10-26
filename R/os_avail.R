#' Overall survival when switching becomes possible during the follow-up
#'
#' OS survival function and density at the follow-up time \code{s} of a
#' patient who can switch only from the follow-up time \code{o} on (switching
#' from a calendar time: \code{o = max(switch_start - entry, 0)}). Before
#' \code{o} the patient does not switch at progression, and a patient who has
#' progressed dies with \code{death_prog} only; from \code{o} on, a patient
#' who progressed before \code{o} and is alive can switch with
#' \code{switch_rate}, and progressions after \code{o} follow the model. With
#' u the time of progression, the post-progression part is
#' \code{.post_trans()} from the clock value at progression (u not before o)
#' or from the clock value at \code{o} times the survival with
#' \code{death_prog} up to \code{o} (u before o).
#'
#' The integral over u uses Gauss-Legendre rules with 10 nodes on pieces that
#' end at the break points of the integrand: \code{time}, u = o, the points
#' where the post-progression clock at \code{s} (and, with the reset clock, at
#' \code{o}) crosses a break point of \code{time_post}, the further points of
#' \code{.os_conv()} for the switch clock and the AFT effect and, with the AFT
#' effect and the reset clock, the points where the clock value at \code{o}
#' meets a break of the integral over the switch time in
#' \code{.post_general()}. Pieces are at most 2 / lambda long, lambda being
#' the largest intensity. The pairs are handled in chunks.
#'
#' @param model An \code{ms_model} object without switching before
#'   progression.
#' @param s,o Numeric vectors of equal length: follow-up times and the
#'   follow-up times from which switching is possible.
#' @param chunk Number of pairs handled at once.
#' @return A list with \code{surv} and \code{dens}.
#' @noRd
.os_avail <- function(model, s, o, chunk = 4000L) {
  ns <- length(s)
  tm <- model$time
  p11 <- exp(-.pw_cumhaz(s, tm, model$prog + model$death))
  surv <- p11
  dens <- .pw_value(pmax(s, 0), tm, model$death) * p11
  if (ns == 0L || all(model$prog == 0)) {
    return(list(surv = surv, dens = dens))
  }
  tp <- model$time_post
  lam <- max(c(model$prog + model$death, model$switch_rate + model$death_prog,
               model$death_switch), 1e-12)
  h <- 2 / lam
  forward <- model$clock == "forward"
  sw <- identical(model$switch_clock, "switch")
  phi <- if (is.null(model$switch_af) || model$switch_af == 1) NULL else model$switch_af
  rule <- .gl_rule(10L)
  for (a in seq(1L, ns, by = chunk)) {
    i <- a:min(ns, a + chunk - 1L)
    si <- s[i]
    oi <- o[i]
    m <- length(i)
    rowv <- function(v) matrix(v, m, length(v), byrow = TRUE)
    cols <- list(matrix(0, m, 1L), rowv(tm), oi)
    if (forward) {
      cols <- c(cols, list(rowv(tp)))
      if (sw) {
        cols <- c(cols, list(outer(si, tp, function(x, y) x - y)))
      } else if (!is.null(phi)) {
        cols <- c(cols, list(outer(si, tp, function(x, y) (y - phi * x) / (1 - phi))))
      }
    } else {
      cols <- c(cols, list(outer(si, tp, "-"), outer(oi, tp, "-")))
      if (sw) {
        cols <- c(cols, list(outer(si, as.vector(outer(tp, tp, "+")), "-")))
      } else if (!is.null(phi)) {
        ex <- as.vector(outer(tp, tp, function(x, y) (x - (1 - phi) * y) / phi))
        cols <- c(cols, list(outer(si, ex, "-"), outer(oi * (1 - phi) + phi * si, tp, "-")))
      }
    }
    smax <- max(si)
    if (smax > h) {
      cols <- c(cols, list(rowv(seq(h, smax, by = h))))
    }
    br <- do.call(cbind, c(cols, list(si)))
    br <- pmin(pmax(br, 0), si)
    nb <- ncol(br)
    v <- as.vector(base::t(br))
    v <- v[order(rep(seq_len(m), each = nb), v)]
    br <- matrix(v, m, nb, byrow = TRUE)
    lo <- br[, -nb, drop = FALSE]
    hi <- br[, -1L, drop = FALSE]
    keep <- which(hi > lo)
    if (length(keep) == 0L) {
      next
    }
    row_k <- (keep - 1L) %% m + 1L
    half <- (hi[keep] - lo[keep]) / 2
    mid <- (hi[keep] + lo[keep]) / 2
    u <- as.vector(outer(rule$x, half)) + rep(mid, each = 10L)
    w <- as.vector(outer(rule$w, half))
    pr <- rep(row_k, each = 10L)
    sv <- si[pr]
    ov <- oi[pr]
    kap <- exp(-.pw_cumhaz(u, tm, model$prog + model$death)) * .pw_value(u, tm, model$prog)
    after <- u >= ov
    t_u <- if (forward) u else numeric(length(u))
    t_s <- if (forward) sv else sv - u
    gs <- numeric(length(u))
    gd <- numeric(length(u))
    if (any(after)) {
      tr <- .post_trans(model, t_u[after], t_s[after])
      pp <- .pw_value(u[after], tm, model$switch_prob)
      gs[after] <- (1 - pp) * (tr$r22 + tr$r23) + pp * tr$r33
      gd[after] <- (1 - pp) * (tr$r22 * tr$c_to + tr$f23) + pp * tr$r33 * tr$d_to
    }
    b <- !after
    if (any(b)) {
      t_a <- pmin(if (forward) ov[b] else ov[b] - u[b], t_s[b])
      pre <- exp(-(.pw_cumhaz(t_a, tp, model$death_prog) -
                     .pw_cumhaz(t_u[b], tp, model$death_prog)))
      tr <- .post_trans(model, t_a, t_s[b])
      gs[b] <- pre * (tr$r22 + tr$r23)
      gd[b] <- pre * (tr$r22 * tr$c_to + tr$f23)
    }
    sm <- .sum_by(cbind(w * kap * gs, w * kap * gd), pr, m)
    surv[i] <- surv[i] + sm[, 1L]
    dens[i] <- dens[i] + sm[, 2L]
  }
  list(surv = surv, dens = dens)
}
