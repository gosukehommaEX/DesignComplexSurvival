#' Integral over the entry time of the OS given the availability (forward
#' clock)
#'
#' For a model with the forward clock, the part of \code{.cal_group()} over
#' the entry times e from \code{lo} to \code{hi}, where the availability
#' o = tc - e lies between 0 and s, computed without the integral over e of
#' \code{.os_avail()}. With u the time of progression, kap(u) its density
#' before progression, a(e) the accrual density, C(x) the cumulative hazard
#' of \code{death_prog} at the clock value x and R(x, s) the probability of
#' being alive in the progressed or switched state at s after being alive
#' and progressed at x without switching (\code{.post_trans()}), the
#' integral of a(e) S(s | tc - e) over e is p11(s) A plus the integral over
#' u from 0 to s of kap(u) [B(u) G(u, s) + exp(C(u)) H(u)], where A is the
#' accrued fraction between \code{lo} and \code{hi}, B(u) the part of it
#' with o <= u (progression from the availability on: G is the
#' post-progression survival of the model, as in \code{.os_avail()}), and
#' H(u) the integral over o > u of a(tc - o) exp(-C(o)) R(o, s)
#' (progression before the availability: no switch at progression, survival
#' with \code{death_prog} up to o, from where the patient can switch). The
#' factor exp(C(u) - C(o)) of \code{.os_avail()} splits because the clock
#' does not restart at progression. The density is the same with the
#' densities of the death at s. The integrals over u and o use one
#' Gauss-Legendre rule with 10 nodes per piece, pieces ending at 0,
#' \code{time}, \code{time_post}, tc minus the accrual breaks, tc - hi,
#' tc - lo, s, the points of \code{.os_conv()} for the switch clock and the
#' AFT effect and at most 2 / lambda long (lambda the largest intensity);
#' H at the nodes is the cumulative integral of its integrand at the nodes
#' (\code{.gl_cumint()}).
#'
#' @param model An \code{ms_model} object with \code{clock = "forward"}.
#' @param acc Output of \code{.acc_setup()}.
#' @param s Follow-up times.
#' @param lo,hi Entry times (vectors as \code{s}), with
#'   \code{max(tc - s, 0) <= lo < hi <= tc}.
#' @param tc Calendar time from which switching is possible.
#' @param chunk Number of follow-up times handled at once.
#' @return A list with \code{surv} and \code{dens}, the integrals of the OS
#'   survival function and density times the accrual density.
#' @noRd
.cal_mid_fwd <- function(model, acc, s, lo, hi, tc, chunk = 4000L) {
  ns <- length(s)
  tm <- model$time
  tp <- model$time_post
  a <- acc$alpha
  p11 <- exp(-.pw_cumhaz(s, tm, model$prog + model$death))
  a_e <- .acc_cum(acc, hi) - .acc_cum(acc, lo)
  surv <- p11 * a_e
  dens <- .pw_value(pmax(s, 0), tm, model$death) * p11 * a_e
  if (ns == 0L || all(model$prog == 0)) {
    return(list(surv = surv, dens = dens))
  }
  lam <- max(c(model$prog + model$death, model$switch_rate + model$death_prog,
               model$death_switch), 1e-12)
  h <- 2 / lam
  sw <- identical(model$switch_clock, "switch")
  phi <- if (is.null(model$switch_af) || model$switch_af == 1) NULL else model$switch_af
  rule <- .gl_rule(10L)
  qm <- .gl_cumint(10L)
  for (c0 in seq(1L, ns, by = chunk)) {
    i <- c0:min(ns, c0 + chunk - 1L)
    si <- s[i]
    lo_i <- lo[i]
    hi_i <- hi[i]
    m <- length(i)
    rowv <- function(v) matrix(v, m, length(v), byrow = TRUE)
    cols <- list(matrix(0, m, 1L), rowv(tm), rowv(tp), rowv(tc - a), tc - hi_i, tc - lo_i)
    if (sw) {
      cols <- c(cols, list(outer(si, tp, function(x, y) x - y)))
    } else if (!is.null(phi)) {
      cols <- c(cols, list(outer(si, tp, function(x, y) (y - phi * x) / (1 - phi))))
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
    pl <- br[, -nb, drop = FALSE]
    pu <- br[, -1L, drop = FALSE]
    # pieces in the order of the rows, and within a row in increasing order
    kk <- which(base::t(pu > pl)) - 1L
    if (length(kk) == 0L) {
      next
    }
    rk <- kk %/% (nb - 1L) + 1L
    pc <- cbind(rk, kk %% (nb - 1L) + 1L)
    lo_k <- pl[pc]
    hi_k <- pu[pc]
    np <- length(rk)
    half <- (hi_k - lo_k) / 2
    mid <- (hi_k + lo_k) / 2
    u <- as.vector(outer(rule$x, half)) + rep(mid, each = 10L)
    w <- as.vector(outer(rule$w, half))
    pr <- rep(rk, each = 10L)
    kap <- exp(-.pw_cumhaz(u, tm, model$prog + model$death)) * .pw_value(u, tm, model$prog)
    tr <- .post_trans(model, u, si[pr])
    pp <- .pw_value(u, tm, model$switch_prob)
    gs <- (1 - pp) * (tr$r22 + tr$r23) + pp * tr$r33
    gd <- (1 - pp) * (tr$r22 * tr$c_to + tr$f23) + pp * tr$r33 * tr$d_to
    # accrued fraction with o = tc - e <= u
    lo_n <- lo_i[pr]
    hi_n <- hi_i[pr]
    b <- .acc_cum(acc, hi_n) - .acc_cum(acc, pmin(pmax(tc - u, lo_n), hi_n))
    # integrand over o at the nodes (o = u): a(tc - o) on (tc - hi, tc - lo)
    e_o <- tc - u
    inside <- e_o > lo_n & e_o < hi_n
    ao <- numeric(length(u))
    ao[inside] <- .pw_value(e_o[inside], a[-length(a)], acc$rho)
    cdp <- .pw_cumhaz(u, tp, model$death_prog)
    ecd <- exp(-cdp)
    es <- exp(cdp)
    if (any(!is.finite(es))) {
      stop("Numerical overflow in the integral over the entry time; use a larger time unit.",
           call. = FALSE)
    }
    last <- c(rk[-1L] != rk[-np], TRUE)
    # integral from each node to s: the rest of its piece and the later pieces of the row
    upper <- function(hv) {
      hm <- matrix(hv, 10L, np)
      ip <- half * colSums(rule$w * hm)
      cum <- (qm %*% hm) * rep(half, each = 10L)
      cs <- cumsum(ip)
      row_end <- numeric(m)
      row_end[rk[last]] <- cs[last]
      as.vector(rep(row_end[rk] - cs + ip, each = 10L) - cum)
    }
    hs <- upper(ao * ecd * (tr$r22 + tr$r23))
    hd <- upper(ao * ecd * (tr$r22 * tr$c_to + tr$f23))
    sm <- .sum_by(cbind(w * kap * (b * gs + es * hs), w * kap * (b * gd + es * hd)), pr, m)
    surv[i] <- surv[i] + sm[, 1L]
    dens[i] <- dens[i] + sm[, 2L]
  }
  list(surv = surv, dens = dens)
}
