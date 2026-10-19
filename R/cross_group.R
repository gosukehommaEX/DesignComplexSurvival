#' Products of the PFS and OS functionals in one group
#'
#' E[r_P(t_p) r_O(t_o)] for a patient of group \code{grp}, where r_P is a
#' functional of PFS at calendar time \code{t_p} and r_O a functional of OS at
#' \code{t_o} (score, observed variance or event indicator, see
#' \code{.fun_vals()}). Writing r = int rho dN + int beta Y ds, the product
#' expectation is C1 + C2 + C3 + C4 of the covariance note (Section 4.3), with
#' the weight omega(s, u) = G(max(s, u)) Abar(min(t_p - s, t_o - u)) for the
#' dropout survival G and the enrolled fraction Abar. The joint law of PFS and
#' OS consists of death without progression (density S_P(s) b(s) on the
#' diagonal) and progression at v followed by death at u (density
#' S_P(v) a(v) f_post(u | v)). The integrals over the OS time u use an outer
#' Gauss-Legendre grid; for each u, the integrals over the PFS time s (and the
#' progression time v) use an inner grid whose pieces end at u, at
#' t_p - t_o + u, at t_p minus the accrual breaks and (clock reset) at u minus
#' the post-progression breaks. The triple integrals of C3 and C4 are reduced
#' with Gamma(v, u) = int_0^v beta_P(s) Abar(min(t_p - s, t_o - u)) ds,
#' obtained by cumulative integration along the inner grid.
#'
#' @param scenario An \code{ms_scenario} object.
#' @param grp 0 (control group) or 1 (treatment group).
#' @param acc Output of \code{.acc_setup()}.
#' @param t_p,t_o Calendar times of the PFS and OS functionals (positive).
#' @param brk Break points of the intensities (\code{.scn_breaks()}).
#' @param hmax Maximum piece length.
#' @param n Number of nodes per piece.
#' @return A 3 x 3 matrix with rows for the PFS functionals and columns for
#'   the OS functionals (\code{psi}, \code{phi}, \code{dlt}).
#' @noRd
.cross_group <- function(scenario, grp, acc, t_p, t_o, brk, hmax, n = 8L) {
  keys <- c("psi", "phi", "dlt")
  model <- if (grp == 1L) scenario$treatment else scenario$control
  tm <- model$time
  tpost <- model$time_post
  forward <- model$clock == "forward"
  a_br <- acc$alpha
  shift <- t_o - t_p
  pre <- function(x) {
    list(sp = exp(-.pw_cumhaz(x, tm, model$prog + model$death)),
         a = .pw_value(x, tm, model$prog), b = .pw_value(x, tm, model$death),
         g = exp(-.pw_cumhaz(x, tm, model$dropout)),
         pp = .pw_value(x, tm, model$switch_prob))
  }
  # outer grid in the OS follow-up time u
  u_br <- c(brk, t_o - a_br, t_p, t_p - a_br, shift + brk)
  if (!forward) {
    u_br <- c(u_br, t_p + tpost)
  }
  go <- .gl_pieces(0, t_o, u_br, hmax = hmax, n = n)
  u <- go$x
  nu <- length(u)
  # inner grids in the PFS follow-up time s, from 0 to max(u, t_p)
  inner <- lapply(seq_len(nu), function(i) {
    ib <- c(u[i], brk, t_p - a_br, u[i] - shift)
    if (!forward) {
      ib <- c(ib, u[i] - tpost)
    }
    .gl_pieces(0, max(u[i], t_p), ib, hmax = hmax, n = n)
  })
  len <- vapply(inner, function(z) length(z$x), integer(1))
  npc <- vapply(inner, function(z) length(z$lo), integer(1))
  idx <- rep.int(seq_len(nu), len)
  s <- unlist(lapply(inner, `[[`, "x"), use.names = FALSE)
  ws <- unlist(lapply(inner, `[[`, "w"), use.names = FALSE)
  p_lo <- unlist(lapply(inner, `[[`, "lo"), use.names = FALSE)
  p_hi <- unlist(lapply(inner, `[[`, "hi"), use.names = FALSE)
  p_u <- rep.int(seq_len(nu), npc)
  # pieces below u (s < u) and from u on (s >= u)
  p_lt <- p_hi <= u[p_u]
  lt <- rep(p_lt, each = n)
  ge <- !lt
  ab <- .acc_cum(acc, pmin(t_p - s, t_o - u[idx]))
  ms <- pre(s)
  mu <- pre(u)
  kap <- ms$sp * ms$a
  # post-progression survival and death density for s < u
  sl <- s[lt]
  ul <- u[idx[lt]]
  tr <- if (forward) {
    .post_trans(model, sl, ul)
  } else {
    .post_trans(model, numeric(length(sl)), ul - sl)
  }
  ppl <- ms$pp[lt]
  s_post <- (1 - ppl) * (tr$r22 + tr$r23) + ppl * tr$r33
  f_post <- (1 - ppl) * (tr$r22 * tr$c_to + tr$r23 * tr$d_to) + ppl * tr$r33 * tr$d_to
  # cumulative integral along the inner grid of each u
  qm <- .gl_cumint(n)
  wref <- .gl_rule(n)$w
  half <- (p_hi - p_lo) / 2
  gamma_fun <- function(v) {
    fm <- matrix(v, nrow = n)
    tot <- half * colSums(wref * fm)
    start <- stats::ave(tot, p_u, FUN = cumsum) - tot
    list(node = rep(start, each = n) + as.vector(qm %*% fm) * rep(half, each = n),
         below = .sum_by(tot * p_lt, p_u, nu)[, 1L])
  }
  fp <- .fun_vals(scenario, s, "pfs", grp)
  wl <- ws[lt]
  il <- idx[lt]
  wg <- ws[ge]
  ig <- idx[ge]
  pr <- matrix(0, nu, 3L, dimnames = list(NULL, keys))
  pb <- pr
  for (k in keys) {
    rho_s <- fp$rho[, k]
    beta_s <- fp$beta[, k]
    if (k == "dlt") {
      gm_node <- numeric(length(s))
      gm_below <- numeric(nu)
    } else {
      gm <- gamma_fun(beta_s * ab)
      gm_node <- gm$node
      gm_below <- gm$below
    }
    t1 <- .sum_by(wl * rho_s[lt] * kap[lt] * f_post * ab[lt], il, nu)[, 1L]
    t2 <- .sum_by(wl * kap[lt] * f_post * gm_node[lt], il, nu)[, 1L]
    pr[, k] <- mu$g * (t1 + mu$sp * mu$b * gm_below + t2)
    t3 <- .sum_by(wl * rho_s[lt] * kap[lt] * s_post * ab[lt], il, nu)[, 1L]
    t4 <- .sum_by(wl * kap[lt] * s_post * gm_node[lt], il, nu)[, 1L]
    t5 <- .sum_by(wg * (rho_s[ge] * (ms$a[ge] + ms$b[ge]) + beta_s[ge]) *
                    ms$sp[ge] * ms$g[ge] * ab[ge], ig, nu)[, 1L]
    pb[, k] <- mu$g * (t3 + mu$sp * gm_below + t4) + t5
  }
  fo <- .fun_vals(scenario, u, "os", grp)
  res <- crossprod(pr, go$w * fo$rho) + crossprod(pb, go$w * fo$beta)
  # death without progression observed for both endpoints
  tmin <- min(t_p, t_o)
  gd <- .gl_pieces(0, tmin, c(brk, tmin - a_br), hmax = hmax, n = n)
  if (length(gd$x) > 0L) {
    md <- pre(gd$x)
    wd <- gd$w * md$sp * md$b * md$g * .acc_cum(acc, tmin - gd$x)
    res <- res + crossprod(.fun_vals(scenario, gd$x, "pfs", grp)$rho,
                           wd * .fun_vals(scenario, gd$x, "os", grp)$rho)
  }
  dimnames(res) <- list(keys, keys)
  res
}
