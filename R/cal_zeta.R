#' Covariance of the standardized OS log-rank statistics when switching
#' starts at a calendar time
#'
#' The counterpart of \code{.zeta_cov()} for one endpoint (OS) and scenarios
#' with \code{.is_calendar()} (calendar-time switching note, Section 7). At
#' each look the functionals psi, phi and dlt use the proportion
#' \code{q1} and the pooled hazard \code{lbar} of that look
#' (\code{.cal_look()}); their means are integrals against the group
#' functions, their products \code{.cal_pair()}. The standardized statistic
#' has the influence function zeta of equation (2) of the covariance note;
#' the term of a look driven by OS events uses the time derivatives of the
#' mean, the information and the events, obtained by differentiating the
#' upper limit t - s of the entry integral (\code{.cal_dt()}).
#'
#' The variance options give: \code{"studentized"} the covariance of zeta;
#' \code{"null"} its correlation matrix (variances 1); \code{"alternative"}
#' the covariance of psi divided by the square roots of the information.
#' Under the null hypothesis the three coincide; the correlation between
#' looks then differs in general from the square root of the ratio of the
#' information, because the hazard of a patient depends on the entry time.
#'
#' @param scenario An \code{ms_scenario} object with \code{.is_calendar()}.
#' @param acc Output of \code{.acc_setup()}.
#' @param time Increasing calendar times of the looks.
#' @param trigger \code{NULL} or, per look, \code{"os"} for a look driven by
#'   OS events and \code{NA} otherwise.
#' @param variance \code{"null"}, \code{"alternative"} or
#'   \code{"studentized"}.
#' @return A list with \code{sigma}, the per-patient score means \code{mu} and
#'   information \code{v}.
#' @noRd
.cal_zeta <- function(scenario, acc, time, trigger = NULL, variance = "studentized") {
  keys <- c("psi", "phi", "dlt")
  k <- length(time)
  if (is.null(trigger)) {
    trigger <- rep(NA_character_, k)
  }
  trigger <- rep_len(as.character(trigger), k)
  if (any(!is.na(trigger) & trigger != "os")) {
    stop("Internal error: looks driven by PFS events with calendar-time switching.",
         call. = FALSE)
  }
  p <- scenario$prob
  looks <- lapply(time, function(t) .cal_look(scenario, acc, t))
  mu <- vapply(looks, `[[`, numeric(1), "mu")
  v <- vapply(looks, `[[`, numeric(1), "v")
  if (any(!(v > 0))) {
    stop("The information of the log-rank statistic is zero at some times.", call. = FALSE)
  }
  funs <- lapply(looks, function(lk) {
    lapply(1:2, function(j) {
      gj <- (j - 1) - lk$q1
      rho <- cbind(psi = gj, phi = lk$q1 * (1 - lk$q1), dlt = 1)
      beta <- cbind(psi = -gj * lk$lbar, phi = (1 - 2 * lk$q1) * gj * lk$lbar, dlt = 0)
      list(rho = rho, beta = beta, pr = .pw_poly(lk$g, rho, 12L),
           pb = .pw_poly(lk$g, beta, 12L))
    })
  })
  # means of the functionals
  avg <- matrix(0, k, 3L, dimnames = list(NULL, keys))
  for (i in seq_len(k)) {
    lk <- looks[[i]]
    for (j in 1:2) {
      fj <- funs[[i]][[j]]
      wt <- lk$g$w * lk$G[[j]]
      avg[i, ] <- avg[i, ] + p[[j]] * (colSums(fj$rho * (wt * lk$D[[j]])) +
                                         colSums(fj$beta * (wt * lk$Y[[j]])))
    }
  }
  # covariance of the elementary functionals (key within look)
  ne <- 3L * k
  cv <- matrix(0, ne, ne)
  for (i in seq_len(k)) {
    for (l in i:k) {
      pm <- .cal_pair(scenario, acc, looks[[i]], looks[[l]], funs[[i]], funs[[l]])
      blk <- pm - outer(avg[i, ], avg[l, ])
      ri <- 3L * (i - 1L) + 1:3
      rl <- 3L * (l - 1L) + 1:3
      cv[ri, rl] <- blk
      cv[rl, ri] <- base::t(blk)
    }
  }
  # time derivatives for looks driven by OS events
  mdot <- numeric(k)
  ddot <- numeric(k)
  for (i in which(!is.na(trigger))) {
    lk <- looks[[i]]
    s <- lk$g$x
    dt <- .cal_dt(scenario, acc, lk)
    dy <- dt$dy
    dd <- dt$dd
    y <- lapply(1:2, function(j) p[[j]] * lk$G[[j]] * lk$Y[[j]])
    d <- lapply(1:2, function(j) p[[j]] * lk$G[[j]] * lk$D[[j]])
    ys <- y[[1L]] + y[[2L]]
    ok <- ys > 0
    num <- y[[1L]] * d[[2L]] - y[[2L]] * d[[1L]]
    dnum <- dy[[1L]] * d[[2L]] + y[[1L]] * dd[[2L]] - dy[[2L]] * d[[1L]] - y[[2L]] * dd[[1L]]
    wv <- y[[1L]] * y[[2L]] * (d[[1L]] + d[[2L]])
    dwv <- (dy[[1L]] * y[[2L]] + y[[1L]] * dy[[2L]]) * (d[[1L]] + d[[2L]]) +
      y[[1L]] * y[[2L]] * (dd[[1L]] + dd[[2L]])
    du <- numeric(length(s))
    dv <- numeric(length(s))
    du[ok] <- (dnum[ok] * ys[ok] - num[ok] * (dy[[1L]][ok] + dy[[2L]][ok])) / ys[ok]^2
    dv[ok] <- (dwv[ok] * ys[ok] - 2 * wv[ok] * (dy[[1L]][ok] + dy[[2L]][ok])) / ys[ok]^3
    mu_dot <- sum(lk$g$w * du)
    v_dot <- sum(lk$g$w * dv)
    mdot[i] <- -mu_dot / sqrt(v[i]) + mu[i] * v_dot / (2 * v[i]^1.5)
    ddot[i] <- sum(lk$g$w * (dd[[1L]] + dd[[2L]]))
  }
  coef <- matrix(0, k, ne)
  ipsi <- 3L * (seq_len(k) - 1L) + 1L
  for (i in seq_len(k)) {
    coef[i, ipsi[i]] <- -1 / sqrt(v[i])
    coef[i, ipsi[i] + 1L] <- mu[i] / (2 * v[i]^1.5)
    if (!is.na(trigger[i])) {
      coef[i, ipsi[i] + 2L] <- -mdot[i] / ddot[i]
    }
  }
  s_stud <- coef %*% cv %*% base::t(coef)
  s_stud <- (s_stud + base::t(s_stud)) / 2
  sigma <- if (variance == "null") {
    sd <- sqrt(diag(s_stud))
    s_stud / outer(sd, sd)
  } else if (variance == "alternative") {
    cv[ipsi, ipsi, drop = FALSE] / sqrt(outer(v, v))
  } else {
    s_stud
  }
  list(sigma = sigma, mu = mu, v = v)
}
