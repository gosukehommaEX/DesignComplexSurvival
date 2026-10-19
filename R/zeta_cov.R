#' Covariance of standardized log-rank statistics of one or two endpoints
#'
#' Covariance matrix of Z = -U / sqrt(V) for a set of items, each given by an
#' endpoint, a calendar time and the endpoint whose events determine that time
#' (\code{NA} for a calendar time). The influence function of Z for an item
#' of endpoint k at time t is (covariance note, Planning/covariance_note)
#' zeta = -(psi - mu) / sqrt(v) + mu / (2 v^(3/2)) (phi - v)
#'   - (m' / d_D') (dlt_D - d_D),
#' where mu and v are the per-patient mean of the score and information,
#' m = -mu / sqrt(v), d_D the per-patient expected events of the trigger
#' endpoint D and ' the derivative with respect to calendar time; the last
#' term is present only for event-driven times. The moments of psi, phi and
#' dlt come from \code{.same_moments()} and \code{.cross_moments()}, the
#' derivatives from \code{.acc_integral_deriv()}.
#'
#' With \code{variance = "studentized"} the covariance of zeta is returned.
#' With \code{"null"}, items of the same endpoint have unit variances and
#' correlations sqrt(v_j / v_k) (as in the first stage), and items of
#' different endpoints the correlation of zeta. With \code{"alternative"},
#' items of the same endpoint have the covariance of the score of Luo et al.
#' (2019) divided by sqrt(v_j v_k) (\code{.lr_alt_cov()}, as in the first
#' stage), and items of different endpoints Cov(psi_P, psi_O) divided by
#' sqrt(v_P v_O).
#'
#' @param scenario An \code{ms_scenario} object.
#' @param acc Output of \code{.acc_setup()}.
#' @param endpoint Endpoint of each item (\code{"os"} or \code{"pfs"}).
#' @param time Calendar time of each item (positive).
#' @param trigger Endpoint whose events determine the time of each item, or
#'   \code{NA} for a calendar time (\code{NULL}: all calendar times).
#' @param variance \code{"null"}, \code{"alternative"} or
#'   \code{"studentized"}.
#' @return A list with \code{sigma} (covariance matrix of the items),
#'   \code{mu} and \code{v} (per-patient score mean and information of each
#'   item).
#' @noRd
.zeta_cov <- function(scenario, acc, endpoint, time, trigger = NULL,
                      variance = "studentized") {
  keys <- c("psi", "phi", "dlt")
  n <- acc$n
  k <- length(time)
  if (is.null(trigger)) {
    trigger <- rep(NA_character_, k)
  }
  trigger <- rep_len(as.character(trigger), k)
  eps <- intersect(c("pfs", "os"), c(endpoint, trigger[!is.na(trigger)]))
  tset <- lapply(eps, function(e) {
    sort(unique(c(time[endpoint == e], time[!is.na(trigger) & trigger == e])))
  })
  names(tset) <- eps
  brk <- .scn_breaks(scenario)
  # elementary functionals: endpoint, functional and index of the time
  el <- do.call(rbind, lapply(eps, function(e) {
    nt <- length(tset[[e]])
    data.frame(endpoint = e, key = rep(keys, each = nt), it = rep(seq_len(nt), 3L),
               stringsAsFactors = FALSE)
  }))
  el_id <- function(e, key, it) {
    which(el$endpoint == e & el$key == key & el$it == it)
  }
  ne <- nrow(el)
  cv <- matrix(0, ne, ne)
  sm <- list()
  mom <- list()
  for (e in eps) {
    tt <- tset[[e]]
    sm[[e]] <- .same_moments(scenario, acc, tt, e)
    core <- .lr_core(scenario, acc, tt, e, "direct")
    hfun <- function(s) {
      .lr_integrands(scenario, s, e)[, c("ev0", "ev1", "u", "v"), drop = FALSE]
    }
    der <- .acc_integral_deriv(acc, tt, hfun, brk, .scn_hmax(scenario, max(c(tt, 1))))
    mu <- core$score_mean / n
    v <- core$info / n
    if (any(!(v > 0))) {
      stop("The information of the log-rank statistic is zero at some times.",
           call. = FALSE)
    }
    mom[[e]] <- list(mu = mu, v = v,
                     mdot = -der[, "u"] / sqrt(v) + mu * der[, "v"] / (2 * v^1.5),
                     ddot = der[, "ev0"] + der[, "ev1"])
    ii <- which(el$endpoint == e)
    for (a in ii) {
      for (b in ii) {
        cv[a, b] <- sm[[e]]$prod[el$it[a], el$it[b], el$key[a], el$key[b]] -
          sm[[e]]$mean[el$it[a], el$key[a]] * sm[[e]]$mean[el$it[b], el$key[b]]
      }
    }
  }
  if (length(eps) == 2L) {
    for (i1 in seq_along(tset$pfs)) {
      for (i2 in seq_along(tset$os)) {
        cm <- .cross_moments(scenario, acc, tset$pfs[i1], tset$os[i2])
        for (kp in keys) {
          for (ko in keys) {
            a <- el_id("pfs", kp, i1)
            b <- el_id("os", ko, i2)
            cv[a, b] <- cm[kp, ko] - sm$pfs$mean[i1, kp] * sm$os$mean[i2, ko]
            cv[b, a] <- cv[a, b]
          }
        }
      }
    }
  }
  coef <- matrix(0, k, ne)
  mu_i <- numeric(k)
  v_i <- numeric(k)
  ipsi <- integer(k)
  for (i in seq_len(k)) {
    e <- endpoint[i]
    it <- match(time[i], tset[[e]])
    mu_i[i] <- mom[[e]]$mu[it]
    v_i[i] <- mom[[e]]$v[it]
    ipsi[i] <- el_id(e, "psi", it)
    coef[i, ipsi[i]] <- -1 / sqrt(v_i[i])
    coef[i, el_id(e, "phi", it)] <- mu_i[i] / (2 * v_i[i]^1.5)
    if (!is.na(trigger[i])) {
      d <- trigger[i]
      jd <- match(time[i], tset[[d]])
      jdl <- el_id(d, "dlt", jd)
      coef[i, jdl] <- coef[i, jdl] - mom[[e]]$mdot[it] / mom[[d]]$ddot[jd]
    }
  }
  s_stud <- coef %*% cv %*% base::t(coef)
  s_stud <- (s_stud + base::t(s_stud)) / 2
  sigma <- s_stud
  same <- outer(endpoint, endpoint, "==")
  if (variance == "null") {
    sd <- sqrt(diag(s_stud))
    sigma <- s_stud / outer(sd, sd)
    vv <- sqrt(outer(v_i, v_i, pmin) / outer(v_i, v_i, pmax))
    sigma[same] <- vv[same]
  } else if (variance == "alternative") {
    sigma <- cv[ipsi, ipsi, drop = FALSE] / sqrt(outer(v_i, v_i))
    for (e in unique(endpoint)) {
      ie <- which(endpoint == e)
      sigma[ie, ie] <- .lr_alt_cov(scenario, acc, time[ie], e, mu_i[ie]) /
        sqrt(outer(v_i[ie], v_i[ie]))
    }
  }
  list(sigma = sigma, mu = mu_i, v = v_i)
}
