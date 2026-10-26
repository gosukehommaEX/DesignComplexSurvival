# Compare the graphical procedure with recycling of the level
# (ms_test(transition = ...), gs_prob_mtp() with upper_recycled; stage 3,
# candidate 1) with the decisions of lrstat::fseqbon().
#
# For each design below, the statistics of the two endpoints at their looks
# are drawn from the multivariate normal distribution with the mean and
# covariance of ms_power() (alternative), and with the mean set to 0 for both
# endpoints (global null hypothesis, "null") or for one of them ("null_first",
# "null_second"). Their one-sided p-values are passed to lrstat::fseqbon()
# with the weights, the transition matrix, the spending functions, the
# expected numbers of events as the information and the information
# fractions of the package as the spending times, with lookback = FALSE: as in
# the package, a hypothesis is tested again with its bound at the new level
# at the analysis at which the other is rejected, but not at earlier looks.
# fseqbon() has no futility bounds. In design L3 the futility stop of the
# package (the first hypothesis not rejected and its statistic below the
# non-binding futility bound at an interim look) is applied to the decisions
# of fseqbon() afterwards; without lookback, the decisions up to an analysis
# depend only on the statistics up to that analysis.
#
# Judged (4se_p, with the standard error sqrt(p (1 - p) / B) of the package
# value p and B simulated trials): the probabilities of gs_prob_mtp() (and,
# under the alternative, of ms_power()) of rejecting each hypothesis at each
# look, of rejecting the first, the second, both and at least one, and of
# stopping for futility, against the proportions of the simulated trials in
# which fseqbon() rejects. Judged (abs:1e-5): the efficacy bounds of
# ms_power() at the initial and at the recycled levels against
# lrstat::getBound(). Judged (abs:1e-4): the proportion of simulated trials
# in which the decisions of fseqbon() differ from those of the procedure
# applied directly with the bounds of ms_power() (graph_direct() below).
#
# Designs (scenarios of helpers.R):
#   L1  V3, looks at 252 and 336 OS events, both endpoints tested at both
#       looks; OS (O'Brien-Fleming type spending) and PFS (Pocock type),
#       weights 0.8 and 0.2, transition c(1, 1)
#   L2  V6, look 1 at 330 PFS events (PFS and OS tested), look 2 at 250 OS
#       events (OS tested); PFS and OS (both O'Brien-Fleming type), weights
#       0.8 and 0.2, transition c(1, 1)
#   L3  V3 with the looks of L1; OS (O'Brien-Fleming type, non-binding
#       futility bound by O'Brien-Fleming type beta spending with beta 0.2)
#       and PFS (Hwang-Shih-DeCani, gamma -4), weights 0.6 and 0.4,
#       transition c(0.5, 1)

if (!exists("val_scenarios", mode = "function")) {
  source(file.path("inst", "validation", "helpers.R"))
}

val_graph_designs <- function() {
  defs <- val_scenarios()
  looks_v3 <- ms_looks("os", fraction = c(252, 336) / 336)
  list(
    L1 = list(def = defs$V3, looks = looks_v3, events = 336,
              test = ms_test(alpha = 0.025, endpoints = c("os", "pfs"),
                             sf_eff = c("OF", "Pocock"), mtp = "bonferroni",
                             weights = c(0.8, 0.2), transition = c(1, 1))),
    L2 = list(def = defs$V6,
              looks = ms_looks(c("pfs", "os"), fraction = c(1, 1),
                               test = list(c("pfs", "os"), "os")),
              events = c(pfs = 330, os = 250),
              test = ms_test(alpha = 0.025, endpoints = c("pfs", "os"),
                             sf_eff = c("OF", "OF"), mtp = "bonferroni",
                             weights = c(0.8, 0.2), transition = c(1, 1))),
    L3 = list(def = defs$V3, looks = looks_v3, events = 336,
              test = ms_test(alpha = 0.025, endpoints = c("os", "pfs"),
                             sf_eff = c("OF", "HSD"), param_eff = list(NULL, -4),
                             futility = "beta", sf_fut = "OF", beta = 0.2,
                             binding = FALSE, mtp = "bonferroni",
                             weights = c(0.6, 0.4), transition = c(0.5, 1)))
  )
}

# Spending function of the package as type and parameter of lrstat
val_lrstat_sf <- function(sf, param) {
  type <- c(OF = "sfOF", Pocock = "sfP", HSD = "sfHSD", power = "sfKD")[sf]
  if (is.na(type)) {
    stop("No lrstat spending function for sf_eff = \"", sf, "\".", call. = FALSE)
  }
  list(type = unname(type), param = if (is.null(param)) NA_real_ else as.numeric(param))
}

# The procedure applied directly to the statistics z (one row per trial,
# columns: the looks of the first endpoint, then those of the second) with
# the initial bounds up, the recycled bounds upr and the futility bounds
# lower of the first endpoint. Returns the analysis at which each hypothesis
# is rejected (0 for none) and the analysis of a futility stop (0 for none).
graph_direct <- function(z, an, up, upr, lower) {
  nb <- nrow(z)
  col <- list(seq_along(an[[1L]]), length(an[[1L]]) + seq_along(an[[2L]]))
  rej <- matrix(0L, nb, 2L)
  stopped <- rep(FALSE, nb)
  fut <- integer(nb)
  for (a in seq_len(max(unlist(an)))) {
    k <- c(match(a, an[[1L]]), match(a, an[[2L]]))
    act <- !stopped
    newr <- matrix(FALSE, nb, 2L)
    for (pass in 1:2) {
      for (j in 1:2) {
        if (is.na(k[j])) {
          next
        }
        other <- rej[, 3L - j] > 0L | newr[, 3L - j]
        b <- ifelse(other, upr[[j]][k[j]], up[[j]][k[j]])
        newr[, j] <- newr[, j] | (act & rej[, j] == 0L & z[, col[[j]][k[j]]] >= b)
      }
    }
    rej[newr & rej == 0L] <- a
    if (!is.null(lower) && !is.na(k[1L]) && k[1L] < length(an[[1L]])) {
      f <- act & rej[, 1L] == 0L & z[, col[[1L]][k[1L]]] < lower[k[1L]]
      fut[f] <- a
      stopped <- stopped | f
    }
    stopped <- stopped | (rej[, 1L] > 0L & rej[, 2L] > 0L)
  }
  list(rej = rej, fut = fut)
}

compare_graph_lrstat <- function(B = 100000, seed = 20261026) {
  if (!requireNamespace("lrstat", quietly = TRUE)) {
    message("lrstat is not installed: the comparison with lrstat::fseqbon() is not made.")
    return(NULL)
  }
  ver <- as.character(utils::packageVersion("lrstat"))
  if (!all(c("k1", "lookback") %in% names(formals(lrstat::fseqbon)))) {
    message("lrstat ", ver, ": fseqbon() has no arguments 'k1' and 'lookback' ",
            "(lrstat 0.3.4 or later is needed); the comparison is not made.")
    return(NULL)
  }
  message("Comparison with lrstat ", ver, " (fseqbon() with lookback = FALSE, getBound()).")
  des <- val_graph_designs()
  out <- list()
  emit <- function(scn, hyp, look, time, quantity, endpoint, variant, criterion, simulated,
                   package, se = NA_real_, n = NA_real_, note = "") {
    out[[length(out) + 1L]] <<- data.frame(
      id = paste("lrstat_fseqbon", scn, hyp, "design", look, ifelse(is.na(time), "", format(time)),
                 quantity, endpoint, "graph", variant, sep = "|"),
      source = "lrstat_fseqbon", scenario = scn, hypothesis = hyp, look_type = "design",
      look = as.character(look), time = time, quantity = quantity, endpoint = endpoint,
      group = "graph", variant = variant, criterion = criterion, simulated = simulated,
      se = se, n = n, package = package, note = note, stringsAsFactors = FALSE)
  }
  for (id in seq_along(des)) {
    nm <- names(des)[id]
    d <- des[[nm]]
    ts <- d$test
    pw <- ms_power(d$def$scenario, n = d$def$n, looks = d$looks, test = ts, events = d$events)
    eps <- ts$endpoints
    lk <- pw$looks
    rows_j <- lapply(eps, function(e) which(lk$endpoint == e))
    an <- lapply(rows_j, function(r) lk$look[r])
    up <- lapply(rows_j, function(r) lk$z_upper[r])
    upr <- lapply(rows_j, function(r) lk$z_upper_recycled[r])
    lower <- if (ts$futility == "none") NULL else lk$z_lower[rows_j[[1L]]]
    kmax <- max(unlist(an))
    w <- ts$weights
    g <- ts$transition
    lev <- ts$alpha * (w + rev(g * w))
    G <- matrix(c(0, g[1L], g[2L], 0), 2L, 2L, byrow = TRUE)
    sf <- lapply(1:2, function(j) val_lrstat_sf(ts$sf_eff[j], ts$param_eff[[j]]))
    types <- vapply(sf, `[[`, "", "type")
    pars <- vapply(sf, `[[`, 0, "param")
    atime <- lk$time[match(seq_len(kmax), lk$look)]

    # bounds at the initial and at the recycled levels
    for (j in 1:2) {
      r <- rows_j[[j]]
      ir <- lk$events[r] / lk$events[r[length(r)]]
      for (v in c("initial", "recycled")) {
        gb <- lrstat::getBound(k = length(r), informationRates = ir,
                               alpha = if (v == "initial") ts$alpha * w[j] else lev[j],
                               typeAlphaSpending = types[j],
                               parameterAlphaSpending = pars[j], spendingTime = lk$info_frac[r])
        pk <- if (v == "initial") up[[j]] else upr[[j]]
        for (k in seq_along(r)) {
          emit(nm, "design", an[[j]][k], atime[an[[j]][k]],
               if (v == "initial") "bound" else "bound_recycled", eps[j], "getBound", "abs:1e-5",
               gb[k], pk[k], se = 0)
        }
      }
    }

    # inputs of fseqbon() that do not depend on the statistics: information
    # and spending time at the analyses at which an endpoint is not tested are
    # those of its previous look (or of its first look)
    incid <- matrix(FALSE, kmax, 2L)
    info_a <- matrix(NA_real_, kmax, 2L)
    st_a <- matrix(NA_real_, kmax, 2L)
    for (j in 1:2) {
      r <- rows_j[[j]]
      incid[an[[j]], j] <- TRUE
      for (a in seq_len(kmax)) {
        k <- max(c(0L, which(an[[j]] <= a)))
        if (k == 0L) {
          k <- 1L
        }
        info_a[a, j] <- lk$events[r[k]]
        st_a[a, j] <- lk$info_frac[r[k]]
      }
    }
    maxinfo <- vapply(rows_j, function(r) lk$events[r[length(r)]], 0)
    fs <- function(pm) {
      nb <- nrow(pm) / kmax
      lrstat::fseqbon(w = w, G = G, alpha = ts$alpha, kMax = kmax, typeAlphaSpending = types,
                      parameterAlphaSpending = pars, maxInformation = maxinfo,
                      incidenceMatrix = incid, k1 = kmax, p = pm,
                      information = info_a[rep(seq_len(kmax), nb), , drop = FALSE],
                      spendingTime = st_a[rep(seq_len(kmax), nb), , drop = FALSE],
                      lookback = FALSE)
    }

    mu <- pw$mean
    i1 <- seq_along(an[[1L]])
    i2 <- length(an[[1L]]) + seq_along(an[[2L]])
    hyps <- list(alternative = mu, null = 0 * mu, null_first = replace(mu, i1, 0),
                 null_second = replace(mu, i2, 0))
    ch <- chol(pw$sigma)
    for (ih in seq_along(hyps)) {
      h <- names(hyps)[ih]
      gp <- gs_prob_mtp(up, lower = lower, mean = hyps[[h]], sigma = pw$sigma, analysis = an,
                        mtp = "bonferroni", upper_recycled = upr)
      set.seed(seed + 100L * id + ih)
      z <- matrix(stats::rnorm(B * length(mu)), B) %*% ch
      z <- sweep(z, 2L, hyps[[h]], "+")
      pv <- stats::pnorm(z, lower.tail = FALSE)
      # p-values by trial and analysis (rows: trial 1 analyses 1..kmax, trial 2 ...)
      pm <- matrix(1, B * kmax, 2L)
      for (j in 1:2) {
        cj <- if (j == 1L) i1 else i2
        for (k in seq_along(an[[j]])) {
          pm[(seq_len(B) - 1L) * kmax + an[[j]][k], j] <- pv[, cj[k]]
        }
      }
      res <- fs(pm)
      # shape of the result: checked against single-trial calls
      nchk <- min(20L, B)
      chk <- t(vapply(seq_len(nchk), function(b) {
        as.numeric(fs(pm[(b - 1L) * kmax + seq_len(kmax), , drop = FALSE]))
      }, numeric(2)))
      cand <- list(if (is.matrix(res) && nrow(res) == B) res,
                   if (length(res) == 2L * B) matrix(res, ncol = 2L, byrow = TRUE),
                   if (length(res) == 2L * B) matrix(res, ncol = 2L))
      ok <- vapply(cand, function(x) !is.null(x) && all(x[seq_len(nchk), ] == chk), logical(1))
      if (!any(ok)) {
        stop("Cannot match the result of lrstat::fseqbon() for ", B,
             " trials with single-trial calls (design ", nm, ").", call. = FALSE)
      }
      res <- cand[[which(ok)[1L]]]
      r1 <- as.integer(res[, 1L])
      r2 <- as.integer(res[, 2L])
      fut <- integer(B)
      if (!is.null(lower)) {
        for (k in seq_len(length(an[[1L]]) - 1L)) {
          a <- an[[1L]][k]
          f <- fut == 0L & (r1 == 0L | r1 > a) & z[, i1[k]] < lower[k]
          fut[f] <- a
        }
        stop_a <- ifelse(fut > 0L, fut, kmax + 1L)
        r1[r1 > stop_a] <- 0L
        r2[r2 > stop_a] <- 0L
      }
      dr <- graph_direct(z, an, up, upr, lower)
      emit(nm, h, "", NA_real_, "decisions_differ", "", "direct", "abs:1e-4",
           mean(dr$rej[, 1L] != r1 | dr$rej[, 2L] != r2 | dr$fut != fut), 0, se = 0)
      rr <- list(r1, r2)
      for (j in 1:2) {
        for (k in seq_along(an[[j]])) {
          a <- an[[j]][k]
          pk <- gp$looks$prob_upper[gp$looks$endpoint == j & gp$looks$analysis == a]
          emit(nm, h, a, atime[a], "reject_look", eps[j], "gs_prob_mtp", "4se_p",
               mean(rr[[j]] == a), pk, n = B)
        }
      }
      if (!is.null(lower)) {
        for (k in seq_len(length(an[[1L]]) - 1L)) {
          a <- an[[1L]][k]
          pk <- gp$looks$prob_lower[gp$looks$endpoint == 1L & gp$looks$analysis == a]
          emit(nm, h, a, atime[a], "futility_look", eps[1L], "gs_prob_mtp", "4se_p",
               mean(fut == a), pk, n = B)
        }
      }
      sim <- c(mean(r1 > 0L), mean(r2 > 0L), mean(r1 > 0L & r2 > 0L), mean(r1 > 0L | r2 > 0L))
      names(sim) <- c(eps, "both", "any")
      pk <- c(gp$reject[["first"]], gp$reject[["second"]], gp$reject[["both"]], gp$reject[["any"]])
      for (q in seq_along(sim)) {
        emit(nm, h, "", NA_real_, "reject", names(sim)[q], "gs_prob_mtp", "4se_p", sim[[q]], pk[q],
             n = B)
        if (h == "alternative") {
          emit(nm, h, "", NA_real_, "reject", names(sim)[q], "ms_power", "4se_p", sim[[q]],
               pw$power[[names(sim)[q]]], n = B)
        }
      }
      message(sprintf("  %s %s: fseqbon decisions differ from the direct procedure in %d of %d trials",
                      nm, h, sum(dr$rej[, 1L] != r1 | dr$rej[, 2L] != r2 | dr$fut != fut), B))
    }
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
