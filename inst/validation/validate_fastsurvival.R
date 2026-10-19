# Patient-level simulation of V1-V5 with FastSurvival (Phase 3 plan,
# sections 4.1 and 12).
#
# For each scenario and hypothesis, 20000 trials are simulated in 10 batches
# of 2000 (independent dqrng streams of one seed). Only summaries are kept,
# in the format of python/simulate_ms.py, and written to
# results/sim_fastsurvival.csv:
#   calendar looks: group-wise OS and PFS events and dropouts, log-rank Z of
#     OS and PFS (sign changed so that it is positive when the treatment is
#     better), and the correlation of Z between the looks;
#   event looks (calendar time of the d-th OS death in each trial, from
#     cutoff_fast()): look time, group-wise events, Z, correlation and the
#     rejection of the O'Brien-Fleming design;
#   latent quantities (alternative hypothesis, no dropout): OS and PFS
#     survival and the state occupation probabilities at times since
#     randomization. FastSurvival switches only after an observed
#     progression, so these come from a separate simulation without dropout
#     (2000 trials).
# FastSurvival allocates fixed group sizes and places a fixed number of
# patients in each accrual interval.

if (!exists("val_scenarios", mode = "function")) {
  source(file.path("inst", "validation", "helpers.R"))
}

fs_available <- requireNamespace("FastSurvival", quietly = TRUE)
fs_version <- if (fs_available) as.character(utils::packageVersion("FastSurvival")) else NA_character_

# Arguments of simdata_fast() for each scenario; hyp = "null" gives the
# treatment group the control parameters
fs_args <- function(name, hyp) {
  two <- function(ctl, trt) if (hyp == "alternative") list(ctl, trt) else list(ctl, ctl)
  l2 <- log(2)
  switch(name,
    V1 = list(n = 400, alloc = c(1, 1), a.time = c(0, 20), a.rate = 20,
              e.hazard = two(l2 / 12, 0.7 * l2 / 12), d.hazard = list(0.005, 0.005)),
    V2 = list(n = 300, alloc = c(1, 2), a.time = c(0, 3), a.rate = c(10, 30),
              e.hazard = two(c(0.05, 0.04), c(0.05, 0.025)), e.time = c(0, 6, Inf),
              d.hazard = if (hyp == "alternative") list(0.01, 0.012) else list(0.01, 0.01)),
    V3 = list(n = 482, alloc = c(1, 1), a.time = c(0, 23), a.prop = 1,
              h01.hazard = two(l2 / 3 - l2 / 9, l2 / 4.3 - l2 / 12.9),
              h02.hazard = two(l2 / 9, l2 / 12.9), h12.hazard = two(l2 / 9, l2 / 12.9),
              d.hazard = -log(0.95) / 12),
    V4 = list(n = 500, alloc = c(1, 1), a.time = c(0, 6, 18), a.prop = c(0.3, 0.7),
              h01.hazard = two(c(0.10, 0.08), c(0.06, 0.05)), h01.time = c(0, 4, Inf),
              h02.hazard = two(0.01, 0.008),
              h12.hazard = two(c(0.05, 0.06), c(0.045, 0.05)), h12.time = c(0, 3, Inf),
              switch.prop = two(0.3, 0),
              h12.switch.hazard = two(c(0.03, 0.035), c(0.045, 0.05)),
              h12.switch.time = c(0, 3, Inf), d.hazard = 0.005),
    V5 = {
      lc <- log(2)
      he <- c(0.6, 0.8) * lc
      dh <- function(x) -log(1 - x) / 4
      list(n = 292, alloc = c(1, 1), a.time = c(0, 2), a.prop = 1,
           h01.hazard = two(dh(0.1), dh(0.1)),
           h02.hazard = two(c(lc, lc), he), h02.time = c(0, 2, Inf),
           h12.hazard = two(he, c(lc, lc)), h12.time = c(0, 2, Inf),
           d.hazard = dh(0.2))
    })
}

fs_simulate <- function(name, def, reps_batch = 2000L, n_batch = 10L, seed = 20261019L) {
  rows <- list()
  idm <- !is.null(fs_args(name, "alternative")$h01.hazard)
  eps <- if (idm) c("os", "pfs") else "os"
  cols <- list(os = if (idm) c("e2_tte", "e2_event") else c("tte", "event"),
               pfs = c("e1_tte", "e1_event"))
  sid <- match(name, c("V1", "V2", "V3", "V4", "V5"))
  for (hyp in c("alternative", "null")) {
    args <- fs_args(name, hyp)
    st <- list()
    for (b in seq_len(n_batch)) {
      d <- do.call(FastSurvival::simdata_fast,
                   c(list(nsim = reps_batch, seed = as.integer(seed + 10L * sid + (hyp == "null")), stream = b),
                     args))
      trt <- d$group != 1
      add <- function(key, v) st[[key]] <<- c(st[[key]], v)
      looks <- list()
      for (k in seq_along(def$cal)) {
        looks[[length(looks) + 1L]] <- list(tag = paste0("calendar|", k), cut = rep(def$cal[k], reps_batch))
      }
      if (!is.null(def$events)) {
        cutm <- FastSurvival::cutoff_fast(d, event.looks = def$events, tte.col = cols$os[1],
                                          event.col = cols$os[2])
        sims <- as.integer(rownames(cutm))
        for (k in seq_along(def$events)) {
          cut <- cutm[match(seq_len(reps_batch), sims), k]
          add(paste0("event|", k, "|time"), cut)
          looks[[length(looks) + 1L]] <- list(tag = paste0("event|", k), cut = cut)
        }
      }
      for (lk in looks) {
        cut_i <- lk$cut[d$sim]
        cut_i[is.na(cut_i)] <- Inf
        fu <- cut_i - d$accrual_time
        enr <- fu >= 0
        for (ep in eps) {
          tt <- d[[cols[[ep]][1]]]
          ev <- d[[cols[[ep]][2]]] == 1
          obs_ev <- enr & ev & tt <= fu
          obs_do <- enr & !ev & is.finite(tt) & tt <= fu
          m <- rowsum(cbind(ec = obs_ev & !trt, et = obs_ev & trt, dc = obs_do & !trt,
                            dt = obs_do & trt) * 1, d$sim, reorder = TRUE)
          add(paste0(lk$tag, "|events_", ep, "_control"), m[, "ec"])
          add(paste0(lk$tag, "|events_", ep, "_treatment"), m[, "et"])
          add(paste0(lk$tag, "|dropouts_", ep, "_control"), m[, "dc"])
          add(paste0(lk$tag, "|dropouts_", ep, "_treatment"), m[, "dt"])
          de <- data.frame(sim = d$sim, group = d$group, accrual_time = d$accrual_time,
                           tte = tt, event = as.numeric(ev))
          cutv <- matrix(lk$cut, ncol = 1, dimnames = list(seq_len(reps_batch), NULL))
          a <- FastSurvival::analysis_fast(de, control = 1, cutoff.looks = cutv)
          a <- a[order(a$sim), ]
          add(paste0(lk$tag, "|z_", ep), -a$logrank.z)
        }
      }
      message(sprintf("  %s %s batch %d of %d", name, hyp, b, n_batch))
    }
    tags <- c(paste0("calendar|", seq_along(def$cal)),
              if (!is.null(def$events)) paste0("event|", seq_along(def$events)))
    for (tg in tags) {
      lt <- sub("\\|.*", "", tg)
      k <- as.integer(sub(".*\\|", "", tg))
      tval <- if (lt == "calendar") def$cal[k] else NA_real_
      tgt <- if (lt == "event") def$events[k] else NA_real_
      if (lt == "event") {
        rows[[length(rows) + 1L]] <- val_sim_row("fastsurvival", name, hyp, lt, k, NA_real_, tgt, "time",
                                                 "os", "", st[[paste0(tg, "|time")]])
      }
      for (ep in eps) {
        for (g in c("control", "treatment")) {
          for (q in c("events", "dropouts")) {
            rows[[length(rows) + 1L]] <- val_sim_row("fastsurvival", name, hyp, lt, k, tval, tgt, q, ep, g,
                                                     st[[paste0(tg, "|", q, "_", ep, "_", g)]])
          }
        }
        rows[[length(rows) + 1L]] <- val_sim_row("fastsurvival", name, hyp, lt, k, tval, tgt, "z", ep, "",
                                                 st[[paste0(tg, "|z_", ep)]])
      }
    }
    for (lt in unique(sub("\\|.*", "", tags))) {
      for (ep in eps) {
        z1 <- st[[paste0(lt, "|1|z_", ep)]]
        z2 <- st[[paste0(lt, "|2|z_", ep)]]
        ok <- is.finite(z1) & is.finite(z2)
        r <- stats::cor(z1[ok], z2[ok])
        rows[[length(rows) + 1L]] <- data.frame(
          source = "fastsurvival", scenario = name, hypothesis = hyp, look_type = lt, look = "1-2",
          time = NA_real_, target_events = NA_real_, quantity = "z_corr", endpoint = ep, group = "",
          mean = r, sd = NA_real_, se = (1 - r^2) / sqrt(sum(ok)), n = sum(ok), stringsAsFactors = FALSE)
      }
    }
    if (!is.null(def$events)) {
      fr <- def$events / def$events[length(def$events)]
      b <- gs_bounds(fr, alpha = 0.025, sf_eff = "OF")$upper
      z1 <- st[["event|1|z_os"]]
      z2 <- st[["event|2|z_os"]]
      r1 <- !is.na(z1) & z1 >= b[1]
      ra <- r1 | (!is.na(z2) & z2 >= b[2])
      for (q in c("reject_look1", "reject_any")) {
        x <- if (q == "reject_look1") r1 else ra
        p <- mean(x)
        rows[[length(rows) + 1L]] <- data.frame(
          source = "fastsurvival", scenario = name, hypothesis = hyp, look_type = "event", look = "",
          time = NA_real_, target_events = NA_real_, quantity = q, endpoint = "os", group = "",
          mean = p, sd = NA_real_, se = sqrt(p * (1 - p) / length(x)), n = length(x),
          stringsAsFactors = FALSE)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        source = "fastsurvival", scenario = name, hypothesis = hyp, look_type = "event", look = "",
        time = NA_real_, target_events = NA_real_, quantity = "not_reached", endpoint = "os", group = "",
        mean = sum(is.na(st[[paste0("event|", length(def$events), "|time")]])), sd = NA_real_,
        se = NA_real_, n = length(z1), stringsAsFactors = FALSE)
    }
  }
  # latent quantities without dropout (alternative)
  args <- fs_args(name, "alternative")
  args$d.hazard <- 0
  d <- do.call(FastSurvival::simdata_fast, c(list(nsim = reps_batch, seed = as.integer(seed + 10L * sid + 5L)), args))
  for (g in c("control", "treatment")) {
    gi <- if (g == "control") d$group == 1 else d$group != 1
    for (t in def$latent) {
      if (idm) {
        s1 <- d$e1_surv_time[gi] > t
        dead <- d$e2_surv_time[gi] <= t
        swt <- d$switched[gi] == 1 & !s1 & !dead
        vals <- list(p_state1 = s1, p_state2 = !s1 & !dead & !swt, p_state3 = swt,
                     p_state4 = dead, surv_os = !dead, surv_pfs = s1)
      } else {
        vals <- list(surv_os = d$surv_time[gi] > t)
      }
      for (q in names(vals)) {
        p <- mean(vals[[q]])
        ng <- sum(gi)
        rows[[length(rows) + 1L]] <- data.frame(
          source = "fastsurvival", scenario = name, hypothesis = "alternative", look_type = "latent",
          look = "", time = t, target_events = NA_real_, quantity = q,
          endpoint = if (q == "surv_os") "os" else if (q == "surv_pfs") "pfs" else "",
          group = g, mean = p, sd = NA_real_, se = sqrt(p * (1 - p) / ng), n = ng,
          stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

if (fs_available) {
  defs <- val_scenarios()
  out <- list()
  timing <- numeric(0)
  for (nm in c("V1", "V2", "V3", "V4", "V5")) {
    t0 <- proc.time()[["elapsed"]]
    out[[nm]] <- fs_simulate(nm, defs[[nm]])
    timing[nm] <- proc.time()[["elapsed"]] - t0
    message(sprintf("%s done in %.0f s", nm, timing[nm]))
  }
  fs_result <- do.call(rbind, out)
  rownames(fs_result) <- NULL
  utils::write.csv(fs_result, file.path(val_results_dir(), "sim_fastsurvival.csv"), row.names = FALSE)
  writeLines(c(paste("FastSurvival", fs_version), paste("time_s", names(timing), round(timing, 1))),
             file.path(val_results_dir(), "sim_fastsurvival_info.txt"))
} else {
  message("FastSurvival is not installed; V1-V5 are not simulated.")
}
