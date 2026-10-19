# Compare the simulation summaries with the package (Phase 3 plan, sections
# 4.1, 4.2 and 12.4).
#
# Inputs: results/sim_python.csv (python/simulate_ms.py, V6-V8) and
# results/sim_fastsurvival.csv (validate_fastsurvival.R, V1-V5).
# Output: results/compare_simulation.csv.
#
# Judged (criterion 4se or 4se_p): calendar-time group-wise OS and PFS events
# and dropouts under both hypotheses, and the latent state occupation
# probabilities and OS and PFS survival without dropout. Judged (abs:1e-6):
# the O'Brien-Fleming bounds of the event-driven looks computed in Python
# (ref_gsd.py). INFO: log-rank Z (mean, standard deviation, correlation
# between looks), the event-driven look times and events, the rejection rates
# (power and type I error), and the dropouts at the event-driven looks
# (compared with the expected dropouts at the expected look time).

if (!exists("val_scenarios", mode = "function")) {
  source(file.path("inst", "validation", "helpers.R"))
}

compare_simulation <- function() {
  read_sim <- function(f) {
    p <- file.path(val_results_dir(), f)
    if (!file.exists(p)) {
      return(NULL)
    }
    d <- utils::read.csv(p, stringsAsFactors = FALSE,
                         colClasses = c(look = "character", group = "character",
                                        endpoint = "character"))
    d$group[is.na(d$group)] <- ""
    d$endpoint[is.na(d$endpoint)] <- ""
    d
  }
  sims <- rbind(read_sim("sim_python.csv"), read_sim("sim_fastsurvival.csv"))
  if (is.null(sims)) {
    stop("No simulation results in ", val_results_dir(), call. = FALSE)
  }
  defs <- val_scenarios()
  cache <- new.env()
  memo <- function(key, f) {
    if (is.null(cache[[key]])) {
      cache[[key]] <- f()
    }
    cache[[key]]
  }
  scn_of <- function(nm, hyp) {
    s <- defs[[nm]]$scenario
    if (hyp == "null") ms_null(s) else s
  }
  pw_of <- function(nm, hyp, ep, mean, variance) {
    memo(paste("pw", nm, hyp, ep, mean, variance), function() {
      d <- defs[[nm]]
      ms_power(scn_of(nm, hyp), n = d$n,
               looks = ms_looks("os", fraction = d$events / d$events[length(d$events)]),
               test = ms_test(alpha = 0.025, endpoints = ep, sf_eff = "OF", mean = mean,
                              variance = variance),
               events = d$events[length(d$events)])
    })
  }
  cov_of <- function(nm, hyp, ep, times, variance) {
    memo(paste("cov", nm, hyp, ep, paste(signif(times, 10), collapse = ","), variance), function() {
      lr_cov(scn_of(nm, hyp), n = defs[[nm]]$n, calendar_time = times, endpoint = ep,
             variance = variance)
    })
  }
  mean_of <- function(nm, hyp, ep, times, mean) {
    memo(paste("mean", nm, hyp, ep, paste(signif(times, 10), collapse = ","), mean), function() {
      lr_moments(scn_of(nm, hyp), n = defs[[nm]]$n, calendar_time = times, endpoint = ep,
                 mean = mean)$z_mean
    })
  }
  out <- list()
  emit <- function(r, variant, criterion, simulated, se, package, n = r$n) {
    out[[length(out) + 1L]] <<- data.frame(
      id = paste(r$source, r$scenario, r$hypothesis, r$look_type, r$look,
                 ifelse(is.na(r$time), "", format(r$time)), r$quantity, r$endpoint, r$group,
                 variant, sep = "|"),
      source = r$source, scenario = r$scenario, hypothesis = r$hypothesis,
      look_type = r$look_type, look = r$look, time = r$time, quantity = r$quantity,
      endpoint = r$endpoint, group = r$group, variant = variant, criterion = criterion,
      simulated = simulated, se = se, n = n, package = package, stringsAsFactors = FALSE)
  }
  latent_name <- c(p_state1 = "progression_free", p_state2 = "progressed",
                   p_state3 = "switched", p_state4 = "dead")
  for (i in seq_len(nrow(sims))) {
    r <- sims[i, ]
    d <- defs[[r$scenario]]
    nm <- r$scenario
    hyp <- r$hypothesis
    ep <- r$endpoint
    k <- suppressWarnings(as.integer(r$look))
    if (r$look_type == "latent") {
      m <- scn_of(nm, hyp)[[r$group]]
      pk <- if (r$quantity %in% names(latent_name)) {
        ms_prob(m, r$time)[[latent_name[[r$quantity]]]]
      } else {
        ms_surv(m, r$time, endpoint = ep)$surv
      }
      emit(r, "", "4se_p", r$mean, r$se, pk)
      next
    }
    if (r$look_type == "calendar") {
      if (r$quantity %in% c("events", "dropouts")) {
        ev <- memo(paste("ev", nm, hyp, ep, r$time), function() {
          ms_events(scn_of(nm, hyp), n = d$n, calendar_time = r$time, endpoint = ep)
        })
        emit(r, "", "4se", r$mean, r$se, ev[[paste0(r$quantity, "_", r$group)]])
      } else if (r$quantity == "z") {
        for (mn in c("direct", "schoenfeld")) {
          emit(r, paste0("mean=", mn), "info", r$mean, r$se, mean_of(nm, hyp, ep, d$cal, mn)[k])
        }
        for (v in c("null", "alternative")) {
          cv <- cov_of(nm, hyp, ep, d$cal, v)
          emit(r, paste0("sd;variance=", v), "info", r$sd, r$sd / sqrt(2 * (r$n - 1)),
               sqrt(cv$sigma[k, k]))
        }
      } else if (r$quantity == "z_corr") {
        for (v in c("null", "alternative")) {
          s <- cov_of(nm, hyp, ep, d$cal, v)$sigma
          emit(r, paste0("variance=", v), "info", r$mean, r$se, s[1, 2] / sqrt(s[1, 1] * s[2, 2]))
        }
      }
      next
    }
    # event-driven looks
    if (r$quantity == "time") {
      emit(r, "", "info", r$mean, r$se, pw_of(nm, hyp, "os", "direct", "null")$looks$time[k])
    } else if (r$quantity == "events") {
      lk <- pw_of(nm, hyp, ep, "direct", "null")$looks
      emit(r, "", "info", r$mean, r$se, lk[[paste0("events_", r$group)]][k])
    } else if (r$quantity == "dropouts") {
      tt <- pw_of(nm, hyp, "os", "direct", "null")$looks$time[k]
      ev <- ms_events(scn_of(nm, hyp), n = d$n, calendar_time = tt, endpoint = ep)
      emit(r, "at the expected look time", "info", r$mean, r$se,
           ev[[paste0("dropouts_", r$group)]])
    } else if (r$quantity == "z") {
      for (mn in c("direct", "schoenfeld")) {
        emit(r, paste0("mean=", mn), "info", r$mean, r$se,
             pw_of(nm, hyp, ep, mn, "null")$looks$z_mean[k])
      }
      for (v in c("null", "alternative")) {
        emit(r, paste0("sd;variance=", v), "info", r$sd, r$sd / sqrt(2 * (r$n - 1)),
             pw_of(nm, hyp, ep, "direct", v)$looks$z_sd[k])
      }
    } else if (r$quantity == "z_corr") {
      tt <- pw_of(nm, hyp, "os", "direct", "null")$looks$time
      for (v in c("null", "alternative")) {
        s <- cov_of(nm, hyp, ep, tt, v)$sigma
        emit(r, paste0("variance=", v), "info", r$mean, r$se, s[1, 2] / sqrt(s[1, 1] * s[2, 2]))
      }
    } else if (r$quantity %in% c("reject_any", "reject_look1")) {
      for (vv in list(c("direct", "null"), c("direct", "alternative"), c("schoenfeld", "null"))) {
        p <- pw_of(nm, hyp, "os", vv[1], vv[2])
        pk <- if (r$quantity == "reject_any") p$power else p$looks$cum_power[1]
        emit(r, paste0("mean=", vv[1], ";variance=", vv[2]), "info", r$mean, r$se, pk)
      }
    } else if (r$quantity == "bound") {
      emit(r, "", "abs:1e-6", r$mean, 0, pw_of(nm, hyp, "os", "direct", "null")$looks$z_upper[k])
    } else if (r$quantity == "not_reached") {
      emit(r, "", "info", r$mean, NA_real_, 0)
    }
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
