# Compare the simulation of switching from a calendar time
# (python/simulate_calendar.py, scenarios V13-V15; stage 3, candidate 3) with
# the package.
#
# Input: results/sim_python_calendar.csv.
# Output: the rows of results/compare_calendar.csv (source "python_calendar"),
# judged in run_all.R.
#
# Judged (4se): per patient, the OS events by group and in total, the mean of
# the log-rank score and its information at calendar months 24 and 36, from
# 100 trials of 20000 patients, against ms_events() and lr_moments() divided
# by the number of patients. Judged (abs:1e-6): the O'Brien-Fleming bounds of
# the looks at 150 and 250 OS events computed in Python (ref_gsd.py). INFO:
# Z of OS at the calendar looks (mean against lr_moments(), standard
# deviation against variance "null", "alternative" and "studentized") and its
# correlation between the looks (lr_cov() with variance "null" and
# "studentized"), and at the looks at 150 and 250 OS events the look times, Z
# (mean and standard deviation), its correlation and the rejection rates of
# the O'Brien-Fleming design (ms_power()).
#
# Scenarios (python/simulate_calendar.py defines the same scenarios
# independently): scenario C of the tests (helpers.R, val_model_c()) with
# switching possible from calendar month 15 (switch_start = 15)
#   V13  clock forward
#   V14  clock reset
#   V15  clock forward, accelerated failure time effect of the switch
#        (switch_af = 0.6 instead of death_switch)
# The scenarios are defined here rather than in helpers.R, whose MD5 decides
# whether the FastSurvival simulation is reused.

if (!exists("val_scenarios", mode = "function")) {
  source(file.path("inst", "validation", "helpers.R"))
}

val_scenarios_calendar <- function() {
  acc_c <- ms_accrual(time = c(0, 6), prop = c(0.3, 0.7), duration = 18)
  models <- function(clock, af = NULL) {
    list(control = ms_model(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01,
                            switch_prob = 0.3, time_post = c(0, 3),
                            death_prog = c(0.05, 0.06), switch_rate = 0.04,
                            death_switch = if (is.null(af)) c(0.03, 0.035),
                            dropout = 0.005, clock = clock, switch_af = af),
         treatment = ms_model(time = c(0, 4), prog = c(0.06, 0.05), death = 0.008,
                              switch_prob = 0, time_post = c(0, 3),
                              death_prog = c(0.045, 0.05), switch_rate = 0,
                              dropout = 0.005, clock = clock))
  }
  def <- function(m) {
    list(scenario = ms_scenario(m$control, m$treatment, accrual = acc_c, switch_start = 15),
         n = 500, cal = c(24, 36), events = c(150, 250))
  }
  list(V13 = def(models("forward")), V14 = def(models("reset")),
       V15 = def(models("forward", af = 0.6)))
}

compare_calendar <- function() {
  p <- file.path(val_results_dir(), "sim_python_calendar.csv")
  if (!file.exists(p)) {
    return(NULL)
  }
  sims <- utils::read.csv(p, stringsAsFactors = FALSE,
                          colClasses = c(look = "character", group = "character",
                                         endpoint = "character"))
  sims$group[is.na(sims$group)] <- ""
  sims$look[is.na(sims$look)] <- ""
  sims$source <- "python_calendar"
  defs <- val_scenarios_calendar()
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
  # package values that may be outside the scope of the calendar-time
  # switching are reported as NA with the error message
  try_pk <- function(f) {
    tryCatch(list(value = f(), msg = ""), error = function(e) list(value = NULL, msg = conditionMessage(e)))
  }
  num1 <- function(x) if (is.numeric(x) && length(x) == 1L) x else NA_real_
  pw_of <- function(nm, hyp, variance) {
    memo(paste("pw", nm, hyp, variance), function() {
      d <- defs[[nm]]
      try_pk(function() {
        ms_power(scn_of(nm, hyp), n = d$n,
                 looks = ms_looks("os", fraction = d$events / d$events[length(d$events)]),
                 test = ms_test(alpha = 0.025, endpoints = "os", sf_eff = "OF", variance = variance),
                 events = d$events[length(d$events)])
      })
    })
  }
  mom_of <- function(nm, hyp, n, variance) {
    memo(paste("mom", nm, hyp, n, variance), function() {
      try_pk(function() {
        lr_moments(scn_of(nm, hyp), n = n, calendar_time = defs[[nm]]$cal, endpoint = "os",
                   variance = variance)
      })
    })
  }
  cov_of <- function(nm, hyp, times, variance, trigger = NULL) {
    memo(paste("cov", nm, hyp, paste(signif(times, 10), collapse = ","), variance,
               paste(trigger, collapse = ",")), function() {
      try_pk(function() {
        lr_cov(scn_of(nm, hyp), n = defs[[nm]]$n, calendar_time = times, endpoint = "os",
               variance = variance, trigger = trigger)$sigma
      })
    })
  }
  corr12 <- function(s) {
    if (is.matrix(s$value)) s$value[1, 2] / sqrt(s$value[1, 1] * s$value[2, 2]) else NA_real_
  }
  out <- list()
  emit <- function(r, variant, criterion, simulated, se, package, n = r$n, note = "") {
    out[[length(out) + 1L]] <<- data.frame(
      id = paste(r$source, r$scenario, r$hypothesis, r$look_type, r$look,
                 ifelse(is.na(r$time), "", format(r$time)), r$quantity, r$endpoint, r$group,
                 variant, sep = "|"),
      source = r$source, scenario = r$scenario, hypothesis = r$hypothesis,
      look_type = r$look_type, look = r$look, time = r$time, quantity = r$quantity,
      endpoint = r$endpoint, group = r$group, variant = variant, criterion = criterion,
      simulated = simulated, se = se, n = n, package = package, note = note,
      stringsAsFactors = FALSE)
  }
  big_n <- 20000
  for (i in seq_len(nrow(sims))) {
    r <- sims[i, ]
    nm <- r$scenario
    hyp <- r$hypothesis
    k <- suppressWarnings(as.integer(r$look))
    cal <- defs[[nm]]$cal
    if (r$look_type == "perpatient") {
      j <- match(r$time, cal)
      if (r$quantity == "events") {
        ev <- memo(paste("ev", nm, hyp), function() {
          ms_events(scn_of(nm, hyp), n = big_n, calendar_time = cal, endpoint = "os")
        })
        col <- if (nzchar(r$group)) paste0("events_", r$group) else "events"
        emit(r, "", "4se", r$mean, r$se, ev[[col]][j] / big_n)
      } else {
        m <- mom_of(nm, hyp, big_n, "null")
        col <- if (r$quantity == "score") "score_mean" else "info"
        emit(r, "", "4se", r$mean, r$se, num1(m$value[[col]][j]) / big_n, note = m$msg)
      }
      next
    }
    if (r$look_type == "calendar") {
      if (r$quantity == "z") {
        m <- mom_of(nm, hyp, defs[[nm]]$n, "null")
        emit(r, "mean=direct", "info", r$mean, r$se, num1(m$value$z_mean[k]), note = m$msg)
        for (v in c("null", "alternative", "studentized")) {
          m <- mom_of(nm, hyp, defs[[nm]]$n, v)
          emit(r, paste0("sd;variance=", v), "info", r$sd, r$sd / sqrt(2 * (r$n - 1)),
               sqrt(num1(m$value$z_var[k])), note = m$msg)
        }
      } else if (r$quantity == "z_corr") {
        for (v in c("null", "studentized")) {
          s <- cov_of(nm, hyp, cal, v)
          emit(r, paste0("variance=", v), "info", r$mean, r$se, corr12(s), note = s$msg)
        }
      }
      next
    }
    # looks at 150 and 250 OS events
    base <- pw_of(nm, hyp, "null")
    lk <- base$value$looks
    if (r$quantity == "time") {
      emit(r, "", "info", r$mean, r$se, num1(lk$time[k]), note = base$msg)
    } else if (r$quantity == "z") {
      emit(r, "mean=direct", "info", r$mean, r$se, num1(lk$z_mean[k]), note = base$msg)
      for (v in c("null", "studentized")) {
        p <- pw_of(nm, hyp, v)
        emit(r, paste0("sd;variance=", v), "info", r$sd, r$sd / sqrt(2 * (r$n - 1)),
             num1(p$value$looks$z_sd[k]), note = p$msg)
      }
    } else if (r$quantity == "z_corr") {
      # correlation at the expected look times, which are driven by OS events
      for (v in c("null", "studentized")) {
        s <- if (length(lk$time) == 2L) cov_of(nm, hyp, lk$time, v, "os") else base
        emit(r, paste0("variance=", v), "info", r$mean, r$se, corr12(s), note = s$msg)
      }
    } else if (r$quantity %in% c("reject_any", "reject_look1")) {
      for (v in c("null", "studentized")) {
        p <- pw_of(nm, hyp, v)
        pk <- if (r$quantity == "reject_any") p$value$power else p$value$looks$cum_power[1]
        emit(r, paste0("variance=", v), "info", r$mean, r$se, num1(pk), note = p$msg)
      }
    } else if (r$quantity == "bound") {
      emit(r, "", "abs:1e-6", r$mean, 0, num1(lk$z_upper[k]), note = base$msg)
    } else if (r$quantity == "not_reached") {
      emit(r, "", "info", r$mean, NA_real_, 0)
    }
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
