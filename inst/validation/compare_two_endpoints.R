# Compare the simulation of PFS and OS together (python/simulate_two.py,
# scenarios V3, V4, V6, V7, V11 and V12) and the covariance estimator of
# Danzer et al. (2025) (python/danzer_check.py) with the package (stage 2,
# and stage 3 for V11, V12 and the Bonferroni procedure with recycling).
#
# Inputs: results/sim_python_two.csv and results/danzer_check.csv.
# Output: the rows of results/compare_two_endpoints.csv (sources "python_two"
# and "danzer"), judged in run_all.R.
#
# Judged (4se): the per-patient covariances of PFS at calendar look j and OS at
# calendar look k under the alternative hypothesis, estimated from 10 million
# simulated patients, against lr_cov(endpoint = c("pfs", "os")): the
# covariance of the scores divided by sqrt(V_P V_O) (variance =
# "alternative"), and the covariance and correlation of the influence
# functions of Z (variance = "studentized" and "null"). Judged (rel:1e-6): the
# same quantities in V3 computed deterministically by python/ref_zeta.py
# (Fleischer model). Judged (abs:1e-6): the efficacy bounds of the designs
# computed with python/ref_gsd.py. INFO: from the simulated trials, the
# correlation of Z_PFS and Z_OS at the calendar looks, and at the looks of
# the designs the look times, events, mean and standard deviation of Z, the
# correlations of the tested statistics and the rejection rates of each
# procedure (power and, under the global null hypothesis, familywise error
# rate with the bounds of the design); the correlation estimated as in Danzer
# et al. (2025) under the global null hypothesis.
#
# Designs (python/simulate_two.py defines the same designs independently):
#   V3          looks at 252 and 336 OS events, both endpoints tested at both
#               looks; OS first (O'Brien-Fleming type spending), then PFS
#               (Pocock type spending)
#   V4, V6, V7  look 1 at 330 PFS events (PFS and OS tested), look 2 at 250
#               OS events (OS tested); PFS first, then OS, both with
#               O'Brien-Fleming type spending
#   V11, V12    calendar looks at months 20 (PFS and OS tested) and 30 (OS
#               tested); PFS first, then OS, both with O'Brien-Fleming type
#               spending
# Procedures: fixed sequence with hierarchical = "overall" (fs_overall) and
# "stagewise" (fs_stagewise), each endpoint at level 0.025, Bonferroni
# with weights 0.8 (first endpoint) and 0.2 (bonferroni), and the same
# Bonferroni procedure with recycling of the level, transition c(1, 1)
# (bonferroni_recycle). V11 and V12 (stage 3, candidates 4 and 5) are
# scenario C (clock reset) with switching before progression in the control
# group (switch_pre 0.02, prog_switch 0.07 and 0.06, death_switch_pre
# 0.008); in V12 PFS is censored at that switch (pfs_switch = "censor") and
# both groups have a censoring of PFS only (dropout_pfs 0.01). The designs
# and V11 and V12 are defined here rather than in helpers.R, whose MD5
# decides whether the FastSurvival simulation is reused.

if (!exists("val_scenarios", mode = "function")) {
  source(file.path("inst", "validation", "helpers.R"))
}

val_scenarios_two3 <- function() {
  acc_c <- ms_accrual(time = c(0, 6), prop = c(0.3, 0.7), duration = 18)
  def <- function(cens) {
    dp <- if (cens) 0.01 else 0
    ctl <- ms_model(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01, switch_prob = 0.3,
                    time_post = c(0, 3), death_prog = c(0.05, 0.06), switch_rate = 0.04,
                    death_switch = c(0.03, 0.035), dropout = 0.005, clock = "reset",
                    switch_pre = 0.02, prog_switch = c(0.07, 0.06), death_switch_pre = 0.008,
                    dropout_pfs = dp)
    trt <- ms_model(time = c(0, 4), prog = c(0.06, 0.05), death = 0.008, time_post = c(0, 3),
                    death_prog = c(0.045, 0.05), dropout = 0.005, clock = "reset",
                    dropout_pfs = dp)
    list(scenario = ms_scenario(ctl, trt, accrual = acc_c,
                                pfs_switch = if (cens) "censor" else "event"),
         n = 500, cal = c(20, 30), events = NULL, latent = c(2, 5, 10, 20, 30),
         progression = TRUE)
  }
  list(V11 = def(FALSE), V12 = def(TRUE))
}

val_designs <- function() {
  des_c <- list(looks = ms_looks(c("pfs", "os"), fraction = c(1, 1),
                                 test = list(c("pfs", "os"), "os")),
                endpoints = c("pfs", "os"), sf_eff = c("OF", "OF"),
                events = c(pfs = 330, os = 250))
  des_cal <- list(looks = ms_looks("calendar", time = c(20, 30),
                                   test = list(c("pfs", "os"), "os")),
                  endpoints = c("pfs", "os"), sf_eff = c("OF", "OF"), events = NULL)
  list(V3 = list(looks = ms_looks("os", fraction = c(252, 336) / 336),
                 endpoints = c("os", "pfs"), sf_eff = c("OF", "Pocock"), events = 336),
       V4 = des_c, V6 = des_c, V7 = des_c, V11 = des_cal, V12 = des_cal)
}

val_test2 <- function(design, procedure, mean = "direct", variance = "null") {
  if (procedure %in% c("bonferroni", "bonferroni_recycle")) {
    return(ms_test(alpha = 0.025, endpoints = design$endpoints, sf_eff = design$sf_eff,
                   mean = mean, variance = variance, mtp = "bonferroni",
                   weights = c(0.8, 0.2),
                   transition = if (procedure == "bonferroni_recycle") c(1, 1) else c(0, 0)))
  }
  ms_test(alpha = 0.025, endpoints = design$endpoints, sf_eff = design$sf_eff, mean = mean,
          variance = variance, mtp = "fixed_sequence",
          hierarchical = sub("fs_", "", procedure, fixed = TRUE))
}

compare_two_endpoints <- function() {
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
    d$look[is.na(d$look)] <- ""
    d
  }
  sims <- rbind(read_sim("sim_python_two.csv"), read_sim("danzer_check.csv"))
  if (is.null(sims)) {
    return(NULL)
  }
  # rows of simulate_two.py are reported as "python_two", apart from the
  # stage 1 simulation (source "python")
  sims$source[sims$source == "python"] <- "python_two"
  defs <- c(val_scenarios(), val_scenarios_two3())
  des <- val_designs()
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
  pw_of <- function(nm, hyp, proc, mean, variance) {
    memo(paste("pw", nm, hyp, proc, mean, variance), function() {
      ms_power(scn_of(nm, hyp), n = defs[[nm]]$n, looks = des[[nm]]$looks,
               test = val_test2(des[[nm]], proc, mean, variance), events = des[[nm]]$events)
    })
  }
  cov_of <- function(nm, hyp, variance) {
    memo(paste("cov", nm, hyp, variance), function() {
      lr_cov(scn_of(nm, hyp), n = defs[[nm]]$n, calendar_time = defs[[nm]]$cal,
             endpoint = c("pfs", "os"), variance = variance)$sigma
    })
  }
  # row of the looks table of ms_power() for an endpoint and a look
  item <- function(lk, ep, look) {
    which(lk$endpoint == ep & lk$look == look)
  }
  # rejection probabilities under the global null hypothesis with the bounds
  # of the design (computed under the alternative)
  null_reject <- function(nm, proc) {
    memo(paste("null", nm, proc), function() {
      pa <- pw_of(nm, "alternative", proc, "direct", "null")
      p0 <- pw_of(nm, "null", proc, "direct", "null")
      eps <- des[[nm]]$endpoints
      lk <- pa$looks
      up <- lapply(eps, function(e) lk$z_upper[lk$endpoint == e])
      an <- lapply(eps, function(e) lk$look[lk$endpoint == e])
      upr <- if (!is.null(lk$z_upper_recycled)) {
        lapply(eps, function(e) lk$z_upper_recycled[lk$endpoint == e])
      }
      ts <- val_test2(des[[nm]], proc)
      r <- gs_prob_mtp(up, mean = 0, sigma = p0$sigma, analysis = an, mtp = ts$mtp,
                       hierarchical = ts$hierarchical, upper_recycled = upr)$reject
      out <- c(r[["first"]], r[["second"]], r[["both"]], r[["any"]])
      names(out) <- c(eps, "both", "any")
      out
    })
  }
  vars_of <- function(hyp) {
    if (hyp == "null") "null" else c("null", "alternative", "studentized")
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
  for (i in seq_len(nrow(sims))) {
    r <- sims[i, ]
    nm <- r$scenario
    hyp <- r$hypothesis
    nc <- length(defs[[nm]]$cal)
    jk <- suppressWarnings(as.integer(strsplit(r$look, "-", fixed = TRUE)[[1L]]))
    if (r$look_type == "calendar") {
      # correlation of Z_PFS at look j and Z_OS at look k
      for (v in vars_of(hyp)) {
        s <- cov_of(nm, hyp, v)
        emit(r, paste0("variance=", v), "info", r$mean, r$se,
             s[jk[1L], nc + jk[2L]] / sqrt(s[jk[1L], jk[1L]] * s[nc + jk[2L], nc + jk[2L]]))
      }
      next
    }
    if (r$look_type %in% c("perpatient", "reference")) {
      v <- c(cov_alt = "alternative", cov_stud = "studentized", corr_stud = "null")[[r$quantity]]
      crit <- if (r$look_type == "perpatient") "4se" else "rel:1e-6"
      emit(r, paste0("variance=", v), crit, r$mean, r$se, cov_of(nm, "alternative", v)[jk[1L], nc + jk[2L]])
      next
    }
    # looks of the design
    base <- pw_of(nm, hyp, "fs_overall", "direct", "null")
    lk <- base$looks
    k <- jk[1L]
    if (r$quantity == "time") {
      emit(r, "", "info", r$mean, r$se, lk$time[lk$look == k][1L])
    } else if (r$quantity == "events") {
      tt <- lk$time[lk$look == k][1L]
      ev <- memo(paste("ev", nm, hyp, r$endpoint, k), function() {
        ms_events(scn_of(nm, hyp), n = defs[[nm]]$n, calendar_time = tt, endpoint = r$endpoint)
      })
      emit(r, "", "info", r$mean, r$se, ev[[paste0("events_", r$group)]])
    } else if (r$quantity == "z") {
      ii <- item(lk, r$endpoint, k)
      for (mn in c("direct", "schoenfeld")) {
        emit(r, paste0("mean=", mn), "info", r$mean, r$se,
             pw_of(nm, hyp, "fs_overall", mn, "null")$looks$z_mean[ii])
      }
      for (v in vars_of(hyp)) {
        emit(r, paste0("sd;variance=", v), "info", r$sd, r$sd / sqrt(2 * (r$n - 1)),
             pw_of(nm, hyp, "fs_overall", "direct", v)$looks$z_sd[ii])
      }
    } else if (r$quantity == "z_corr") {
      ee <- strsplit(r$endpoint, "-", fixed = TRUE)[[1L]]
      i1 <- item(lk, ee[1L], jk[1L])
      i2 <- item(lk, ee[2L], jk[2L])
      for (v in vars_of(hyp)) {
        s <- pw_of(nm, hyp, "fs_overall", "direct", v)$sigma
        emit(r, paste0("variance=", v), "info", r$mean, r$se, s[i1, i2] / sqrt(s[i1, i1] * s[i2, i2]))
      }
    } else if (r$quantity == "reject") {
      if (hyp == "null") {
        emit(r, "global null, bounds of the design", "info", r$mean, r$se,
             null_reject(nm, r$group)[[r$endpoint]])
        if (r$endpoint == "any") {
          emit(r, "type1 of ms_power", "info", r$mean, r$se,
               pw_of(nm, "alternative", r$group, "direct", "null")$type1[["nonbinding"]])
        }
      } else {
        for (v in vars_of(hyp)) {
          emit(r, paste0("variance=", v), "info", r$mean, r$se,
               pw_of(nm, hyp, r$group, "direct", v)$power[[r$endpoint]])
        }
      }
    } else if (r$quantity == "bound") {
      proc <- switch(r$group, fs = "fs_overall", bonferroni_recycled = "bonferroni_recycle",
                     "bonferroni")
      la <- pw_of(nm, "alternative", proc, "direct", "null")$looks
      col <- if (r$group == "bonferroni_recycled") "z_upper_recycled" else "z_upper"
      emit(r, "", "abs:1e-6", r$mean, 0, la[[col]][item(la, r$endpoint, k)])
    } else if (r$quantity == "not_reached") {
      emit(r, "", "info", r$mean, NA_real_, 0)
    }
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
