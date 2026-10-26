# Scenarios shared by the tests. They are defined identically in
# inst/validation/python/expected_values.py, which computes the expected
# values with an independent implementation (scipy expm and quad).

# A: exponential, hazard ratio 0.7, uniform accrual of 20 per month, 1:1
scn_a <- function() {
  ms_scenario(ms_model(death = log(2) / 12, dropout = 0.005),
              ms_model(death = 0.7 * log(2) / 12, dropout = 0.005),
              allocation = 1, accrual = ms_accrual(rate = 20))
}

# B: piecewise exponential with a delayed effect, piecewise accrual, 2:1
scn_b <- function() {
  ms_scenario(ms_model(time = c(0, 6), death = c(0.05, 0.04), dropout = 0.01),
              ms_model(time = c(0, 6), death = c(0.05, 0.025), dropout = 0.012),
              allocation = 2,
              accrual = ms_accrual(time = c(0, 3), rate = c(10, 30)))
}

# C: progression with switching in the control group, accrual by proportions
models_c <- function(clock = "reset") {
  list(control = ms_model(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01,
                          switch_prob = 0.3, time_post = c(0, 3),
                          death_prog = c(0.05, 0.06), switch_rate = 0.04,
                          death_switch = c(0.03, 0.035), dropout = 0.005,
                          clock = clock),
       treatment = ms_model(time = c(0, 4), prog = c(0.06, 0.05), death = 0.008,
                            switch_prob = 0, time_post = c(0, 3),
                            death_prog = c(0.045, 0.05), switch_rate = 0,
                            dropout = 0.005, clock = clock))
}

scn_c <- function(clock = "reset") {
  m <- models_c(clock)
  ms_scenario(m$control, m$treatment, allocation = 1,
              accrual = ms_accrual(time = c(0, 6), prop = c(0.3, 0.7),
                                   duration = 18))
}

# F: Fleischer-type model (exponential PFS and OS, death hazard unchanged by
# progression), 2:1, piecewise accrual, group-specific dropout; defined
# identically in inst/validation/python/expected_zeta.py
scn_f <- function() {
  b0 <- log(2) / 9
  b1 <- log(2) / 12.9
  ms_scenario(ms_model(prog = log(2) / 3 - b0, death = b0, death_prog = b0, dropout = 0.004),
              ms_model(prog = log(2) / 4.3 - b1, death = b1, death_prog = b1, dropout = 0.006),
              allocation = 2, accrual = ms_accrual(time = c(0, 4), rate = c(10, 25)))
}

# Stage 3: models and scenarios with the options of ms_model() and ms_scenario() of version
# 0.3.0, defined identically in inst/validation/python/expected_model5.py (model(),
# model_post() and scenario()). Further arguments of models_s3() are passed to ms_model() and
# replace those of the model.
models_s3 <- function(name, clock = "reset", ...) {
  ctl <- list(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01, switch_prob = 0.3,
              time_post = c(0, 3), death_prog = c(0.05, 0.06), switch_rate = 0.04,
              death_switch = c(0.03, 0.035), dropout = 0.005)
  trt <- list(time = c(0, 4), prog = c(0.06, 0.05), death = 0.008, switch_prob = 0,
              time_post = c(0, 3), death_prog = c(0.045, 0.05), switch_rate = 0,
              dropout = 0.005)
  args <- switch(name,
    P = c(ctl, list(switch_pre = c(0.03, 0.02), prog_switch = c(0.07, 0.06),
                    death_switch_pre = 0.007)),
    P2 = c(ctl, list(switch_pre = c(0.05, 0.01), prog_switch = c(0.04, 0.09),
                     death_switch_pre = c(0.012, 0.004))),
    Z = list(time = c(0, 4), death = 0.01, switch_pre = c(0.03, 0.02),
             prog_switch = c(0.07, 0.06), time_post = c(0, 3),
             death_switch = c(0.03, 0.035), dropout = 0.005),
    A06 = c(ctl[names(ctl) != "death_switch"], list(switch_af = 0.6)),
    A15 = c(ctl[names(ctl) != "death_switch"], list(switch_af = 1.5)),
    W = c(ctl, list(switch_clock = "switch")),
    T = c(trt, list(switch_pre = 0.01)),
    T0 = trt)
  do.call(ms_model, utils::modifyList(args, list(clock = clock, ...)))
}

# post-progression parameters with three intervals of time_post
model_post_s3 <- function(name) {
  args <- list(death = 0.01, time_post = c(0, 3, 8), death_prog = c(0.05, 0.08, 0.04),
               switch_rate = c(0.04, 0.06, 0.02))
  extra <- switch(name,
    switch = list(death_switch = c(0.03, 0.05, 0.02), switch_clock = "switch"),
    aft0.6 = list(switch_af = 0.6),
    aft1.5 = list(switch_af = 1.5))
  do.call(ms_model, c(args, extra))
}

scn_s3 <- function(name) {
  m <- switch(name,
    PF = list(models_s3("P", "forward"), models_s3("T", "forward"), "event"),
    PFC = list(models_s3("P", "forward", dropout_pfs = 0.02), models_s3("T", "forward"), "censor"),
    PR = list(models_s3("P"), models_s3("T"), "event"),
    PRC = list(models_s3("P"), models_s3("T"), "censor"),
    AR = list(models_s3("A06", dropout_pfs = c(0.02, 0.03)),
              models_s3("T0", dropout_pfs = 0.015), "event"),
    WF = list(models_s3("W", "forward", switch_pre = 0.02), models_s3("T", "forward"), "censor"))
  ms_scenario(m[[1L]], m[[2L]], allocation = 1,
              accrual = ms_accrual(time = c(0, 6), prop = c(0.3, 0.7), duration = 18),
              pfs_switch = m[[3L]])
}

# Stage 3, switching from a calendar time: scenario C with switching of the control group
# possible from the calendar time tc (switch_start), defined identically in
# inst/validation/python/expected_calendar.py (design()); null = TRUE gives the control model
# in both groups
scn_cal <- function(tc = 15, clock = "forward", null = FALSE) {
  m <- models_c(clock)
  ms_scenario(m$control, if (null) m$control else m$treatment, allocation = 1,
              accrual = ms_accrual(time = c(0, 6), prop = c(0.3, 0.7), duration = 18),
              switch_start = tc)
}

# Post-progression models of calendar_os.csv (models_post() of expected_calendar.py), named
# "<clock>_<variant>" with clock "reset" or "forward" and variant "same", "switch", "aft0.6"
# or "aft1.5"
model_post_cal <- function(name) {
  parts <- strsplit(name, "_", fixed = TRUE)[[1L]]
  args <- list(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01, switch_prob = 0.3,
               time_post = c(0, 3, 8), death_prog = c(0.05, 0.08, 0.04),
               switch_rate = c(0.04, 0.06, 0.02), dropout = 0.005, clock = parts[1L])
  extra <- switch(parts[2L],
    same = list(death_switch = c(0.03, 0.05, 0.02)),
    switch = list(death_switch = c(0.03, 0.05, 0.02), switch_clock = "switch"),
    aft0.6 = list(switch_af = 0.6),
    aft1.5 = list(switch_af = 1.5))
  do.call(ms_model, c(args, extra))
}

# Relative difference
rel_diff <- function(x, y) max(abs(x - y) / pmax(abs(y), 1e-300))
