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

# Relative difference
rel_diff <- function(x, y) max(abs(x - y) / pmax(abs(y), 1e-300))
