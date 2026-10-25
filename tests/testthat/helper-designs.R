# Designs with two endpoints or with looks driven by the events of different
# endpoints, defined identically in inst/validation/python/expected_design.py
# (scenario F of helper-scenarios.R, n = 400). Design D5 (stage 3) is design D2
# with recycling of the level and has no values in design2.csv.

# Scenario F with the same OS ("os") or the same PFS ("pfs") in both groups
scn_f_null <- function(which) {
  b0 <- log(2) / 9
  b1 <- log(2) / 12.9
  ctl <- ms_model(prog = log(2) / 3 - b0, death = b0, death_prog = b0, dropout = 0.004)
  trt <- switch(which,
    os = ms_model(prog = log(2) / 4.3 - b0, death = b0, death_prog = b0, dropout = 0.006),
    pfs = ms_model(prog = log(2) / 3 - b1, death = b1, death_prog = b1, dropout = 0.006))
  ms_scenario(ctl, trt, allocation = 2, accrual = ms_accrual(time = c(0, 4), rate = c(10, 25)))
}

design_def <- function(name) {
  switch(name,
    D1 = list(looks = ms_looks(c("pfs", "pfs", "os"), fraction = c(0.6, 1, 1),
                               test = list("pfs", c("pfs", "os"), "os")),
              test = ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence",
                             variance = "studentized"),
              events = c(pfs = 250, os = 220), followup = NULL,
              null = list(os = scn_f_null("os"))),
    D2 = list(looks = ms_looks(c("calendar", "os"), fraction = c(NA, 1), time = c(15, NA)),
              test = ms_test(endpoints = c("pfs", "os"), sf_eff = c("OF", "Pocock"),
                             futility = "beta", sf_fut = "HSD", param_fut = -2, beta = 0.1,
                             mtp = "bonferroni", weights = c(0.8, 0.2)),
              events = NULL, followup = 12, null = list(pfs = scn_f_null("pfs"))),
    D5 = list(looks = ms_looks(c("calendar", "os"), fraction = c(NA, 1), time = c(15, NA)),
              test = ms_test(endpoints = c("pfs", "os"), sf_eff = c("OF", "Pocock"),
                             futility = "beta", sf_fut = "HSD", param_fut = -2, beta = 0.1,
                             mtp = "bonferroni", weights = c(0.8, 0.2), transition = c(1, 1)),
              events = NULL, followup = 12, null = list(pfs = scn_f_null("pfs"))),
    D3 = list(looks = ms_looks("os", fraction = c(0.5, 1)),
              test = ms_test(endpoints = c("pfs", "os"), mtp = "fixed_sequence",
                             hierarchical = "stagewise", variance = "alternative"),
              events = 220, followup = NULL, null = NULL),
    D4 = list(looks = ms_looks(c("pfs", "os"), fraction = c(1, 1)),
              test = ms_test(variance = "studentized"),
              events = c(pfs = 250, os = 220), followup = NULL, null = NULL))
}

run_design <- function(name) {
  d <- design_def(name)
  ms_power(scn_f(), n = 400, looks = d$looks, test = d$test, events = d$events,
           followup = d$followup, null_scenarios = d$null)
}

# Expected value of a quantity of design 'name' (design2.csv); endpoint and
# look select the rows (NA: not used)
design_expected <- function(name, quantity, endpoint = NA, look = NA) {
  ex <- expected_csv("design2.csv")
  e <- ex[ex$config == name & ex$quantity == quantity, ]
  if (!is.na(endpoint)) {
    e <- e[e$endpoint == endpoint, ]
  }
  if (!is.na(look)) {
    e <- e[e$look == look, ]
  }
  e
}

# Values of a column of res$looks in the order of the rows of the expected
# values e (matched by endpoint and look)
looks_value <- function(res, e, column) {
  lk <- res$looks
  ep <- if ("endpoint" %in% names(lk)) lk$endpoint else e$endpoint[1L]
  key <- paste(ep, lk$look)
  lk[[column]][match(paste(e$endpoint, e$look), key)]
}
