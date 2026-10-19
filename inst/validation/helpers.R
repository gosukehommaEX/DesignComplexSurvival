# Helpers for the validation scripts in inst/validation.
#
# The scripts are run from the package root after devtools::load_all() (or
# with the package installed). Scenario definitions V1-V8 (Phase 3 plan,
# sections 4.1, 4.2 and 12) are given here once and used by
# validate_fastsurvival.R and compare_simulation.R; the Python simulation
# (python/simulate_ms.py) defines V6-V8 independently.

if (!exists("ms_size", mode = "function")) {
  library(DesignComplexSurvival)
}

val_dir <- function() {
  d <- file.path("inst", "validation")
  if (dir.exists(d)) {
    return(d)
  }
  d <- system.file("validation", package = "DesignComplexSurvival")
  if (!nzchar(d)) {
    stop("Cannot find inst/validation. Set the working directory to the package root.",
         call. = FALSE)
  }
  d
}

val_results_dir <- function() {
  d <- file.path(val_dir(), "results")
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
  d
}

# Scenarios -----------------------------------------------------------------
#
# Each scenario: scenario (ms_scenario), n, calendar looks (cal), OS event
# targets of the event-driven looks (events, NULL for none), times since
# randomization for the latent quantities (latent), and whether the model
# has progression (two endpoints).

val_model_c <- function(clock = "reset", switch_rate = 0.04) {
  list(control = ms_model(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01,
                          switch_prob = 0.3, time_post = c(0, 3),
                          death_prog = c(0.05, 0.06), switch_rate = switch_rate,
                          death_switch = c(0.03, 0.035), dropout = 0.005,
                          clock = clock),
       treatment = ms_model(time = c(0, 4), prog = c(0.06, 0.05), death = 0.008,
                            switch_prob = 0, time_post = c(0, 3),
                            death_prog = c(0.045, 0.05), switch_rate = 0,
                            dropout = 0.005, clock = clock))
}

val_scenarios <- function() {
  acc_c <- ms_accrual(time = c(0, 6), prop = c(0.3, 0.7), duration = 18)
  def_c <- function(m) {
    list(scenario = ms_scenario(m$control, m$treatment, accrual = acc_c), n = 500,
         cal = c(20, 30), events = c(150, 250), latent = c(2, 5, 10, 20, 30),
         progression = TRUE)
  }
  # V8: Luo et al. (2019), Example 1, delayed effect with treatment
  # discontinuation (dilution, Markov)
  tm <- c(0, 6, 12, 24)
  lam0 <- c(0.09, 0.09, 0.07, 0.05) / 12
  drop <- c(0, 0, 0.015, 0.015) / 12
  disc <- c(0.04, 0.04, 0.03, 0.03) / 12
  hr_d <- c(1, 0.75, 0.75, 0.75)
  # V5: Barthel et al. (2006), Table V, row 9 (years)
  lc <- log(2)
  he <- c(0.6, 0.8) * lc
  dh <- function(x) -log(1 - x) / 4
  dr3 <- -log(0.95) / 12
  list(
    V1 = list(scenario = ms_scenario(ms_model(death = log(2) / 12, dropout = 0.005),
                                     ms_model(death = 0.7 * log(2) / 12, dropout = 0.005),
                                     accrual = ms_accrual(rate = 20)),
              n = 400, cal = c(24, 36), events = c(200, 300), latent = c(6, 12, 24, 36),
              progression = FALSE),
    V2 = list(scenario = ms_scenario(ms_model(time = c(0, 6), death = c(0.05, 0.04), dropout = 0.01),
                                     ms_model(time = c(0, 6), death = c(0.05, 0.025), dropout = 0.012),
                                     allocation = 2,
                                     accrual = ms_accrual(time = c(0, 3), rate = c(10, 30))),
              n = 300, cal = c(15, 25), events = c(75, 150), latent = c(3, 6, 12, 24),
              progression = FALSE),
    V3 = list(scenario = ms_scenario(ms_fleischer(median_pfs = 3, median_os = 9, dropout = dr3),
                                     ms_fleischer(median_pfs = 4.3, median_os = 12.9, dropout = dr3),
                                     accrual = ms_accrual(time = 0, prop = 1, duration = 23)),
              n = 482, cal = c(27, 35), events = c(252, 336), latent = c(3, 6, 12, 24),
              progression = TRUE),
    V4 = def_c(val_model_c("reset", switch_rate = 0)),
    V5 = list(scenario = ms_scenario(
                ms_model(death = lc, prog = dh(0.1), death_prog = he, time_post = c(0, 2),
                         dropout = dh(0.2), clock = "reset"),
                ms_model(time = c(0, 2), death = he, prog = dh(0.1), death_prog = lc,
                         dropout = dh(0.2), clock = "forward"),
                accrual = ms_accrual(time = 0, prop = 1, duration = 2)),
              n = 292, cal = c(3, 4), events = NULL, latent = c(1, 2, 3, 4),
              progression = TRUE),
    V6 = def_c(val_model_c("reset")),
    V7 = def_c(val_model_c("forward")),
    V8 = list(scenario = ms_scenario(
                ms_model(time = tm, death = lam0, dropout = drop),
                ms_model(time = tm, death = hr_d * lam0, dropout = drop, prog = disc,
                         death_prog = lam0, time_post = tm, clock = "forward"),
                accrual = ms_accrual(time = 0:13,
                                     rate = c(100, 200, 300, 300, rep(400, 8), 300, 200))),
              n = 4600, cal = c(30, 40), events = c(361, 722), latent = c(6, 12, 24, 36),
              progression = TRUE)
  )
}

# Judgement -------------------------------------------------------------------
#
# Criteria:
#   4se     |simulated - package| <= 4 se (se of the simulation mean)
#   4se_p   as 4se, with se = sqrt(p (1 - p) / n) from the package value p
#           (proportions of simulated patients)
#   abs:tol |simulated - package| <= tol
#   rel:tol |simulated / package - 1| <= tol
#   info    not judged

val_se <- function(criterion, se, package, n) {
  if (criterion == "4se_p") sqrt(pmax(package * (1 - package), 0) / n) else se
}

val_pass <- function(criterion, simulated, package, se) {
  if (criterion == "info") {
    return(NA)
  }
  if (is.na(simulated) || is.na(package)) {
    return(FALSE)
  }
  eps <- 1e-9 * max(1, abs(package))
  if (criterion %in% c("4se", "4se_p")) {
    return(abs(simulated - package) <= 4 * se + eps)
  }
  parts <- strsplit(criterion, ":", fixed = TRUE)[[1L]]
  tol <- as.numeric(parts[2L])
  if (parts[1L] == "abs") {
    return(abs(simulated - package) <= tol + eps)
  }
  if (parts[1L] == "rel") {
    return(abs(simulated / package - 1) <= tol + 1e-12)
  }
  stop("Unknown criterion: ", criterion, call. = FALSE)
}

# Judge a data frame with columns id, criterion, simulated, package, se and n.
# With perturb = TRUE the package values are multiplied by 1.1. Rows for which
# a 10% change of the package value exceeds the tolerance plus the observed
# difference |simulated - package| must FAIL after the change; they are marked
# in the column sensitive.
val_judge <- function(d, explained = character(0), perturb = FALSE) {
  f <- if (perturb) 1.1 else 1
  pk <- d$package * f
  se <- vapply(seq_len(nrow(d)), function(i) val_se(d$criterion[i], d$se[i], d$package[i], d$n[i]),
               numeric(1))
  ok <- vapply(seq_len(nrow(d)), function(i) val_pass(d$criterion[i], d$simulated[i], pk[i], se[i]),
               logical(1))
  j <- ifelse(d$criterion == "info", "INFO", ifelse(ok, "PASS", "FAIL"))
  j[j == "FAIL" & d$id %in% explained] <- "EXPLAINED"
  tol <- vapply(seq_len(nrow(d)), function(i) {
    cr <- d$criterion[i]
    if (cr %in% c("4se", "4se_p")) return(4 * se[i])
    if (cr == "info") return(NA_real_)
    p <- strsplit(cr, ":", fixed = TRUE)[[1L]]
    # rel: the tolerance after the change is tol * 1.1 * |package|
    if (p[1L] == "abs") as.numeric(p[2L]) else as.numeric(p[2L]) * 1.1 * abs(d$package[i])
  }, numeric(1))
  dev <- abs(d$simulated - d$package)
  margin <- 1e-9 * pmax(1, 1.1 * abs(d$package))
  data.frame(judgement = j, se_used = se, diff = d$simulated - pk,
             sensitive = d$criterion != "info" & !is.na(dev) & 0.1 * abs(d$package) > tol + dev + margin,
             stringsAsFactors = FALSE)
}

val_explained <- function() {
  p <- file.path(val_dir(), "explained.csv")
  if (!file.exists(p)) {
    return(data.frame(id = character(0), reason = character(0)))
  }
  utils::read.csv(p, colClasses = "character", encoding = "UTF-8")
}

# Summary rows of a simulation (common format of python/simulate_ms.py and
# validate_fastsurvival.R)
val_sim_row <- function(source, scenario, hypothesis, look_type, look, time, target, quantity,
                        endpoint, group, x, se = NULL, n = NULL) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  m <- mean(x[ok])
  s <- if (sum(ok) > 1) stats::sd(x[ok]) else 0
  data.frame(source = source, scenario = scenario, hypothesis = hypothesis,
             look_type = look_type, look = as.character(look), time = time,
             target_events = target, quantity = quantity, endpoint = endpoint, group = group,
             mean = m, sd = s, se = if (is.null(se)) s / sqrt(sum(ok)) else se,
             n = if (is.null(n)) sum(ok) else n, stringsAsFactors = FALSE)
}
