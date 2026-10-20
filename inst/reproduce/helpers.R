# Helpers for the reproduction scripts in inst/reproduce.
#
# The scripts are run from the package root after devtools::load_all() (or
# with the package installed). Each script reads its published values from
# published_values.csv by id, computes them with the exported functions of
# the package, and returns its results through repro_finish().

if (!exists("ms_size", mode = "function")) {
  library(DesignComplexSurvival)
}

# Folders -------------------------------------------------------------------

repro_dir <- function() {
  d <- file.path("inst", "reproduce")
  if (dir.exists(d)) {
    return(d)
  }
  d <- system.file("reproduce", package = "DesignComplexSurvival")
  if (!nzchar(d)) {
    stop("Cannot find inst/reproduce. Set the working directory to the package root.",
         call. = FALSE)
  }
  d
}

repro_results_dir <- function() {
  d <- file.path(repro_dir(), "results")
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
  d
}

# Published values and registered explanations -------------------------------

repro_published <- function(source = NULL) {
  pv <- utils::read.csv(file.path(repro_dir(), "published_values.csv"),
                        colClasses = "character", encoding = "UTF-8")
  pv$value <- as.numeric(pv$value)
  if (!is.null(source)) {
    pv <- pv[pv$source == source, , drop = FALSE]
  }
  pv
}

repro_explained <- function() {
  utils::read.csv(file.path(repro_dir(), "explained.csv"), colClasses = "character",
                  encoding = "UTF-8")
}

# Common model parts ------------------------------------------------------

# Simultaneous entry, approximated by uniform accrual over a short period
repro_simult <- function(t_end, eps = 1e-6) {
  ms_accrual(time = 0, prop = 1, duration = eps * max(t_end, 1))
}

repro_unif <- function(duration) {
  ms_accrual(time = 0, prop = 1, duration = duration)
}

# Hazards per period from the probability p over length(ratio) periods and
# the ratios of the hazards between the periods
repro_per_period <- function(p, ratio) {
  -log(1 - p) / sum(ratio) * ratio
}

# Two groups with switching in both directions (Lakatos 1988, SIZE), with at
# most two switches. Periods start at 0, 1, 2, ...; hc is the control hazard
# in each period and k the hazard ratio. din: switching from control to the
# experimental treatment; dout: from the experimental treatment to control;
# loss: loss to follow-up. After a switch the hazard is that of the other
# treatment, measured from randomization (clock = "forward").
repro_switch_arms <- function(hc, k, din = 0, dout = 0, loss = 0) {
  np <- length(hc)
  tm <- seq(0, np - 1)
  he <- rep_len(k, np) * hc
  din <- rep_len(din, np)
  dout <- rep_len(dout, np)
  loss <- rep_len(loss, np)
  one <- function(h0, h1, cross1, cross2) {
    if (!any(cross1 > 0)) {
      return(ms_model(time = tm, death = h0, dropout = loss))
    }
    ms_model(time = tm, prog = cross1, death = h0, death_prog = h1,
             switch_rate = cross2, death_switch = h0, dropout = loss,
             time_post = tm, clock = "forward")
  }
  list(control = one(hc, he, din, dout), treatment = one(he, hc, dout, din))
}

# Fixed design: required total sample size (unrounded) for a one-sided test
# at level alpha with the analysis at calendar time t_end
repro_n_fixed <- function(scn, t_end, alpha, power, mean = "direct") {
  ms_size(scn, looks = ms_looks("calendar", time = t_end),
          test = ms_test(alpha = alpha, mean = mean), power = power,
          solve_for = "n")$solution
}

# Fixed design: power of the one-sided test at level alpha at calendar time t
repro_power_fixed <- function(scn, n, t, alpha, variance = "null",
                              endpoint = "os") {
  m <- lr_moments(scn, n = n, calendar_time = t, endpoint = endpoint,
                  variance = variance)
  stats::pnorm((m$z_mean - stats::qnorm(1 - alpha)) / sqrt(m$z_var))
}

# Classical bounds of a design with two looks at information fractions t1
# and 1 and one-sided level alpha: O'Brien-Fleming (c / sqrt(t1), c) or
# Pocock (c, c)
repro_classical2 <- function(type = c("OF", "PO"), alpha, t1) {
  type <- match.arg(type)
  corr <- matrix(c(1, sqrt(t1), sqrt(t1), 1), 2)
  bnd <- function(cc) if (type == "OF") c(cc / sqrt(t1), cc) else c(cc, cc)
  f <- function(cc) sum(gs_prob(c(-Inf, -Inf), bnd(cc), 0, corr)$prob_upper) - alpha
  bnd(stats::uniroot(f, c(1, 5), tol = 1e-12)$root)
}

# Rounding half up (for printed values rounded by hand)
repro_round_half_up <- function(x, digits) {
  floor(x * 10^digits + 0.5 + 1e-9) / 10^digits
}

# SIZE: discrete-time Markov chain of Lakatos (1988) as used by SIZE ------
#
# Transcribed from Planning/phase3_feasibility/size_chain.py. States: on the
# experimental treatment (AE), on control (AC), event, lost. Each period is
# split into sbdv intervals with probabilities 1 - (1 - x)^(1 / sbdv). With
# theta = "prob" the hazard ratio of an interval is the ratio of the interval
# event probabilities (SAS code in the appendix of Lakatos 1988), with
# theta = "hazard" the ratio of the interval hazards (Table 1 of Lakatos
# 1988). accrual gives the proportions of patients entering in each period
# (NULL for simultaneous entry).

size_per_int <- function(hp, sbdv) {
  1 - exp(-rep(hp, each = sbdv) / sbdv)
}

size_chain_run <- function(start, qc, qe, qdin, qdout, ql, qadm, max_cross = Inf) {
  nk <- length(qc)
  ncr <- if (is.finite(max_cross)) as.integer(max_cross) else 50L
  # rows: number of switches so far (0 to ncr); columns: AE, AC
  act <- matrix(0, ncr + 1L, 2L)
  act[1L, if (start == "AE") 1L else 2L] <- 1
  a_risk <- numeric(nk)
  d_ev <- numeric(nk)
  last <- ncr + 1L
  for (i in seq_len(nk)) {
    a_risk[i] <- sum(act)
    d_ev[i] <- sum(act[, 1L] * qe[i] + act[, 2L] * qc[i])
    m_ec <- act[, 1L] * qdout[i]
    m_ce <- act[, 2L] * qdin[i]
    new <- act
    new[, 1L] <- new[, 1L] - (act[, 1L] * (qe[i] + ql[i] + qadm[i]) + m_ec)
    new[, 2L] <- new[, 2L] - (act[, 2L] * (qc[i] + ql[i] + qadm[i]) + m_ce)
    new[-1L, 2L] <- new[-1L, 2L] + m_ec[-last]
    new[-1L, 1L] <- new[-1L, 1L] + m_ce[-last]
    # no further switch after max_cross switches
    new[last, 1L] <- new[last, 1L] + m_ec[last]
    new[last, 2L] <- new[last, 2L] + m_ce[last]
    act <- new
  }
  list(A = a_risk, D = d_ev)
}

size_calc <- function(hc, k, sbdv, din = 0, dout = 0, loss = 0, ratio_ce = 1,
                      accrual = NULL, max_cross = Inf, n = NULL,
                      za = stats::qnorm(0.975), zb = stats::qnorm(0.8),
                      theta = c("prob", "hazard")) {
  theta <- match.arg(theta)
  npd <- length(hc)
  nk <- npd * sbdv
  bc <- function(v) {
    v <- as.numeric(v)
    if (!(length(v) %in% c(1L, npd))) {
      stop("Inputs must have length 1 or the number of periods.", call. = FALSE)
    }
    rep_len(v, npd)
  }
  qc <- size_per_int(bc(hc), sbdv)
  qe <- size_per_int(bc(k) * bc(hc), sbdv)
  qdin <- size_per_int(bc(din), sbdv)
  qdout <- size_per_int(bc(dout), sbdv)
  ql <- size_per_int(bc(loss), sbdv)
  qadm <- numeric(nk)
  if (!is.null(accrual)) {
    # Shih (1995): censoring with probability p_k / sum_{i <= k} p_i in
    # follow-up interval K - k + 1
    p <- rep(as.numeric(accrual) / sbdv, each = sbdv)
    cs <- cumsum(p)
    for (m in seq_len(nk)) {
      kk <- nk - m + 1L
      qadm[m] <- if (cs[kk] > 0) p[kk] / cs[kk] else 0
    }
  }
  rc <- size_chain_run("AC", qc, qe, qdin, qdout, ql, qadm, max_cross)
  re <- size_chain_run("AE", qc, qe, qdin, qdout, ql, qadm, max_cross)
  wc <- ratio_ce / (1 + ratio_ce)
  we <- 1 - wc
  phi <- (wc * rc$A) / (we * re$A)
  th <- if (theta == "prob") {
    (rc$D / rc$A) / (re$D / re$A)
  } else {
    log(1 - rc$D / rc$A) / log(1 - re$D / re$A)
  }
  dd <- wc * rc$D + we * re$D
  rho <- dd / sum(dd)
  gam <- phi * th / (1 + phi * th) - phi / (1 + phi)
  eta <- phi / (1 + phi)^2
  pc <- sum(rc$D)
  pe <- sum(re$D)
  pbar <- wc * pc + we * pe
  if (is.null(n)) {
    d <- (za + zb)^2 * sum(rho * eta) / sum(rho * gam)^2
    return(list(n = d / pbar, d = d, PC = pc, PE = pe))
  }
  d <- n * pbar
  list(power = stats::pnorm(sqrt(d) * abs(sum(rho * gam)) / sqrt(sum(rho * eta)) - za),
       d = d, PC = pc, PE = pe)
}

# Judgement -----------------------------------------------------------------
#
# Criteria (published value P, computed value x):
#   ceil                    P - 1.05 < x <= P + 0.05
#   even_ceil               P - 2.05 < x <= P + 0.05
#   ceil_to:step            P - 1.05 step < x <= P + 0.05 step (rounded up to a multiple of step)
#   abs:tol                 |x - P| <= tol
#   power_at_time:target:tol  x is the power at the published time P, |x - target| <= tol
#   round_n_from_d          x rounded to the nearest integer equals P
#   info                    not judged

repro_check <- function(criterion, published, computed) {
  if (is.na(computed) || is.na(published)) {
    return(list(pass = FALSE, diff = NA_real_))
  }
  eps <- 1e-9
  parts <- strsplit(criterion, ":", fixed = TRUE)[[1L]]
  type <- parts[1L]
  if (type == "ceil") {
    return(list(pass = computed > published - 1.05 && computed <= published + 0.05 + eps,
                diff = computed - published))
  }
  if (type == "even_ceil") {
    return(list(pass = computed > published - 2.05 && computed <= published + 0.05 + eps,
                diff = computed - published))
  }
  if (type == "ceil_to") {
    step <- as.numeric(parts[2L])
    return(list(pass = computed > published - 1.05 * step &&
                  computed <= published + 0.05 * step + eps,
                diff = computed - published))
  }
  if (type == "abs") {
    tol <- as.numeric(parts[2L])
    return(list(pass = abs(computed - published) <= tol + eps,
                diff = computed - published))
  }
  if (type == "power_at_time") {
    target <- as.numeric(parts[2L])
    tol <- as.numeric(parts[3L])
    return(list(pass = abs(computed - target) <= tol + eps, diff = computed - target))
  }
  if (type == "round_n_from_d") {
    return(list(pass = abs(floor(computed + 0.5) - published) < eps,
                diff = computed - published))
  }
  stop("Unknown criterion: ", criterion, call. = FALSE)
}

# One result row. computed is compared with the published value; for the
# criterion power_at_time it is the power at the published time.
# evidence: value used when the row fails and has an explained_id (see
# explained.csv); evidence_criterion defaults to the criterion of the row;
# evidence_ref_id gives another row whose published value is the reference
# (defaults to the row itself); evidence_extra is an additional logical
# condition (for example, agreement of the package with a fine
# discretization). The *_perturbed values are used by the perturbation check
# (published value times 1.1): for power_at_time they are the powers at 1.1
# times the published time.
repro_row <- function(id, computed, note = "", evidence = NA_real_,
                      evidence_criterion = NA_character_,
                      evidence_ref_id = NA_character_, evidence_extra = NA,
                      computed_perturbed = computed,
                      evidence_perturbed = evidence) {
  data.frame(id = id, computed = as.numeric(computed), note = note,
             evidence = as.numeric(evidence),
             evidence_criterion = evidence_criterion,
             evidence_ref_id = evidence_ref_id,
             evidence_extra = as.logical(evidence_extra),
             computed_perturbed = as.numeric(computed_perturbed),
             evidence_perturbed = as.numeric(evidence_perturbed),
             stringsAsFactors = FALSE)
}

repro_diag <- function(id, value, reference, tol, note = "") {
  rel <- abs(value / reference - 1)
  data.frame(id = id, value = value, reference = reference, rel_diff = rel,
             tol = tol, status = if (rel <= tol) "OK" else "NG", note = note,
             stringsAsFactors = FALSE)
}

repro_judge_one <- function(criterion, explained_id, published, computed,
                            evidence, evidence_criterion, evidence_ref,
                            evidence_extra) {
  if (criterion == "info") {
    return(list(judgement = "INFO", diff = computed - published, evidence_pass = NA))
  }
  main <- repro_check(criterion, published, computed)
  if (main$pass) {
    return(list(judgement = "PASS", diff = main$diff, evidence_pass = NA))
  }
  if (!nzchar(explained_id)) {
    return(list(judgement = "FAIL", diff = main$diff, evidence_pass = NA))
  }
  crit <- if (is.na(evidence_criterion)) criterion else evidence_criterion
  ev_pass <- !is.na(evidence) && repro_check(crit, evidence_ref, evidence)$pass &&
    (is.na(evidence_extra) || evidence_extra)
  list(judgement = if (ev_pass) "EXPLAINED" else "FAIL", diff = main$diff,
       evidence_pass = ev_pass)
}

# Judge all rows. With perturb = TRUE, each published value is multiplied by
# 1.1 (one row at a time; a reference value taken from another row is not
# changed) and the perturbed computed values are used.
repro_judge <- function(res, pv, perturb = FALSE) {
  idx <- match(res$id, pv$id)
  f <- if (perturb) 1.1 else 1
  pub <- pv$value[idx] * f
  ref <- ifelse(is.na(res$evidence_ref_id), pub,
                pv$value[match(res$evidence_ref_id, pv$id)])
  comp <- if (perturb) res$computed_perturbed else res$computed
  ev <- if (perturb) res$evidence_perturbed else res$evidence
  out <- lapply(seq_len(nrow(res)), function(i) {
    repro_judge_one(pv$criterion[idx[i]], pv$explained_id[idx[i]], pub[i], comp[i],
                    ev[i], res$evidence_criterion[i], ref[i], res$evidence_extra[i])
  })
  data.frame(id = res$id, published = pub, computed = comp,
             diff = vapply(out, function(z) z$diff, numeric(1)),
             evidence_pass = vapply(out, function(z) z$evidence_pass, logical(1)),
             judgement = vapply(out, function(z) z$judgement, character(1)),
             stringsAsFactors = FALSE)
}

# Finish a script: check that every published value of the source has
# exactly one row, judge the rows and write results/<source>.csv.
repro_finish <- function(source, rows, diagnostics = NULL, helpers = NULL) {
  res <- do.call(rbind, rows)
  pv <- repro_published(source)
  dup <- res$id[duplicated(res$id)]
  missing <- setdiff(pv$id, res$id)
  extra <- setdiff(res$id, pv$id)
  if (length(dup) || length(missing) || length(extra)) {
    stop(sprintf("%s: duplicated rows: %s; missing: %s; unknown ids: %s", source,
                 paste(dup, collapse = ", "), paste(missing, collapse = ", "),
                 paste(extra, collapse = ", ")), call. = FALSE)
  }
  res <- res[match(pv$id, res$id), , drop = FALSE]
  bad_ref <- setdiff(stats::na.omit(res$evidence_ref_id), pv$id)
  if (length(bad_ref)) {
    stop(source, ": unknown evidence_ref_id ", paste(bad_ref, collapse = ", "),
         call. = FALSE)
  }
  j <- repro_judge(res, pv)
  out <- data.frame(id = res$id, source = source, quantity = pv$quantity,
                    criterion = pv$criterion, explained_id = pv$explained_id,
                    published = pv$value, computed = res$computed, diff = j$diff,
                    evidence = res$evidence,
                    evidence_criterion = res$evidence_criterion,
                    evidence_ref_id = res$evidence_ref_id,
                    evidence_extra = res$evidence_extra,
                    evidence_pass = j$evidence_pass, judgement = j$judgement,
                    computed_perturbed = res$computed_perturbed,
                    evidence_perturbed = res$evidence_perturbed,
                    note = trimws(paste(pv$note, res$note, sep = " | "), whitespace = "[ |]"),
                    stringsAsFactors = FALSE)
  rd <- repro_results_dir()
  utils::write.csv(out, file.path(rd, paste0(source, ".csv")), row.names = FALSE,
                   fileEncoding = "UTF-8")
  if (!is.null(helpers)) {
    utils::write.csv(helpers, file.path(rd, paste0("helpers_", source, ".csv")),
                     row.names = FALSE, fileEncoding = "UTF-8")
  }
  tab <- table(factor(out$judgement, levels = c("PASS", "EXPLAINED", "FAIL", "INFO")))
  message(source, ": ", paste(names(tab), tab, sep = " ", collapse = ", "))
  list(results = out, diagnostics = diagnostics)
}
