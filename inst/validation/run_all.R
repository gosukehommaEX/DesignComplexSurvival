# Run the validation against independent routes and summarize the results.
#
# Usage (working directory: package root):
#   devtools::load_all()
#   source("inst/validation/run_all.R")
#
# 1. validate_fastsurvival.R simulates V1-V5 with FastSurvival (about 30
#    minutes). Its result results/sim_fastsurvival.csv is reused when
#    helpers.R and validate_fastsurvival.R are unchanged since it was written
#    (checked by MD5); set rerun_fastsurvival <- TRUE before sourcing to force
#    a new simulation.
# 2. compare_simulation.R compares results/sim_python.csv (V6-V8, written by
#    python/simulate_ms.py) and results/sim_fastsurvival.csv with the package.
# 3. compare_pweall.R compares with PWEALL.
# Outputs in inst/validation/results/: compare_simulation.csv,
# compare_pweall.csv, summary.csv, summary.md and perturb.csv. The run stops
# with an error (after writing the outputs) when the Python results are
# missing, when the perturbation check does not give the expected judgements,
# or when an id of explained.csv is not judged EXPLAINED.

local({
  dir <- if (dir.exists(file.path("inst", "validation"))) {
    file.path("inst", "validation")
  } else {
    system.file("validation", package = "DesignComplexSurvival")
  }
  henv <- new.env(parent = globalenv())
  sys.source(file.path(dir, "helpers.R"), envir = henv)
  rd <- henv$val_results_dir()
  problems <- character(0)
  notes <- character(0)
  timing <- numeric(0)

  # 1. FastSurvival
  md5_in <- unname(tools::md5sum(file.path(dir, c("helpers.R", "validate_fastsurvival.R"))))
  side <- file.path(rd, "sim_fastsurvival_inputs_md5.txt")
  fs_file <- file.path(rd, "sim_fastsurvival.csv")
  force <- exists("rerun_fastsurvival", envir = globalenv()) && isTRUE(get("rerun_fastsurvival", envir = globalenv()))
  reuse <- !force && file.exists(fs_file) && file.exists(side) && identical(readLines(side), md5_in)
  if (reuse) {
    message("Reusing results/sim_fastsurvival.csv (inputs unchanged).")
    notes <- c(notes, "FastSurvival results reused (helpers.R and validate_fastsurvival.R unchanged).")
  } else {
    t0 <- proc.time()[["elapsed"]]
    fenv <- new.env(parent = henv)
    sys.source(file.path(dir, "validate_fastsurvival.R"), envir = fenv)
    timing["validate_fastsurvival"] <- proc.time()[["elapsed"]] - t0
    if (isTRUE(fenv$fs_available)) {
      writeLines(md5_in, side)
    } else {
      notes <- c(notes, "FastSurvival is not installed: V1-V5 are not simulated.")
    }
  }
  if (!file.exists(file.path(rd, "sim_python.csv"))) {
    problems <- c(problems, "results/sim_python.csv (python/simulate_ms.py) is missing.")
  }

  # 2. simulations versus the package
  t0 <- proc.time()[["elapsed"]]
  senv <- new.env(parent = henv)
  sys.source(file.path(dir, "compare_simulation.R"), envir = senv)
  cs <- senv$compare_simulation()
  cs$note <- ""
  timing["compare_simulation"] <- proc.time()[["elapsed"]] - t0

  # 3. PWEALL
  t0 <- proc.time()[["elapsed"]]
  penv <- new.env(parent = henv)
  sys.source(file.path(dir, "compare_pweall.R"), envir = penv)
  cp <- penv$compare_pweall()
  timing["compare_pweall"] <- proc.time()[["elapsed"]] - t0
  if (is.null(cp)) {
    notes <- c(notes, "PWEALL is not installed: the comparison with PWEALL is not made.")
  }
  all <- rbind(cs, cp)
  ex <- henv$val_explained()
  j <- henv$val_judge(all, explained = ex$id)
  all <- cbind(all, j)
  utils::write.csv(all[all$source != "pweall", ], file.path(rd, "compare_simulation.csv"), row.names = FALSE)
  if (!is.null(cp)) {
    utils::write.csv(all[all$source == "pweall", ], file.path(rd, "compare_pweall.csv"), row.names = FALSE)
  }
  utils::write.csv(all, file.path(rd, "summary.csv"), row.names = FALSE)
  unused <- setdiff(ex$id, all$id[all$judgement == "EXPLAINED"])
  if (length(unused)) {
    problems <- c(problems, paste("ids of explained.csv not judged EXPLAINED:", paste(unused, collapse = ", ")))
  }

  # perturbation check: package values times 1.1. A sensitive row must not
  # PASS after the change; the other rows are not required to change.
  jd <- all[all$criterion != "info", , drop = FALSE]
  pj <- henv$val_judge(jd, explained = ex$id, perturb = TRUE)
  pt <- data.frame(id = jd$id, judgement = jd$judgement, sensitive = jd$sensitive,
                   judgement_perturbed = pj$judgement,
                   expected = ifelse(jd$sensitive, "FAIL or EXPLAINED", "any"),
                   ok = !jd$sensitive | pj$judgement != "PASS", stringsAsFactors = FALSE)
  utils::write.csv(pt, file.path(rd, "perturb.csv"), row.names = FALSE)
  if (any(!pt$ok)) {
    problems <- c(problems, paste("perturbation check failed for:", paste(pt$id[!pt$ok], collapse = ", ")))
  }

  # summary.md
  fmt <- function(x, d = 6) ifelse(is.na(x), "", formatC(x, digits = d, format = "g"))
  ver <- function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "not installed"
  desc <- if (file.exists("DESCRIPTION")) read.dcf("DESCRIPTION", fields = "Version")[1, 1] else ver("DesignComplexSurvival")
  md5 <- function(f) if (file.exists(f)) unname(tools::md5sum(f)) else "missing"
  lev <- c("PASS", "EXPLAINED", "FAIL", "INFO")
  key <- paste(all$source, all$scenario, sep = " ")
  tab <- table(factor(key, levels = unique(key)), factor(all$judgement, levels = lev))
  l <- c("# Validation against independent routes", "",
         "Generated by `inst/validation/run_all.R`. Do not edit by hand.", "",
         sprintf("- Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
         sprintf("- %s; DesignComplexSurvival %s; FastSurvival %s; PWEALL %s", R.version.string, desc,
                 ver("FastSurvival"), ver("PWEALL")),
         sprintf("- sim_python.csv md5 %s; sim_fastsurvival.csv md5 %s", md5(file.path(rd, "sim_python.csv")),
                 md5(fs_file)),
         sprintf("- Status: %s", if (length(problems)) "CHECKS FAILED" else "all checks passed"),
         if (length(notes)) paste("- Note:", notes), "")
  if (length(problems)) {
    l <- c(l, "## Problems", "", paste("-", problems), "")
  }
  l <- c(l, "## Judgements", "", "| Source and scenario | PASS | EXPLAINED | FAIL | INFO |", "|---|---|---|---|---|",
         sprintf("| %s | %d | %d | %d | %d |", rownames(tab), tab[, "PASS"], tab[, "EXPLAINED"], tab[, "FAIL"],
                 tab[, "INFO"]),
         sprintf("| Total | %d | %d | %d | %d |", sum(tab[, "PASS"]), sum(tab[, "EXPLAINED"]), sum(tab[, "FAIL"]),
                 sum(tab[, "INFO"])), "",
         sprintf("Run times (s): %s.", paste(names(timing), round(timing, 1), sep = " ", collapse = "; ")), "")
  line <- function(d) {
    sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |", d$source, d$scenario, d$hypothesis,
            paste(d$look_type, d$look), fmt(d$time, 4), paste(d$quantity, d$endpoint, d$group),
            d$variant, fmt(d$simulated), fmt(d$package), ifelse(d$se_used > 0, fmt(d$diff / d$se_used, 3), fmt(d$diff, 3)))
  }
  hdr <- c("| Source | Scenario | Hypothesis | Look | Time | Quantity | Variant | Simulated / PWEALL | Package | Diff / SE (or diff) |",
           "|---|---|---|---|---|---|---|---|---|---|")
  f <- all[all$judgement == "FAIL", , drop = FALSE]
  l <- c(l, "## FAIL", "", if (nrow(f)) c(hdr, line(f)) else "None.", "")
  e <- all[all$judgement == "EXPLAINED", , drop = FALSE]
  if (nrow(e)) {
    l <- c(l, "## EXPLAINED", "", hdr, line(e), "")
  }
  pa <- all[all$source == "pweall", , drop = FALSE]
  if (nrow(pa)) {
    l <- c(l, "## PWEALL", "", hdr, line(pa), "")
  }
  inf <- all[all$judgement == "INFO" & all$source != "pweall", , drop = FALSE]
  l <- c(l, "## INFO (asymptotic approximations)", "",
         "Simulated mean (or standard deviation, correlation, rejection rate) and package value; Diff / SE is the difference divided by the Monte Carlo standard error.", "",
         hdr, line(inf), "")
  l <- c(l, "## Perturbation check", "",
         sprintf("Package values of the %d judged rows multiplied by 1.1: %d sensitive rows (the 10%% change exceeds the tolerance plus the observed difference), all expected to FAIL; %d other rows (not required to change).",
                 nrow(pt), sum(pt$sensitive), sum(!pt$sensitive)),
         sprintf("Sensitive rows judged FAIL after the change: %d of %d.", sum(pt$sensitive & pt$judgement_perturbed != "PASS"),
                 sum(pt$sensitive)), "")
  writeLines(l, file.path(rd, "summary.md"), useBytes = TRUE)
  if (length(problems)) {
    stop(paste(c("Validation checks failed:", problems), collapse = "\n  "), call. = FALSE)
  }
  message("Validation finished. Results in ", rd)
  invisible(all)
})
