# Run all reproduction scripts and summarize the results.
#
# Usage (working directory: package root):
#   devtools::load_all()
#   source("inst/reproduce/run_all.R")
#
# Writes to inst/reproduce/results/: <source>.csv for each script,
# summary.csv, summary.md, perturb.csv, diagnostics.csv and helpers_*.csv.
# The run stops with an error (after writing the outputs) when
#   - a published value has no result row or more than one,
#   - an explained_id of published_values.csv is not in explained.csv,
#   - an id of explained.csv never gives the judgement EXPLAINED,
#   - the perturbation check does not give the expected judgements, or
#   - a diagnostic is outside its tolerance.

local({
  dir <- if (dir.exists(file.path("inst", "reproduce"))) {
    file.path("inst", "reproduce")
  } else {
    system.file("reproduce", package = "DesignComplexSurvival")
  }
  henv <- new.env(parent = globalenv())
  sys.source(file.path(dir, "helpers.R"), envir = henv)

  reproduce_all <- function(perturb = TRUE) {
    sources <- c("Barthel2006", "Glimm2010", "Lakatos1988", "Lakatos2002", "Lu2021",
                 "Luo2019", "Porcher2002", "Shih1995", "Tamhane2010", "Tang2022",
                 "Vergote2024", "Xia2016")
    res <- list()
    diag <- list()
    timing <- numeric(0)
    for (s in sources) {
      env <- new.env(parent = henv)
      t0 <- proc.time()[["elapsed"]]
      sys.source(file.path(dir, paste0("reproduce_", s, ".R")), envir = env)
      timing[s] <- proc.time()[["elapsed"]] - t0
      res[[s]] <- env$result$results
      if (!is.null(env$result$diagnostics)) {
        diag[[s]] <- env$result$diagnostics
      }
    }
    out <- do.call(rbind, res)
    rownames(out) <- NULL
    pv <- henv$repro_published()
    ex <- henv$repro_explained()
    rd <- henv$repro_results_dir()
    problems <- character(0)

    # one result row per published value
    dup <- out$id[duplicated(out$id)]
    miss <- setdiff(pv$id, out$id)
    extra <- setdiff(out$id, pv$id)
    if (length(dup) || length(miss) || length(extra)) {
      problems <- c(problems, sprintf("ids: duplicated %s; missing %s; unknown %s",
                                      paste(dup, collapse = ", "), paste(miss, collapse = ", "),
                                      paste(extra, collapse = ", ")))
    }
    # explained ids
    used <- unique(pv$explained_id[nzchar(pv$explained_id)])
    unreg <- setdiff(used, ex$id)
    if (length(unreg)) {
      problems <- c(problems, paste("explained_id not in explained.csv:",
                                    paste(unreg, collapse = ", ")))
    }
    shown <- unique(out$explained_id[out$judgement == "EXPLAINED"])
    unshown <- setdiff(ex$id, shown)
    if (length(unshown)) {
      problems <- c(problems, paste("ids of explained.csv without EXPLAINED rows:",
                                    paste(unshown, collapse = ", ")))
    }

    # perturbation check (published values times 1.1, one row at a time)
    pt <- NULL
    if (perturb) {
      judged <- out[out$criterion != "info", , drop = FALSE]
      pj <- henv$repro_judge(judged, pv, perturb = TRUE)
      expect <- rep("FAIL", nrow(judged))
      has_ex <- nzchar(judged$explained_id) & judged$judgement == "EXPLAINED"
      ie <- match(judged$explained_id, ex$id)
      expect[has_ex & !is.na(ie)] <- ex$perturb_expect[ie[has_ex & !is.na(ie)]]
      pt <- data.frame(id = judged$id, source = judged$source, judgement = judged$judgement,
                       published_perturbed = pj$published,
                       computed_perturbed = pj$computed,
                       judgement_perturbed = pj$judgement, expected = expect,
                       ok = pj$judgement == expect, stringsAsFactors = FALSE)
      utils::write.csv(pt, file.path(rd, "perturb.csv"), row.names = FALSE)
      if (any(!pt$ok)) {
        problems <- c(problems, paste("perturbation check failed for:",
                                      paste(pt$id[!pt$ok], collapse = ", ")))
      }
    }
    dg <- if (length(diag)) do.call(rbind, diag) else NULL
    if (!is.null(dg)) {
      rownames(dg) <- NULL
      utils::write.csv(dg, file.path(rd, "diagnostics.csv"), row.names = FALSE)
      if (any(dg$status != "OK")) {
        problems <- c(problems, paste("diagnostics outside tolerance:",
                                      paste(dg$id[dg$status != "OK"], collapse = ", ")))
      }
    }
    utils::write.csv(out, file.path(rd, "summary.csv"), row.names = FALSE,
                     fileEncoding = "UTF-8")
    write_summary_md(out, pv, ex, pt, dg, timing, problems, file.path(rd, "summary.md"))
    if (length(problems)) {
      stop(paste(c("Reproduction checks failed:", problems), collapse = "\n  "), call. = FALSE)
    }
    message("All checks passed. Results in ", rd)
    invisible(out)
  }

  write_summary_md <- function(out, pv, ex, pt, dg, timing, problems, path) {
    fmt <- function(x) ifelse(is.na(x), "", formatC(x, digits = 6, format = "g"))
    lev <- c("PASS", "EXPLAINED", "FAIL", "INFO")
    tab <- table(factor(out$source, levels = unique(out$source)),
                 factor(out$judgement, levels = lev))
    desc <- file.path(if (file.exists("DESCRIPTION")) "." else dirname(dirname(dir)), "DESCRIPTION")
    ver <- if (file.exists(desc)) read.dcf(desc, fields = "Version")[1, 1] else
      as.character(utils::packageVersion("DesignComplexSurvival"))
    md5 <- tools::md5sum(file.path(dir, c("published_values.csv", "explained.csv")))
    l <- c("# Reproduction of published results",
           "",
           "Generated by `inst/reproduce/run_all.R`. Do not edit by hand.",
           "",
           sprintf("- Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
           sprintf("- %s; DesignComplexSurvival %s", R.version.string, ver),
           sprintf("- published_values.csv md5 %s; explained.csv md5 %s", md5[[1]], md5[[2]]),
           sprintf("- Status: %s", if (length(problems)) "CHECKS FAILED" else "all checks passed"),
           "")
    if (length(problems)) {
      l <- c(l, "## Problems", "", paste("-", problems), "")
    }
    l <- c(l, "## Judgements by source", "",
           "| Source | PASS | EXPLAINED | FAIL | INFO | Total | Time (s) |",
           "|---|---|---|---|---|---|---|")
    for (s in rownames(tab)) {
      l <- c(l, sprintf("| %s | %d | %d | %d | %d | %d | %.1f |", s, tab[s, "PASS"],
                        tab[s, "EXPLAINED"], tab[s, "FAIL"], tab[s, "INFO"], sum(tab[s, ]),
                        timing[[s]]))
    }
    tot <- colSums(tab)
    l <- c(l, sprintf("| Total | %d | %d | %d | %d | %d | %.1f |", tot[["PASS"]],
                      tot[["EXPLAINED"]], tot[["FAIL"]], tot[["INFO"]], sum(tot), sum(timing)),
           "")
    row_line <- function(d) {
      sprintf("| %s | %s | %s | %s | %s | %s |", d$id, d$criterion, fmt(d$published),
              fmt(d$computed), fmt(d$evidence), gsub("|", "/", d$note, fixed = TRUE))
    }
    f <- out[out$judgement == "FAIL", , drop = FALSE]
    l <- c(l, "## FAIL", "")
    if (nrow(f)) {
      l <- c(l, "| id | criterion | published | computed | evidence | note |",
             "|---|---|---|---|---|---|", row_line(f), "")
    } else {
      l <- c(l, "None.", "")
    }
    e <- out[out$judgement == "EXPLAINED", , drop = FALSE]
    l <- c(l, "## EXPLAINED", "")
    for (k in unique(e$explained_id)) {
      ek <- e[e$explained_id == k, , drop = FALSE]
      exk <- ex[ex$id == k, , drop = FALSE]
      l <- c(l, sprintf("### %s (%d rows)", k, nrow(ek)), "",
             paste("Reason:", exk$reason), "", paste("Evidence:", exk$evidence), "",
             "| id | criterion | published | computed | evidence | note |",
             "|---|---|---|---|---|---|", row_line(ek), "")
    }
    inf <- out[out$judgement == "INFO", , drop = FALSE]
    l <- c(l, "## INFO", "",
           "| id | published | computed | computed / published - 1 | note |",
           "|---|---|---|---|---|",
           sprintf("| %s | %s | %s | %s | %s |", inf$id, fmt(inf$published), fmt(inf$computed),
                   fmt(inf$computed / inf$published - 1),
                   gsub("|", "/", inf$note, fixed = TRUE)), "")
    if (!is.null(pt)) {
      n_ok <- sum(pt$ok)
      l <- c(l, "## Perturbation check", "",
             sprintf("Each judged published value (%d rows) was multiplied by 1.1 and judged again.", nrow(pt)),
             sprintf("Rows with the expected judgement: %d of %d.", n_ok, nrow(pt)),
             sprintf("Expected FAIL: %d; expected EXPLAINED (evidence taken from another row): %d (%s).",
                     sum(pt$expected == "FAIL"), sum(pt$expected == "EXPLAINED"),
                     paste(pt$id[pt$expected == "EXPLAINED"], collapse = ", ")), "")
      if (any(!pt$ok)) {
        l <- c(l, "Rows without the expected judgement:", "", paste("-", pt$id[!pt$ok]), "")
      }
    }
    if (!is.null(dg)) {
      l <- c(l, "## Diagnostics", "",
             "| id | value | reference | relative difference | tolerance | status | note |",
             "|---|---|---|---|---|---|---|",
             sprintf("| %s | %s | %s | %s | %s | %s | %s |", dg$id,
                     formatC(dg$value, digits = 10, format = "g"),
                     formatC(dg$reference, digits = 10, format = "g"),
                     formatC(dg$rel_diff, digits = 3, format = "e"),
                     formatC(dg$tol, digits = 1, format = "e"), dg$status, dg$note), "")
    }
    writeLines(l, path, useBytes = TRUE)
  }

  reproduce_all()
})
