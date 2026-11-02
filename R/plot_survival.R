#' Plot of the survival functions of the groups of the design scenario
#'
#' @param x An \code{ms_power} object.
#' @param max_time Largest time since randomization, or \code{NULL} for the
#'   calendar time of the final analysis.
#' @param dots Named list of arguments of \code{plot()}.
#' @return The plotted values (see \code{\link{plot.ms_power}}).
#' @noRd
.plot_survival <- function(x, max_time, dots) {
  if (is.null(max_time)) {
    max_time <- x$study_duration
  }
  if (!is.numeric(max_time) || length(max_time) != 1L || !is.finite(max_time) ||
      !(max_time > 0)) {
    stop("'max_time' must be a positive number.", call. = FALSE)
  }
  scn <- x$scenario
  tt <- seq(0, max_time, length.out = 201L)
  grp <- c("control", "treatment")
  os <- lapply(grp, function(g) ms_surv(scn[[g]], tt, "os")$surv)
  pfs <- lapply(grp, function(g) ms_surv(scn[[g]], tt, "pfs")$surv)
  eps <- if (max(abs(unlist(pfs) - unlist(os))) > 1e-10) c("os", "pfs") else "os"
  d <- do.call(rbind, lapply(eps, function(e) {
    do.call(rbind, lapply(1:2, function(j) {
      lab <- sprintf("%s, %s", c("Control", "Treatment")[j], toupper(e))
      data.frame(time = tt, group = grp[j], endpoint = e,
                 surv = if (e == "os") os[[j]] else pfs[[j]], label = lab,
                 stringsAsFactors = FALSE)
    }))
  }))
  rownames(d) <- NULL
  if (.is_calendar(scn)) {
    attr(d, "note") <- "OS of the patients who enter after the start of switching"
  }
  col <- .plot_colors()
  args <- .plot_args(list(x = c(0, max_time), y = c(0, 1), type = "n",
                          xlab = "Time since randomization", ylab = "Survival probability"),
                     dots)
  do.call(graphics::plot, args)
  key <- unique(d[c("group", "endpoint", "label")])
  for (i in seq_len(nrow(key))) {
    s <- d[d$group == key$group[i] & d$endpoint == key$endpoint[i], ]
    graphics::lines(s$time, s$surv, col = col[[key$group[i]]],
                    lty = if (key$endpoint[i] == "os") 1 else 2)
  }
  graphics::legend("topright", legend = key$label, col = col[key$group],
                   lty = ifelse(key$endpoint == "os", 1, 2), bty = "n")
  if (!is.null(attr(d, "note"))) {
    graphics::mtext(attr(d, "note"), side = 3, line = 0.25, cex = 0.8)
  }
  d
}
