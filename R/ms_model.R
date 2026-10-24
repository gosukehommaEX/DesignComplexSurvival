#' Multistate model for one treatment group
#'
#' Defines the multistate model of one group: progression-free and alive
#' (state 1), progressed without switch (state 2), switched (state 3) and
#' dead (state 4), and, when \code{switch_pre} is positive, progression-free
#' after a switch before progression (state 5). Progression-free survival
#' (PFS) is the time to leave the progression-free states and overall
#' survival (OS) is the time to enter state 4.
#'
#' @details
#' Before progression, the intensities are piecewise constant in the time
#' since randomization, with intervals starting at \code{time}: progression
#' (\code{prog}), death before progression (\code{death}), and the
#' probability \code{switch_prob} of switching immediately at progression
#' (evaluated at the time of progression). After progression without
#' switch, the patient dies with intensity \code{death_prog} or switches
#' with intensity \code{switch_rate}; after a switch, the patient dies with
#' intensity \code{death_switch}. The post-progression intensities are
#' piecewise constant with intervals starting at \code{time_post}.
#'
#' With \code{clock = "reset"} (default) the post-progression intensities
#' are functions of the time since progression (a semi-Markov model). With
#' \code{clock = "forward"} they are functions of the time since
#' randomization (a Markov model). The two clocks give the same model when
#' the post-progression intensities are constant.
#'
#' The death hazard after a switch at the post-progression clock value
#' \eqn{\sigma} is, at the clock value \eqn{\tau}, \code{death_switch} at
#' \eqn{\tau} with \code{switch_clock = "same"} (default), and
#' \code{death_switch} at \eqn{\tau - \sigma} (the time since the switch)
#' with \code{switch_clock = "switch"}. With \code{clock = "reset"} the two
#' differ only for switches after progression (\code{switch_rate}); with
#' \code{clock = "forward"} they differ also for switches at progression.
#' With \code{switch_af} given, the effect of a switch is an accelerated
#' failure time effect: the death hazard after the switch is
#' \eqn{\phi c(\sigma + \phi (\tau - \sigma))}, where \eqn{\phi} is
#' \code{switch_af} and \eqn{c} is \code{death_prog}, so that the remaining
#' time is that without switch divided by \eqn{\phi} (the rank preserving
#' structural failure time model of Robins and Tsiatis 1991; \eqn{\phi < 1}
#' prolongs survival).
#'
#' A switch before progression moves a patient from state 1 to state 5 with
#' hazard \code{switch_pre}. In state 5 the patient progresses with hazard
#' \code{prog_switch} and dies with hazard \code{death_switch_pre}; at
#' progression the patient enters state 3 and then dies with the hazard
#' \code{death_switch} on the post-progression clock, whatever
#' \code{switch_clock} (with \code{switch_af}, the hazard
#' \code{switch_af * death_prog}). Whether a switch before progression
#' censors PFS is set in \code{\link{ms_scenario}} (\code{pfs_switch}).
#'
#' Dropout is piecewise exponential with intervals starting at \code{time},
#' independent of the event process, and censors both endpoints. The hazard
#' \code{dropout_pfs} adds a censoring that acts on PFS only (for example
#' the end of tumour assessments), independent of everything else.
#'
#' A model without progression (\code{prog = 0}) is a piecewise
#' exponential model for a single endpoint with hazard \code{death}. With
#' \code{death_prog = death}, \code{switch_prob = 0} and
#' \code{switch_rate = 0}, the model is that of Fleischer et al. (2009); see
#' \code{\link{ms_fleischer}}.
#'
#' @param time Interval start times for the pre-progression intensities and
#'   the dropout hazard. Must start at 0 and be strictly increasing.
#' @param prog Progression hazard in each interval (length 1 or
#'   \code{length(time)}).
#' @param death Hazard of death before progression in each interval.
#' @param death_prog Hazard of death after progression without switch, in
#'   each interval of \code{time_post}. Required when \code{prog} is
#'   positive in some interval.
#' @param switch_prob Probability of switching immediately at progression,
#'   in each interval of \code{time}.
#' @param switch_rate Hazard of switching after progression, in each
#'   interval of \code{time_post}.
#' @param death_switch Hazard of death after a switch, in each interval of
#'   \code{time_post}. Defaults to \code{death_prog}, that is, switching has
#'   no effect on survival.
#' @param dropout Dropout hazard in each interval of \code{time}.
#' @param time_post Interval start times for the post-progression
#'   intensities. Must start at 0 and be strictly increasing.
#' @param clock Clock of the post-progression intensities: \code{"reset"}
#'   (time since progression) or \code{"forward"} (time since
#'   randomization).
#' @param switch_clock Clock of \code{death_switch}: \code{"same"} (the
#'   clock of the other post-progression intensities) or \code{"switch"}
#'   (time since the switch).
#' @param switch_af Acceleration factor of an accelerated failure time effect
#'   of the switch (positive), or \code{NULL} (default) for the hazard
#'   \code{death_switch}. When given, \code{death_switch} must not be given
#'   and \code{switch_clock} must be \code{"same"}.
#' @param switch_pre Hazard of switching before progression, in each interval
#'   of \code{time}.
#' @param prog_switch Hazard of progression after a switch before
#'   progression, in each interval of \code{time}. Defaults to \code{prog}.
#' @param death_switch_pre Hazard of death before progression after a switch,
#'   in each interval of \code{time}. Defaults to \code{death}.
#' @param dropout_pfs Hazard of a censoring of PFS only, in each interval of
#'   \code{time}.
#'
#' @return An object of class \code{ms_model}: a list with the intensity
#'   vectors expanded to the lengths of \code{time} and \code{time_post},
#'   \code{clock}, \code{switch_clock}, \code{switch_af} and
#'   \code{pfs_switch} (\code{"event"}; set by \code{\link{ms_scenario}}).
#'   With \code{switch_af}, the element \code{death_switch} holds
#'   \code{switch_af * death_prog}, the death hazard after a switch when
#'   \code{death_prog} is constant.
#'
#' @references
#' Fleischer F, Gaschler-Markefski B, Bluhmki E (2009). A statistical model
#' for the dependence between progression-free survival and overall
#' survival. Statistics in Medicine 28:2669-2686.
#'
#' Robins JM, Tsiatis AA (1991). Correcting for non-compliance in randomized
#' trials using rank preserving structural failure time models.
#' Communications in Statistics - Theory and Methods 20:2609-2631.
#'
#' @seealso \code{\link{ms_fleischer}}, \code{\link{ms_scenario}},
#'   \code{\link{ms_surv}}, \code{\link{ms_prob}}
#'
#' @examples
#' # Single endpoint: exponential survival with a median of 12 months
#' ms_model(death = log(2) / 12)
#'
#' # Progression, with 30 percent of the patients switching at progression
#' ms_model(time = c(0, 6), prog = c(0.10, 0.08), death = 0.01,
#'          death_prog = 0.05, switch_prob = 0.3, death_switch = 0.03)
#'
#' # Switching after progression with an accelerated failure time effect
#' ms_model(prog = 0.1, death = 0.01, death_prog = 0.05, switch_rate = 0.04,
#'          switch_af = 0.7)
#'
#' # Switching before progression, and censoring of PFS only
#' ms_model(prog = 0.1, death = 0.01, death_prog = 0.05, death_switch = 0.03,
#'          switch_pre = 0.02, dropout_pfs = 0.01)
#'
#' @export
ms_model <- function(time = 0, prog = 0, death, death_prog = NULL,
                     switch_prob = 0, switch_rate = 0, death_switch = NULL,
                     dropout = 0, time_post = 0, clock = c("reset", "forward"),
                     switch_clock = c("same", "switch"), switch_af = NULL,
                     switch_pre = 0, prog_switch = NULL, death_switch_pre = NULL,
                     dropout_pfs = 0) {
  clock <- match.arg(clock)
  switch_clock <- match.arg(switch_clock)
  if (missing(death)) {
    stop("'death' must be specified.", call. = FALSE)
  }
  .check_breaks(time, "time")
  .check_breaks(time_post, "time_post")
  nt <- length(time)
  np <- length(time_post)
  prog <- .check_rates(prog, nt, "prog")
  death <- .check_rates(death, nt, "death")
  switch_prob <- .check_rates(switch_prob, nt, "switch_prob", upper = 1)
  dropout <- .check_rates(dropout, nt, "dropout")
  switch_pre <- .check_rates(switch_pre, nt, "switch_pre")
  prog_switch <- if (is.null(prog_switch)) prog else .check_rates(prog_switch, nt, "prog_switch")
  death_switch_pre <- if (is.null(death_switch_pre)) {
    death
  } else {
    .check_rates(death_switch_pre, nt, "death_switch_pre")
  }
  dropout_pfs <- .check_rates(dropout_pfs, nt, "dropout_pfs")
  if (is.null(death_prog)) {
    if (any(prog > 0)) {
      stop("'death_prog' must be specified when 'prog' is positive.", call. = FALSE)
    }
    if (any(switch_pre > 0) && any(prog_switch > 0) && is.null(death_switch)) {
      stop(paste0("'death_prog' or 'death_switch' must be specified when 'switch_pre' ",
                  "and 'prog_switch' are positive."), call. = FALSE)
    }
    death_prog <- 0
  }
  death_prog <- .check_rates(death_prog, np, "death_prog")
  switch_rate <- .check_rates(switch_rate, np, "switch_rate")
  if (!is.null(switch_af)) {
    if (!is.numeric(switch_af) || length(switch_af) != 1L || is.na(switch_af) ||
        !is.finite(switch_af) || switch_af <= 0) {
      stop("'switch_af' must be NULL or one positive number.", call. = FALSE)
    }
    if (!is.null(death_switch)) {
      stop("'death_switch' must not be given with 'switch_af'.", call. = FALSE)
    }
    if (switch_clock != "same") {
      stop("'switch_clock' must be \"same\" with 'switch_af'.", call. = FALSE)
    }
    switch_af <- as.numeric(switch_af)
    death_switch <- switch_af * death_prog
  } else if (is.null(death_switch)) {
    death_switch <- death_prog
  } else {
    death_switch <- .check_rates(death_switch, np, "death_switch")
  }
  structure(list(time = as.numeric(time), prog = prog, death = death,
                 switch_prob = switch_prob, dropout = dropout,
                 time_post = as.numeric(time_post), death_prog = death_prog,
                 switch_rate = switch_rate, death_switch = death_switch,
                 clock = clock, switch_clock = switch_clock, switch_af = switch_af,
                 switch_pre = switch_pre, prog_switch = prog_switch,
                 death_switch_pre = death_switch_pre, dropout_pfs = dropout_pfs,
                 pfs_switch = "event"),
            class = "ms_model")
}
