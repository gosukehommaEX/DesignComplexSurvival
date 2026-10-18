#' Multistate model for one treatment group
#'
#' Defines the four-state model of one group: progression-free and alive
#' (state 1), progressed without switch (state 2), switched (state 3) and
#' dead (state 4). Progression-free survival (PFS) is the time to leave
#' state 1 and overall survival (OS) is the time to enter state 4.
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
#' Dropout is piecewise exponential with intervals starting at \code{time},
#' independent of the event process, and censors both endpoints.
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
#'
#' @return An object of class \code{ms_model}: a list with the intensity
#'   vectors expanded to the lengths of \code{time} and \code{time_post},
#'   and \code{clock}.
#'
#' @references
#' Fleischer F, Gaschler-Markefski B, Bluhmki E (2009). A statistical model
#' for the dependence between progression-free survival and overall
#' survival. Statistics in Medicine 28:2669-2686.
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
#' @export
ms_model <- function(time = 0, prog = 0, death, death_prog = NULL,
                     switch_prob = 0, switch_rate = 0, death_switch = NULL,
                     dropout = 0, time_post = 0, clock = c("reset", "forward")) {
  clock <- match.arg(clock)
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
  if (is.null(death_prog)) {
    if (any(prog > 0)) {
      stop("'death_prog' must be specified when 'prog' is positive.", call. = FALSE)
    }
    death_prog <- 0
  }
  death_prog <- .check_rates(death_prog, np, "death_prog")
  switch_rate <- .check_rates(switch_rate, np, "switch_rate")
  death_switch <- if (is.null(death_switch)) {
    death_prog
  } else {
    .check_rates(death_switch, np, "death_switch")
  }
  structure(list(time = as.numeric(time), prog = prog, death = death,
                 switch_prob = switch_prob, dropout = dropout,
                 time_post = as.numeric(time_post), death_prog = death_prog,
                 switch_rate = switch_rate, death_switch = death_switch,
                 clock = clock),
            class = "ms_model")
}
