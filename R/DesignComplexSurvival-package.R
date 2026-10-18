#' DesignComplexSurvival: Simulation-Free Design of Complex Time-to-Event Trials
#'
#' Computes the operating characteristics of two-group time-to-event trials
#' analysed with the log-rank test, without patient-level simulation. Each
#' group is described by a four-state model (progression-free, progressed,
#' switched, dead) with piecewise constant intensities, so that
#' progression-free survival, overall survival and treatment switching after
#' progression are handled within one model. Accrual is piecewise uniform and
#' dropout is piecewise exponential.
#'
#' Expected numbers of events, the mean, variance and covariance of the
#' log-rank statistic over calendar time, and the probabilities of group
#' sequential designs with efficacy and futility bounds are obtained from
#' closed-form expressions, Gauss-Legendre quadrature and multivariate normal
#' probabilities.
#'
#' @section Workflow:
#' \enumerate{
#'   \item Define each group with \code{\link{ms_model}} (or
#'     \code{\link{ms_fleischer}}) and the accrual with
#'     \code{\link{ms_accrual}}, and combine them with
#'     \code{\link{ms_scenario}}.
#'   \item Define the analysis times with \code{\link{ms_looks}} and the test
#'     with \code{\link{ms_test}}.
#'   \item Compute the power of a given design with \code{\link{ms_power}} or
#'     solve for the number of events, the sample size, the accrual duration
#'     or the follow-up with \code{\link{ms_size}}.
#' }
#' The building blocks \code{\link{ms_surv}}, \code{\link{ms_prob}},
#' \code{\link{ms_events}}, \code{\link{ms_event_time}},
#' \code{\link{lr_moments}}, \code{\link{lr_cov}}, \code{\link{gs_spend}},
#' \code{\link{gs_bounds}} and \code{\link{gs_prob}} are exported so that each
#' step can be checked separately.
#'
#' @keywords internal
#' @importFrom stats pnorm qnorm uniroot
#' @importFrom mvtnorm pmvnorm
"_PACKAGE"

# Cache for quadrature rules
.dcs_env <- new.env(parent = emptyenv())
