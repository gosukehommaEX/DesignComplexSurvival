#' DesignComplexSurvival: Simulation-Free Design of Complex Time-to-Event Trials
#'
#' Computes the operating characteristics of two-group time-to-event trials
#' analyzed with the log-rank test, without patient-level simulation. Each
#' group is described by a multistate model with piecewise constant
#' intensities: progression-free, progressed, switched and dead, and, when
#' patients may switch before progression, progression-free after a switch.
#' Progression-free survival, overall survival and treatment switching
#' before, at or after progression are handled within one model. The effect
#' of a switch is a change of the death hazard or an accelerated failure time
#' effect, and switching may start at a calendar time. Accrual is piecewise
#' uniform, dropout is piecewise exponential, and progression-free survival
#' may have a censoring of its own.
#'
#' Expected numbers of events, the mean, variance and covariance of the
#' log-rank statistic over calendar time, and the probabilities of group
#' sequential designs with efficacy and futility bounds are obtained from
#' closed-form expressions, Gauss-Legendre quadrature and multivariate normal
#' probabilities.
#'
#' A design may test one endpoint, progression-free survival (PFS) or overall
#' survival (OS), or both with a fixed-sequence or Bonferroni procedure; the
#' Bonferroni procedure may pass the level of a rejected hypothesis on to the
#' other (Maurer and Bretz 2013).
#' With two endpoints, the analyses may be driven by the events of either
#' endpoint or by calendar time, and the joint distribution of the PFS and
#' OS statistics gives the power of each endpoint, the probability of
#' rejecting both or at least one hypothesis, and the familywise error rate.
#'
#' @section Workflow:
#' \enumerate{
#'   \item Define each group with \code{\link{ms_model}} (or
#'     \code{\link{ms_fleischer}}) and the accrual with
#'     \code{\link{ms_accrual}}, and combine them with
#'     \code{\link{ms_scenario}}.
#'   \item Define the analysis times (and, with two endpoints, the endpoints
#'     tested at each analysis) with \code{\link{ms_looks}}, and the test and
#'     the multiple testing procedure with \code{\link{ms_test}}.
#'   \item Compute the power of a given design with \code{\link{ms_power}} or
#'     solve for the number of events, the sample size, the accrual duration
#'     or the follow-up with \code{\link{ms_size}}.
#' }
#' The building blocks \code{\link{ms_surv}}, \code{\link{ms_prob}},
#' \code{\link{ms_events}}, \code{\link{ms_event_time}},
#' \code{\link{lr_moments}}, \code{\link{lr_cov}}, \code{\link{gs_spend}},
#' \code{\link{gs_bounds}}, \code{\link{gs_prob}} and
#' \code{\link{gs_prob_mtp}} are exported so that each step can be checked
#' separately.
#'
#' @references
#' Maurer W, Bretz F (2013). Multiple testing in group sequential trials
#' using graphical approaches. Statistics in Biopharmaceutical Research
#' 5:311-320.
#'
#' @keywords internal
#' @importFrom stats pnorm qnorm uniroot
#' @importFrom mvtnorm pmvnorm
"_PACKAGE"

# Cache for quadrature rules
.dcs_env <- new.env(parent = emptyenv())
