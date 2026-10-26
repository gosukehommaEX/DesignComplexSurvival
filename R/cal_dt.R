#' Time derivatives of the group functions at the nodes of a look
#'
#' For the look \code{lk} of \code{.cal_look()} at calendar time t, the
#' derivatives with respect to t of p_j G_j(s) Y_j(t, s) and
#' p_j G_j(s) D_j(t, s) at the nodes s of its rule: differentiating the
#' upper limit t - s of the integral over the entry time gives the OS
#' survival function and density of a patient enrolled at t - s (who can
#' switch from the follow-up time max(tc - t + s, 0) on) times the accrual
#' density there. Used for looks driven by OS events (\code{.cal_zeta()})
#' and for the times of given numbers of OS events
#' (\code{.cal_event_times()}).
#'
#' @param scenario An \code{ms_scenario} object with \code{.is_calendar()}.
#' @param acc Output of \code{.acc_setup()}.
#' @param lk Output of \code{.cal_look()}.
#' @return A list with \code{dy} and \code{dd}, lists by group (control,
#'   treatment) of the derivatives at the nodes.
#' @noRd
.cal_dt <- function(scenario, acc, lk) {
  tc <- scenario$switch_start
  a <- acc$alpha
  p <- scenario$prob
  models <- list(scenario$control, scenario$treatment)
  s <- lk$g$x
  e <- lk$t - s
  inside <- e < acc$duration
  al <- numeric(length(s))
  al[inside] <- .pw_value(e[inside], a[-length(a)], acc$rho)
  o <- pmax(tc - e, 0)
  dy <- vector("list", 2L)
  dd <- vector("list", 2L)
  for (j in 1:2) {
    os <- .cal_os(models[[j]], s, o)
    dy[[j]] <- p[[j]] * lk$G[[j]] * os$surv * al
    dd[[j]] <- p[[j]] * lk$G[[j]] * os$dens * al
  }
  list(dy = dy, dd = dd)
}
