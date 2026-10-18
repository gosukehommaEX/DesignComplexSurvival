#' Survival function, hazard and density of an endpoint
#'
#' Evaluates the distribution of progression-free survival or overall
#' survival implied by an \code{ms_model}, without dropout.
#'
#' @details
#' PFS has the hazard \code{prog + death}. OS combines death before
#' progression with death after progression, with or without switch; its
#' hazard is not piecewise constant when \code{prog} is positive.
#'
#' @param model An \code{ms_model} object.
#' @param time Times since randomization (nonnegative).
#' @param endpoint \code{"os"} or \code{"pfs"}.
#'
#' @return A data frame with \code{time}, the survival function
#'   \code{surv}, the hazard \code{haz} (\code{NA} where \code{surv} is 0)
#'   and the density \code{dens}.
#'
#' @seealso \code{\link{ms_model}}, \code{\link{ms_prob}}
#'
#' @examples
#' m <- ms_fleischer(median_pfs = 6, median_os = 15)
#' ms_surv(m, time = c(6, 15), endpoint = "os")
#'
#' @export
ms_surv <- function(model, time, endpoint = c("os", "pfs")) {
  endpoint <- match.arg(endpoint)
  if (!inherits(model, "ms_model")) {
    stop("'model' must be an 'ms_model' object.", call. = FALSE)
  }
  .check_times(time, "time")
  ep <- .ep_eval(model, time, endpoint)
  haz <- rep(NA_real_, length(time))
  pos <- ep$surv > 0
  haz[pos] <- ep$dens[pos] / ep$surv[pos]
  data.frame(time = time, surv = ep$surv, haz = haz, dens = ep$dens)
}
