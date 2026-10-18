#' Exponential model for PFS and OS of Fleischer et al. (2009)
#'
#' Creates the \code{ms_model} in which the time to progression and the time
#' to death are exponential and the hazard of death does not change at
#' progression. Then OS is exponential with median \code{median_os} and PFS
#' is exponential with median \code{median_pfs}.
#'
#' @details
#' With death hazard \code{lambda_d = log(2) / median_os} and PFS hazard
#' \code{log(2) / median_pfs}, the progression hazard is
#' \code{log(2) / median_pfs - lambda_d}. This is the model of Section 2 of
#' Fleischer et al. (2009) with the hazard of death before and after
#' progression equal to \code{lambda_d}.
#'
#' @param median_pfs Median PFS.
#' @param median_os Median OS (not smaller than \code{median_pfs}).
#' @param dropout Dropout hazard.
#'
#' @return An object of class \code{ms_model}.
#'
#' @references
#' Fleischer F, Gaschler-Markefski B, Bluhmki E (2009). A statistical model
#' for the dependence between progression-free survival and overall
#' survival. Statistics in Medicine 28:2669-2686.
#'
#' @seealso \code{\link{ms_model}}
#'
#' @examples
#' m <- ms_fleischer(median_pfs = 6, median_os = 15)
#' ms_surv(m, time = c(6, 15), endpoint = "pfs")
#' ms_surv(m, time = c(6, 15), endpoint = "os")
#'
#' @export
ms_fleischer <- function(median_pfs, median_os, dropout = 0) {
  .check_positive(median_pfs, "median_pfs")
  .check_positive(median_os, "median_os")
  if (median_os < median_pfs) {
    stop("'median_os' must not be smaller than 'median_pfs'.", call. = FALSE)
  }
  death <- log(2) / median_os
  prog <- log(2) / median_pfs - death
  ms_model(prog = prog, death = death, death_prog = death, dropout = dropout)
}
