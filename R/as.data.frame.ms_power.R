#' Per-look results of an ms_power object as a data frame
#'
#' @param x An \code{ms_power} object.
#' @param row.names,optional Not used.
#' @param ... Not used.
#'
#' @return The data frame \code{x$looks}.
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' res <- ms_power(scn, n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
#'                 test = ms_test(), events = 300)
#' as.data.frame(res)
#'
#' @export
as.data.frame.ms_power <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$looks
}
