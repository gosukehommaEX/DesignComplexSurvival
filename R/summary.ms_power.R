#' Summarize an ms_power object
#'
#' @param object An \code{ms_power} object.
#' @param ... Not used.
#'
#' @return An object of class \code{summary.ms_power}: a list with the
#'   design overview (\code{overview}), the per-look table (\code{looks}) and
#'   the settings (\code{settings}).
#'
#' @examples
#' ctl <- ms_model(death = log(2) / 12)
#' trt <- ms_model(death = 0.7 * log(2) / 12)
#' scn <- ms_scenario(ctl, trt, accrual = ms_accrual(rate = 20))
#' res <- ms_power(scn, n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
#'                 test = ms_test(), events = 300)
#' summary(res)
#'
#' @export
summary.ms_power <- function(object, ...) {
  overview <- c(n = object$n, accrual_duration = object$accrual_duration,
                study_duration = object$study_duration, events = object$events,
                power = object$power, type1_nonbinding = object$type1[["nonbinding"]],
                type1_binding = object$type1[["binding"]])
  structure(list(overview = overview, looks = object$looks,
                 settings = object$settings),
            class = "summary.ms_power")
}
