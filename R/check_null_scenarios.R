#' Check the null scenarios of a design with two endpoints
#'
#' @param null_scenarios \code{NULL}, or a list of \code{ms_scenario} objects
#'   named by the endpoint whose hypothesis is true in the scenario.
#' @param scenario The design scenario.
#' @param endpoints The two tested endpoints.
#' @return \code{TRUE} invisibly; an error otherwise.
#' @noRd
.check_null_scenarios <- function(null_scenarios, scenario, endpoints) {
  if (is.null(null_scenarios)) {
    return(invisible(TRUE))
  }
  nm <- names(null_scenarios)
  if (!is.list(null_scenarios) || inherits(null_scenarios, "ms_scenario") ||
      length(null_scenarios) == 0L || is.null(nm) || !all(nm %in% endpoints) ||
      !all(vapply(null_scenarios, inherits, logical(1), what = "ms_scenario"))) {
    stop(paste0("'null_scenarios' must be a list of 'ms_scenario' objects, each named by ",
                "the endpoint (\"os\" or \"pfs\") whose hypothesis is true in it."),
         call. = FALSE)
  }
  same <- vapply(null_scenarios, function(s) {
    isTRUE(all.equal(s$prob, scenario$prob)) && isTRUE(all.equal(s$accrual, scenario$accrual))
  }, logical(1))
  if (!all(same)) {
    stop("The scenarios of 'null_scenarios' must have the allocation and accrual of 'scenario'.",
         call. = FALSE)
  }
  same_pfs <- vapply(null_scenarios, function(s) {
    identical(s$pfs_switch, scenario$pfs_switch)
  }, logical(1))
  if (!all(same_pfs)) {
    stop("The scenarios of 'null_scenarios' must have the 'pfs_switch' of 'scenario'.",
         call. = FALSE)
  }
  invisible(TRUE)
}
