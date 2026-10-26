#' Whether a group model can switch after progression
#'
#' @param model An \code{ms_model} object.
#' @return \code{TRUE} if the model has progression and a positive
#'   probability of switching at progression or a positive switch hazard after
#'   progression.
#' @noRd
.has_switch <- function(model) {
  any(model$prog > 0) && (any(model$switch_prob > 0) || any(model$switch_rate > 0))
}
