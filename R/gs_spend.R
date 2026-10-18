#' Error spending functions
#'
#' Cumulative error spent at information fractions \code{t} for a total
#' (one-sided) error \code{alpha}.
#'
#' @details
#' The available types are
#' \itemize{
#'   \item \code{"OF"}: Lan-DeMets spending function of O'Brien-Fleming type,
#'     \code{2 - 2 pnorm(qnorm(1 - alpha / 2) / sqrt(t))}.
#'   \item \code{"Pocock"}: Lan-DeMets spending function of Pocock type,
#'     \code{alpha log(1 + (exp(1) - 1) t)}.
#'   \item \code{"HSD"}: Hwang-Shih-DeCani family with parameter
#'     \code{param} (gamma), \code{alpha (1 - exp(-gamma t)) / (1 - exp(-gamma))},
#'     and \code{alpha t} for gamma = 0.
#'   \item \code{"power"}: Kim-DeMets power family with parameter
#'     \code{param} (rho), \code{alpha t^rho}.
#'   \item \code{"user"}: cumulative proportions of \code{alpha} given in
#'     \code{param}, one value per element of \code{t}.
#' }
#' The same functions are used for beta spending, with \code{alpha} replaced
#' by the total type II error.
#'
#' @param t Information fractions between 0 and 1.
#' @param alpha Total error to spend.
#' @param type Type of spending function: \code{"OF"}, \code{"Pocock"},
#'   \code{"HSD"}, \code{"power"} or \code{"user"}.
#' @param param Parameter of the spending function (required for
#'   \code{"HSD"}, \code{"power"} and \code{"user"}).
#'
#' @return Numeric vector of cumulative error spent.
#'
#' @references
#' Lan KKG, DeMets DL (1983). Discrete sequential boundaries for clinical
#' trials. Biometrika 70:659-663.
#'
#' @seealso \code{\link{gs_bounds}}
#'
#' @examples
#' gs_spend(c(0.5, 0.75, 1), alpha = 0.025, type = "OF")
#' gs_spend(c(0.5, 0.75, 1), alpha = 0.1, type = "HSD", param = -2)
#'
#' @export
gs_spend <- function(t, alpha = 0.025,
                     type = c("OF", "Pocock", "HSD", "power", "user"),
                     param = NULL) {
  type <- match.arg(type)
  if (!is.numeric(t) || anyNA(t) || any(t < 0) || any(t > 1 + 1e-12)) {
    stop("'t' must contain values between 0 and 1.", call. = FALSE)
  }
  .check_prob(alpha, "alpha")
  t <- pmin(t, 1)
  if (type %in% c("HSD", "power", "user") && is.null(param)) {
    stop(sprintf("'param' must be given for type = \"%s\".", type), call. = FALSE)
  }
  switch(type,
    OF = 2 * stats::pnorm(stats::qnorm(1 - alpha / 2) / sqrt(t), lower.tail = FALSE),
    Pocock = alpha * log(1 + (exp(1) - 1) * t),
    HSD = {
      if (!is.numeric(param) || length(param) != 1L || !is.finite(param)) {
        stop("'param' must be a finite number for type = \"HSD\".", call. = FALSE)
      }
      if (abs(param) < 1e-12) alpha * t else alpha * (1 - exp(-param * t)) / (1 - exp(-param))
    },
    power = {
      if (!is.numeric(param) || length(param) != 1L || !is.finite(param) || param <= 0) {
        stop("'param' must be a positive number for type = \"power\".", call. = FALSE)
      }
      alpha * t^param
    },
    user = {
      if (!is.numeric(param) || length(param) != length(t) || anyNA(param) ||
          any(param < 0) || any(param > 1) || any(diff(param) < 0)) {
        stop(paste0("'param' must be nondecreasing cumulative proportions between 0 and 1, ",
                    "one per element of 't'."), call. = FALSE)
      }
      alpha * param
    }
  )
}
