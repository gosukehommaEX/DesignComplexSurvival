#' Lines of the (s, o) plane where the OS given the availability is not smooth
#'
#' For a group that can switch, the OS survival function and density of
#' \code{.os_avail()} at the follow-up time s, with switching possible from
#' the follow-up time o on, are smooth in (s, o) except at the break points
#' s of \code{.scn_breaks()} and on lines o = c0 + c1 s. On these lines a
#' break point of the integral over the time of progression u that depends
#' on o meets one that does not, or two of them meet. The break points that
#' depend on o are u = o and, with the reset clock, u = o - b and, with the
#' AFT effect (factor phi) and the reset clock, u = o (1 - phi) + phi s - b,
#' for the break points b of \code{time_post}; the line o = s separates the
#' follow-up times before and after the availability. The lines are
#' \itemize{
#'   \item forward clock: o = s, o = a for the break points a of \code{time}
#'     and \code{time_post}, o = s - b (switch clock \code{"switch"}) and
#'     o = (b - phi s) / (1 - phi) (AFT effect);
#'   \item reset clock: o = a and o = a + b, o = s + b - b', o = s + b - x
#'     for the points x of \code{.os_conv()} (b + b' for the switch clock
#'     \code{"switch"}, (b - (1 - phi) b') / phi for the AFT effect), and with
#'     the AFT effect o = (b + c0 + (c1 - phi) s) / (1 - phi) for the break
#'     points c0 + c1 s in 0, a, s - b', s - x and s, and
#'     o = s + (b - b') / phi.
#' }
#'
#' @param model An \code{ms_model} object.
#' @return A matrix with columns \code{c0} and \code{c1}, one row per line
#'   (no duplicates); no rows for a model that cannot switch.
#' @noRd
.cal_lines <- function(model) {
  if (!.has_switch(model)) {
    return(matrix(numeric(0), 0L, 2L, dimnames = list(NULL, c("c0", "c1"))))
  }
  tm <- model$time
  tp <- model$time_post
  phi <- if (is.null(model$switch_af) || model$switch_af == 1) NULL else model$switch_af
  sw <- identical(model$switch_clock, "switch")
  ln <- function(c0, c1) cbind(c0 = c0, c1 = rep_len(c1, length(c0)))
  if (model$clock == "forward") {
    out <- rbind(ln(c(tm, tp), 0), ln(0, 1))
    if (sw) {
      out <- rbind(out, ln(-tp, 1))
    } else if (!is.null(phi)) {
      out <- rbind(out, ln(tp / (1 - phi), -phi / (1 - phi)))
    }
  } else {
    dd <- as.vector(outer(tp, tp, "-"))
    out <- rbind(ln(c(tm, as.vector(outer(tm, tp, "+"))), 0), ln(dd, 1))
    ex <- NULL
    if (sw) {
      ex <- as.vector(outer(tp, tp, "+"))
    } else if (!is.null(phi)) {
      ex <- as.vector(outer(tp, tp, function(x, y) (x - (1 - phi) * y) / phi))
    }
    if (!is.null(ex)) {
      out <- rbind(out, ln(as.vector(outer(tp, ex, "-")), 1))
    }
    if (!is.null(phi)) {
      cs <- rbind(ln(c(0, tm), 0), ln(-c(0, tp, ex), 1))
      out <- rbind(out,
                   ln(as.vector(outer(cs[, 1L], tp, "+")) / (1 - phi),
                      rep((cs[, 2L] - phi) / (1 - phi), length(tp))),
                   ln(dd / phi, 1))
    }
  }
  unique(out)
}
