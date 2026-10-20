# Evaluates 'code' with the absolute error tolerance 'eps' for the rectangles
# of gs_prob_mtp() computed by the Genz-Bretz algorithm (.pmvn_genz()), so
# that results can be compared with exact values at a tolerance of 1e-5
with_abseps <- function(code, eps = 1e-6) {
  op <- options(DesignComplexSurvival.mtp_abseps = eps)
  on.exit(options(op))
  force(code)
}
