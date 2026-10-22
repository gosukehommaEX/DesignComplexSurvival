# Design V2_K5_T2_alternative of tools/miwa_accuracy.R in the project folder:
# scenario B (scn_b()) with n = 300, five looks at 30, 60, 90, 120 and 150 OS
# events, O'Brien-Fleming type efficacy bounds, non-binding futility bounds by
# beta spending (Hwang-Shih-DeCani, gamma = -2, beta = 0.1) and
# variance = "alternative". The information fractions, the mean and the
# covariance of Z computed by the package, and as limits the bounds of
# version 0.2.0, all rounded to 10 significant digits. Expected values computed
# independently by inst/validation/python/expected_pmvn.py (results/pmvn.csv):
# exact recursion for a Markov control plus a quasi-Monte Carlo correction,
# with standard errors of at most 2.3e-9.
#   gs_prob_*: stopping probabilities at the limits lower and upper
#   bounds_*: bounds of gs_bounds() with these inputs, the stopping
#     probabilities at these bounds and the power
pmvn_v2k5 <- list(
  info_frac = c(0.2, 0.4, 0.6, 0.8, 1),
  mean = c(0.03136500376, 0.1636847573, 0.4297142694, 0.9535298331, 1.416778275),
  sigma = matrix(c(1.001408579, 0.7031236609, 0.5674164188, 0.4841443532, 0.4270453018,
                   0.7031236609, 1.004578143, 0.8117024007, 0.6913800771, 0.6085326028,
                   0.5674164188, 0.8117024007, 1.007772618, 0.8570058855, 0.7525497376,
                   0.4841443532, 0.6913800771, 0.8570058855, 1.006900231, 0.8811053347,
                   0.4270453018, 0.6085326028, 0.7525497376, 0.8811053347, 0.9915623001), 5),
  lower = c(-2.393271892, -2.031302568, -1.523812071, -0.7462330574, 2.031032013),
  upper = c(4.876884949, 3.357012328, 2.680280071, 2.289816764, 2.031032013),
  gs_prob_upper = c(6.423351420599231e-07, 0.0007207284932924066, 0.01190135340444286,
                    0.07970886468088609, 0.18008500596565455),
  gs_prob_lower = c(0.007697924225430667, 0.01148395322603465, 0.017132045208887847,
                    0.025632022361803, 0.6656374600984258),
  bounds_upper = c(4.876884948790766, 3.3570119215746455, 2.6802800645449927,
                   2.289816767712357, 2.031032043491168),
  bounds_lower = c(-2.393271891685609, -2.031302559783144, -1.5238120639071102,
                   -0.7472325815608224, 2.031032043491168),
  bounds_prob_upper = c(6.423351427359951e-07, 0.0007207295039655468,
                        0.011901352891859004, 0.07970886370709977, 0.18008499710979017),
  bounds_prob_lower = c(0.007697924232087868, 0.011483953489999895,
                        0.017132045457208218, 0.02555800884360692, 0.6657114824292399),
  bounds_power = 0.27241658554785725)
