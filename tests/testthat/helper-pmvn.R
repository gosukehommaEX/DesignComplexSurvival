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

# Rectangles of 9 and 10 coordinates of tools/miwa_accuracy.R in the project
# folder: looks at equal information fractions with the correlation of
# independent increments, mean 2.8 sqrt(t), O'Brien-Fleming type efficacy
# bounds and, for nine coordinates, continuation intervals from mean - 2 at
# the first eight looks and crossing of the efficacy bound at the ninth;
# for ten coordinates, crossing of the efficacy bound at the tenth look
# without futility bounds. Limits and means rounded to 10 significant
# digits. Exact probabilities by the Markov recursion of
# inst/validation/python/expected_pmvn.py (results/pmvn.csv, cases K9 and
# K10).
pmvn_corr <- function(k) {
  tt <- seq_len(k) / k
  sqrt(outer(tt, tt, pmin) / outer(tt, tt, pmax))
}
pmvn_k9 <- list(
  lower = c(-1.066666667, -0.6800673418, -0.3834192463, -0.1333333333, 0.086996779,
            0.2861904266, 0.4693678903, 0.6398653164, 2.074017588),
  upper = c(6.622531106, 4.61271216, 3.712056021, 3.178855076, 2.821696693, 2.562388399,
            2.363415838, 2.204586649, Inf),
  mean = c(0.9333333333, 1.319932658, 1.616580754, 1.866666667, 2.086996779, 2.286190427,
           2.46936789, 2.639865316, 2.8),
  corr = pmvn_corr(9),
  exact = 0.089061312288082)
pmvn_k10 <- list(
  lower = c(rep(-Inf, 9), 2.081175108),
  upper = c(6.991351707, 4.87680179, 3.929682519, 3.367079025, 2.989329537, 2.714808617,
            2.504076879, 2.335828578, 2.197502745, Inf),
  mean = c(0.8854377448, 1.252198067, 1.533623161, 1.77087549, 1.979898987, 2.168870674,
           2.342648074, 2.504396135, 2.656313235, 2.8),
  corr = pmvn_corr(10),
  exact = 0.08654380846621315)
