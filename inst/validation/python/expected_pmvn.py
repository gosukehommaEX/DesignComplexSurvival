"""Expected values for the tests of rectangle probabilities with a covariance
matrix that is not Markov (tests/testthat/test-pmvn.R, test-is_markov.R,
test-gs_prob.R, test-gs_bounds.R and test-ms_power.R).

Case V2K5: design V2_K5_T2_alternative of tools/miwa_accuracy.R (scenario B
of the tests with n = 300, five looks at 30, 60, 90, 120 and 150 OS events,
O'Brien-Fleming type efficacy bounds, non-binding futility bounds by beta
spending with Hwang-Shih-DeCani gamma = -2 and beta = 0.1,
variance = "alternative"). The mean and covariance of Z are those computed
by the package (tools/logs/miwa_accuracy_*.csv), rounded to 10 significant
digits; the limits of the gs_prob() case are the bounds of version 0.2.0,
rounded in the same way. Values (ref_markov.py, independent of the R code):
- gs_prob: stopping probabilities at the given limits (reference());
- bounds: efficacy and futility bounds of gs_bounds() with these inputs
  (ref_markov.bounds()) and the stopping probabilities at these bounds;
- is_markov: largest element of the inverse correlation matrix outside the
  tridiagonal band, relative to the largest element.

Case CAP: gs_bounds(c(0.25, 0.5, 0.75, 1), sf_eff = "Pocock",
futility = "beta", sf_fut = "HSD", param_fut = -2, beta = 0.1,
binding = TRUE, mean = c(2, 5, 7, 9)) with sigma_alt equal to the null
correlation, and the same with binding = FALSE: the bounds of look 1, the
efficacy bound of look 2 and the beta spent at look 2 by a futility bound
equal to it, which is below the beta to spend at look 2 (b_inc2), so that
the futility bound of look 2 is set to the efficacy bound (exact Markov
recursion).

Cases K9 and K10: rectangles of 9 and 10 coordinates of tools/miwa_accuracy.R
(looks at equal information fractions with the correlation of independent
increments, mean 2.8 sqrt(t), limits and means rounded to 10 significant
digits; tests/testthat/helper-pmvn.R). K9: continuation intervals at the
first eight looks and crossing of the efficacy bound at the ninth; K10:
crossing of the efficacy bound at the tenth look without futility bounds.
Exact probabilities by the Markov recursion; for K9 also the value that the
inclusion-exclusion over the limits gives without the terms of nine or more
coordinates (miwa_without_terms), which is what the Miwa algorithm of
mvtnorm returns.

Output: results/pmvn.csv (a few minutes).

Usage: python3 expected_pmvn.py
"""
import csv
import itertools
import os

import numpy as np
from scipy import optimize, special

import ref_gsd
import ref_markov as RM

INF = np.inf
T5 = np.array([0.2, 0.4, 0.6, 0.8, 1.0])
M5 = np.array([0.03136500376, 0.1636847573, 0.4297142694, 0.9535298331, 1.416778275])
S5 = np.array([
    [1.001408579, 0.7031236609, 0.5674164188, 0.4841443532, 0.4270453018],
    [0.7031236609, 1.004578143, 0.8117024007, 0.6913800771, 0.6085326028],
    [0.5674164188, 0.8117024007, 1.007772618, 0.8570058855, 0.7525497376],
    [0.4841443532, 0.6913800771, 0.8570058855, 1.006900231, 0.8811053347],
    [0.4270453018, 0.6085326028, 0.7525497376, 0.8811053347, 0.9915623001]])
LO5 = np.array([-2.393271892, -2.031302568, -1.523812071, -0.7462330574, 2.031032013])
UP5 = np.array([4.876884949, 3.357012328, 2.680280071, 2.289816764, 2.031032013])
LO9 = np.array([-1.066666667, -0.6800673418, -0.3834192463, -0.1333333333, 0.086996779,
                0.2861904266, 0.4693678903, 0.6398653164, 2.074017588])
UP9 = np.array([6.622531106, 4.61271216, 3.712056021, 3.178855076, 2.821696693, 2.562388399,
                2.363415838, 2.204586649, INF])
M9 = np.array([0.9333333333, 1.319932658, 1.616580754, 1.866666667, 2.086996779, 2.286190427,
               2.46936789, 2.639865316, 2.8])
LO10 = np.array([-INF] * 9 + [2.081175108])
UP10 = np.array([6.991351707, 4.87680179, 3.929682519, 3.367079025, 2.989329537, 2.714808617,
                 2.504076879, 2.335828578, 2.197502745, INF])
M10 = np.array([0.8854377448, 1.252198067, 1.533623161, 1.77087549, 1.979898987, 2.168870674,
                2.342648074, 2.504396135, 2.656313235, 2.8])


def corr_of(k):
    t = np.arange(1, k + 1) / k
    return np.sqrt(np.minimum.outer(t, t) / np.maximum.outer(t, t))


def big_cases():
    c9 = corr_of(9)
    # upper crossing at look 9 after the continuation intervals of looks 1-8
    exact9 = RM.markov_stop(LO9, np.append(UP9[:-1], LO9[-1]), M9, c9)[0][-1]
    # inclusion-exclusion with the infinite limit replaced by mean + 1000
    # (as in .pmvn()): sum over the subsets S of coordinates set to their
    # lower limits of (-1)^|S| P(X < c_S), here the terms of |S| >= 9
    hi9 = np.where(np.isinf(UP9), M9 + 1000.0, UP9)
    omitted = 0.0
    for size in range(9, 10):
        for sub in itertools.combinations(range(9), size):
            c = hi9.copy()
            c[list(sub)] = LO9[list(sub)]
            lo_last = np.full(9, -INF)
            lo_last[-1] = c[-1]
            omitted += (-1) ** size * RM.markov_stop(lo_last, c, M9, c9)[1][-1]
    c10 = corr_of(10)
    exact10 = RM.markov_stop(LO10, np.append(UP10[:-1], LO10[-1]), M10, c10)[0][-1]
    return dict(k9_exact=exact9, k9_miwa_without_terms=exact9 - omitted, k10_exact=exact10)


def cap_case():
    t = np.array([0.25, 0.5, 0.75, 1.0])
    mu = np.array([2.0, 5.0, 7.0, 9.0])
    corr = np.sqrt(np.minimum.outer(t, t) / np.maximum.outer(t, t))
    a_inc = np.diff(np.concatenate([[0.0], ref_gsd.spend(t, 0.025, "Pocock")]))
    b_inc = np.diff(np.concatenate([[0.0], ref_gsd.spend(t, 0.1, "HSD", -2.0)]))
    u1 = special.ndtri(1 - a_inc[0])
    l1 = min(mu[0] + special.ndtri(b_inc[0]), u1)
    out = dict(u1=u1, l1=l1, b_inc2=b_inc[1])
    for tag, lp in (("binding", l1), ("nonbinding", -INF)):
        def fe(c):
            return RM.markov_stop([lp, -INF], [u1, c], [0.0, 0.0], corr[:2, :2])[0][1] - a_inc[1]

        u2 = optimize.brentq(fe, 0.0, 12.0, xtol=1e-12)
        # beta spent at look 2 if the futility bound were the efficacy bound
        # (the futility bound is set to the efficacy bound when this is at
        # most b_inc2)
        out["u2_" + tag] = u2
        out["beta_at_u2_" + tag] = RM.markov_stop([l1, u2], [u1, u2], mu[:2], corr[:2, :2])[1][1]
    return out


if __name__ == "__main__":
    rows = []
    mk, off = RM.is_markov(S5)
    prec = np.linalg.inv(S5 / np.outer(np.sqrt(np.diag(S5)), np.sqrt(np.diag(S5))))
    rows.append(("V2K5", "is_markov_offband_relative", 0, off / np.abs(prec).max(), ""))
    mk3, off3 = RM.is_markov(S5[:3, :3])
    prec3 = np.linalg.inv(S5[:3, :3] / np.outer(np.sqrt(np.diag(S5[:3, :3])), np.sqrt(np.diag(S5[:3, :3]))))
    rows.append(("V2K5", "is_markov_offband_relative_3", 0, off3 / np.abs(prec3).max(), ""))
    print("is_markov", mk, off, mk3, off3, flush=True)
    (pu, pl), (su, sl), meth = RM.reference(LO5, UP5, M5, S5, se_target=1e-9)
    print("gs_prob", meth, pu, pl, su.max(), sl.max(), flush=True)
    for k in range(5):
        rows.append(("V2K5", "gs_prob_upper", k + 1, pu[k], su[k]))
        rows.append(("V2K5", "gs_prob_lower", k + 1, pl[k], sl[k]))
    lo, up, se_b = RM.bounds(T5, M5, S5)
    print("bounds", lo, up, se_b, flush=True)
    (bu, bl), (bsu, bsl), meth = RM.reference(lo, up, M5, S5, se_target=1e-9)
    print("probs at bounds", meth, bu, bl, bu.sum(), flush=True)
    for k in range(5):
        rows.append(("V2K5", "bounds_upper", k + 1, up[k], ""))
        rows.append(("V2K5", "bounds_lower", k + 1, lo[k], se_b))
        rows.append(("V2K5", "bounds_prob_upper", k + 1, bu[k], bsu[k]))
        rows.append(("V2K5", "bounds_prob_lower", k + 1, bl[k], bsl[k]))
    rows.append(("V2K5", "bounds_power", 0, bu.sum(), float(np.sqrt((bsu ** 2).sum()))))
    cap = cap_case()
    print("cap", cap, flush=True)
    for key, v in cap.items():
        rows.append(("CAP", key, 0, v, ""))
    big = big_cases()
    print("big", big, flush=True)
    rows.append(("K9", "exact", 0, big["k9_exact"], ""))
    rows.append(("K9", "miwa_without_terms", 0, big["k9_miwa_without_terms"], ""))
    rows.append(("K10", "exact", 0, big["k10_exact"], ""))
    os.makedirs("results", exist_ok=True)
    with open(os.path.join("results", "pmvn.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["case", "quantity", "look", "value", "se"])
        for r in rows:
            w.writerow([r[0], r[1], r[2], repr(float(r[3])), "" if r[4] == "" else "%.3g" % r[4]])
    for r in rows:
        print(r)
