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

Output: results/pmvn.csv (a few minutes).

Usage: python3 expected_pmvn.py
"""
import csv
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
    os.makedirs("results", exist_ok=True)
    with open(os.path.join("results", "pmvn.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["case", "quantity", "look", "value", "se"])
        for r in rows:
            w.writerow([r[0], r[1], r[2], repr(float(r[3])), "" if r[4] == "" else "%.3g" % r[4]])
    for r in rows:
        print(r)
