"""Exact probabilities of rectangles over both endpoints for tests/testthat/test-pmvn_genz.R.

The rectangles come from designs with two endpoints of scenario F with n = 400 (D1: PFS rejected at look 2 and OS at
look 3 of design D1 of helper-designs.R; F6: three looks of each endpoint at 0.4, 0.7 and 1 of 220 OS events,
studentized variance, both continued at looks 1 and 2 and rejected at look 3). The parameters are those written in
test-pmvn_genz.R. Methods, independent of the R code:
* recursive conditioning with adaptive quadrature (ref_gsd.mvn_box()) for D1 and for the three-coordinate rectangle of
  the gs_prob_mtp() test;
* quasi-Monte Carlo with the separation of variables of Genz (1992) and 24 independent scramblings of 2^20 Sobol points
  for F6 (and for D1 as a check), with three standard errors over the scramblings.
Output: results/pmvn_genz.csv (about 3 minutes).

Usage: python3 expected_genz.py
"""
import csv
import os

import numpy as np
from scipy import special
from scipy.stats import qmc

import ref_gsd

INF = np.inf
D1 = dict(m=[2.141636271, 2.748018643, 2.034474618, 2.609699421],
          S=[1.0669695172, 0.8174879172, 0.3779545106, 0.2913956653,
             0.8174879172, 1.0538823953, 0.4843409699, 0.3728319907,
             0.3779545106, 0.4843409699, 1.0880711119, 0.8376117154,
             0.2913956653, 0.3728319907, 0.8376117154, 1.0767925675],
          lo=[-INF, 1.980965029, -INF, 1.981224804], hi=[2.668630143, INF, 2.664375831, INF])
F6 = dict(m=[2.3266834358, 2.91543213425, 3.21678539962, 1.66187457073, 2.19259585036, 2.60969942072],
          S=[1.0629926776, 0.832669469667, 0.748564986909, 0.502349318637, 0.378467094341, 0.316332203304,
             0.832669469667, 1.04791190754, 0.934245939148, 0.395728002013, 0.473693925547, 0.395540252994,
             0.748564986909, 0.934245939148, 1.02402167357, 0.355455009638, 0.42309965197, 0.432352726809,
             0.502349318637, 0.395728002013, 0.355455009638, 1.08808430447, 0.818750399697, 0.683601996344,
             0.378467094341, 0.473693925547, 0.42309965197, 0.818750399697, 1.08404799311, 0.903700658043,
             0.316332203304, 0.395540252994, 0.432352726809, 0.683601996344, 0.903700658043, 1.07679256746],
          lo=[-INF, -INF, 2.02976766455, -INF, -INF, 2.00053906808],
          hi=[2.93857575243, 2.25655974097, INF, 3.35686935609, 2.44454228117, INF])


def arrays(b):
    m = np.array(b["m"], float)
    return m, np.array(b["S"], float).reshape(len(m), len(m)), np.array(b["lo"], float), np.array(b["hi"], float)


def genz_integrand(lo, hi, L, u):
    """Separation of variables (Genz 1992) for zero-mean limits lo, hi and Cholesky factor L."""
    n = u.shape[0]
    d = L.shape[0]
    y = np.zeros((n, d))
    a = np.full(n, special.ndtr(lo[0] / L[0, 0]))
    b = np.full(n, special.ndtr(hi[0] / L[0, 0]))
    f = b - a
    for i in range(1, d):
        y[:, i - 1] = special.ndtri(np.clip(a + u[:, i - 1] * (b - a), 1e-300, 1 - 1e-16))
        s = y[:, :i] @ L[i, :i]
        a = special.ndtr((lo[i] - s) / L[i, i])
        b = special.ndtr((hi[i] - s) / L[i, i])
        f = f * (b - a)
    return f


def qmc_box(b, log2n=20, scramblings=24):
    m, S, lo, hi = arrays(b)
    sd = np.sqrt(np.diag(S))
    # most constrained coordinate first
    order = np.argsort(special.ndtr((hi - m) / sd) - special.ndtr((lo - m) / sd))
    L = np.linalg.cholesky(S[np.ix_(order, order)])
    lo0, hi0 = (lo - m)[order], (hi - m)[order]
    est = []
    for r in range(scramblings):
        sob = qmc.Sobol(len(m) - 1, scramble=True, seed=1000 + r)
        tot = 0.0
        for c in range(2 ** max(log2n - 18, 0)):
            u = sob.random_base2(min(log2n, 18)) if c == 0 else sob.random(2 ** 18)
            tot += genz_integrand(lo0, hi0, L, u).sum()
        est.append(tot / 2 ** log2n)
    est = np.array(est)
    return est.mean(), 3 * est.std(ddof=1) / np.sqrt(scramblings)


if __name__ == "__main__":
    out = []
    m, S, lo, hi = arrays(D1)
    out.append(("D1", "recursive_quadrature", ref_gsd.mvn_box(lo, hi, m, S), ""))
    i = [0, 2, 3]
    # gs_prob_mtp() test: first endpoint rejected at look 1, second continued at its look 1 and rejected at look 2
    out.append(("D1_first_look1", "recursive_quadrature",
                ref_gsd.mvn_box([D1["hi"][0], -INF, D1["lo"][3]], [INF, D1["hi"][2], INF], m[i], S[np.ix_(i, i)]),
                ""))
    for name, b in (("D1", D1), ("F6", F6)):
        v, e = qmc_box(b)
        out.append((name, "quasi_monte_carlo", v, e))
        print(name, repr(v), e, flush=True)
    os.makedirs("results", exist_ok=True)
    with open(os.path.join("results", "pmvn_genz.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["rectangle", "method", "value", "three_se"])
        for r in out:
            w.writerow([r[0], r[1], repr(float(r[2])), "" if r[3] == "" else "%.3g" % r[3]])
    for r in out:
        print(r)
