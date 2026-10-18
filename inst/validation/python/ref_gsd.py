"""Reference implementation of spending functions, bounds and stopping
probabilities of group sequential designs.

Written independently of the R code. Rectangle probabilities of the
multivariate normal distribution are computed by conditioning on the first
coordinate and integrating with adaptive quadrature (exact up to the
quadrature tolerance), not with the Miwa algorithm used in R.
"""

import numpy as np
from scipy import integrate, optimize, stats

QUAD = dict(epsabs=1e-13, epsrel=1e-11, limit=200)


def spend(t, alpha, kind, param=None):
    t = np.asarray(t, dtype=float)
    if kind == "OF":
        z = stats.norm.ppf(1 - alpha / 2)
        with np.errstate(divide="ignore"):
            return 2 * stats.norm.sf(z / np.sqrt(t))
    if kind == "Pocock":
        return alpha * np.log(1 + (np.e - 1) * t)
    if kind == "HSD":
        g = param
        if abs(g) < 1e-12:
            return alpha * t
        return alpha * (1 - np.exp(-g * t)) / (1 - np.exp(-g))
    if kind == "power":
        return alpha * t ** param
    raise ValueError(kind)


def mvn_box(lower, upper, mean, cov):
    """P(lower < X < upper) for X ~ N(mean, cov) by recursive conditioning."""
    lower = np.asarray(lower, float)
    upper = np.asarray(upper, float)
    mean = np.asarray(mean, float)
    cov = np.asarray(cov, float)
    k = lower.size
    s1 = np.sqrt(cov[0, 0])
    m1 = mean[0]
    if k == 1:
        return stats.norm.cdf((upper[0] - m1) / s1) - stats.norm.cdf((lower[0] - m1) / s1)
    b = cov[1:, 0] / cov[0, 0]
    rc = cov[1:, 1:] - np.outer(cov[1:, 0], cov[0, 1:]) / cov[0, 0]
    lo = max(lower[0], m1 - 12 * s1)
    hi = min(upper[0], m1 + 12 * s1)
    if hi <= lo:
        return 0.0

    def f(z):
        return (stats.norm.pdf((z - m1) / s1) / s1
                * mvn_box(lower[1:], upper[1:], mean[1:] + b * (z - m1), rc))

    return integrate.quad(f, lo, hi, **QUAD)[0]


def bounds(info_frac, alpha=0.025, sf_eff="OF", param_eff=None, futility="none",
           fut_bound=None, sf_fut="OF", param_fut=None, beta=0.1, binding=False,
           mean=None, corr=None, sigma_alt=None):
    t = np.asarray(info_frac, float)
    K = t.size
    if corr is None:
        corr = np.sqrt(np.minimum.outer(t, t) / np.maximum.outer(t, t))
    corr = np.asarray(corr, float)
    if sigma_alt is None:
        sigma_alt = corr
    sigma_alt = np.asarray(sigma_alt, float)
    a_inc = np.diff(np.concatenate([[0.0], spend(t, alpha, sf_eff, param_eff)]))
    lower = np.full(K, -np.inf)
    upper = np.full(K, np.inf)
    if futility == "fixed":
        lower[:K - 1] = fut_bound
    if futility == "beta":
        b_inc = np.diff(np.concatenate([[0.0], spend(t, beta, sf_fut, param_fut)]))
    for k in range(K):
        lp = lower[:k] if binding else np.full(k, -np.inf)
        if k == 0:
            upper[0] = stats.norm.isf(a_inc[0])
        else:
            f = lambda c: mvn_box(np.append(lp, c), np.append(upper[:k], np.inf),
                                  np.zeros(k + 1), corr[:k + 1, :k + 1]) - a_inc[k]
            upper[k] = optimize.brentq(f, -10, 10, xtol=1e-13)
        if k < K - 1:
            if futility == "beta":
                mu = np.asarray(mean[:k + 1], float)
                sg = sigma_alt[:k + 1, :k + 1]
                g = lambda a: mvn_box(np.append(lower[:k], -np.inf), np.append(upper[:k], a),
                                      mu, sg) - b_inc[k]
                if g(upper[k]) <= 0:
                    lower[k] = upper[k]
                else:
                    lower[k] = optimize.brentq(g, -15, upper[k], xtol=1e-13)
            if lower[k] > upper[k]:
                lower[k] = upper[k]
    lower[K - 1] = upper[K - 1]
    return lower, upper


def stop_probs(lower, upper, mean, cov):
    K = len(upper)
    pu = np.zeros(K)
    pl = np.zeros(K)
    mean = np.asarray(mean, float)
    for k in range(K):
        sg = cov[:k + 1, :k + 1]
        mu = mean[:k + 1]
        if np.isfinite(upper[k]):
            pu[k] = mvn_box(np.append(lower[:k], upper[k]), np.append(upper[:k], np.inf), mu, sg)
        if np.isfinite(lower[k]):
            pl[k] = mvn_box(np.append(lower[:k], -np.inf), np.append(upper[:k], lower[k]), mu, sg)
    return pu, pl
