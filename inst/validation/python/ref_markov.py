"""Reference stopping probabilities of group sequential designs for covariance
matrices that need not be Markov.

Written independently of the R code (which uses mvtnorm::pmvnorm()).

- markov_stop(): exact stopping probabilities when the covariance matrix is
  that of a Markov sequence (tridiagonal inverse, as for statistics with
  independent increments), by recursion over the looks on composite
  Gauss-Legendre grids.
- cv_stop(): stopping probabilities for any covariance matrix: the exact
  values for a Markov control (by default the Markov completion, with the
  same variances and covariances of neighbouring looks) plus the difference
  between the two covariance matrices, estimated by the separation of
  variables of Genz (1992) with scrambled Sobol points (the same points for
  both matrices) and 8 independent scramblings, which give the standard
  error.
- reference(): markov_stop() for a Markov covariance, otherwise cv_stop()
  with 2^18 points per scrambling, multiplied by 4 until the largest standard
  error is at most se_target (at most 2^22).
- gl_stop(): tensor-product Gauss-Legendre quadrature in the Cholesky
  coordinates, deterministic, for at most 4 looks (used to check cv_stop()).
- bounds(): efficacy bounds by alpha spending (null correlation, Markov) and
  futility bounds by beta spending (any covariance) with the rules of
  gs_bounds().

The reference was checked in tools/check_miwa_accuracy.py (designs with 4
to 8 looks): against gl_stop() with 4 looks and, for Markov matrices, with
another Markov control against the exact recursion.
"""

import numpy as np
from scipy import optimize, special
from scipy.stats import qmc

import ref_gsd

ndtr = special.ndtr
ndtri = special.ndtri


def is_markov(cov, tol=1e-8):
    """Whether the inverse of the correlation matrix is tridiagonal (relative tolerance)."""
    cov = np.asarray(cov, float)
    k = cov.shape[0]
    sd = np.sqrt(np.diag(cov))
    prec = np.linalg.inv(cov / np.outer(sd, sd))
    off = max([abs(prec[i, j]) for i in range(k) for j in range(k) if abs(i - j) >= 2] + [0.0])
    return off <= tol * np.abs(prec).max(), off


def _gl_grid(lo, hi, h, nodes=12):
    x0, w0 = np.polynomial.legendre.leggauss(nodes)
    npan = max(1, int(np.ceil((hi - lo) / h)))
    edges = np.linspace(lo, hi, npan + 1)
    half = np.diff(edges) / 2
    mid = (edges[:-1] + edges[1:]) / 2
    return ((mid[:, None] + half[:, None] * x0[None, :]).ravel(),
            (half[:, None] * w0[None, :]).ravel())


def markov_stop(lower, upper, mean, cov, width=14.0, nodes=12, pan=0.5):
    """Stopping probabilities (upper, lower) at each look for a Markov covariance.

    The sub-density of continuing to look k is propagated on composite
    Gauss-Legendre grids over the continuation interval (truncated at the
    marginal mean +/- 14 standard deviations), with panels of half the
    smallest conditional standard deviation and 12 nodes per panel.
    """
    lower = np.asarray(lower, float)
    upper = np.asarray(upper, float)
    mean = np.asarray(mean, float)
    cov = np.asarray(cov, float)
    k_max = upper.size
    sd = np.sqrt(np.diag(cov))
    bet = np.zeros(k_max)
    cs = np.zeros(k_max)
    cs[0] = sd[0]
    for k in range(1, k_max):
        bet[k] = cov[k, k - 1] / cov[k - 1, k - 1]
        cs[k] = np.sqrt(cov[k, k] - cov[k, k - 1] ** 2 / cov[k - 1, k - 1])
    pu = np.zeros(k_max)
    pl = np.zeros(k_max)
    pu[0] = ndtr(-(upper[0] - mean[0]) / sd[0])
    pl[0] = ndtr((lower[0] - mean[0]) / sd[0])
    h = pan * cs[1:].min() if k_max > 1 else 1.0
    x = w = g = None
    for k in range(k_max - 1):
        lo = max(lower[k], mean[k] - width * sd[k])
        hi = min(upper[k], mean[k] + width * sd[k])
        if hi <= lo:
            break
        xk, wk = _gl_grid(lo, hi, h, nodes)
        if k == 0:
            gk = np.exp(-0.5 * ((xk - mean[0]) / sd[0]) ** 2) / (np.sqrt(2 * np.pi) * sd[0])
        else:
            m = mean[k] + bet[k] * (x - mean[k - 1])
            z = (xk[:, None] - m[None, :]) / cs[k]
            gk = (np.exp(-0.5 * z ** 2) / (np.sqrt(2 * np.pi) * cs[k])) @ (w * g)
        x, w, g = xk, wk, gk
        m = mean[k + 1] + bet[k + 1] * (x - mean[k])
        pu[k + 1] = np.sum(w * g * ndtr(-(upper[k + 1] - m) / cs[k + 1]))
        pl[k + 1] = np.sum(w * g * ndtr((lower[k + 1] - m) / cs[k + 1]))
    return pu, pl


def _markov_from(sd, r):
    k_max = sd.size
    out = np.diag(sd ** 2)
    for j in range(k_max):
        for k in range(j + 1, k_max):
            out[j, k] = out[k, j] = sd[j] * sd[k] * np.prod(r[j:k])
    return out


def markov_completion(cov):
    """Markov covariance with the same variances and neighbouring covariances."""
    cov = np.asarray(cov, float)
    sd = np.sqrt(np.diag(cov))
    r = np.array([cov[k, k + 1] / (sd[k] * sd[k + 1]) for k in range(cov.shape[0] - 1)])
    return _markov_from(sd, r)


def markov_shrunk(cov, factor=0.998):
    """Markov covariance with the variances of cov and its neighbouring
    correlations multiplied by factor (used to check cv_stop())."""
    cov = np.asarray(cov, float)
    sd = np.sqrt(np.diag(cov))
    r = np.array([cov[k, k + 1] / (sd[k] * sd[k + 1]) for k in range(cov.shape[0] - 1)])
    return _markov_from(sd, r * factor)


def _trunc(lo, hi, w):
    """N(0,1) probability of (lo, hi) and the inverse-cdf point at w in (lo, hi).

    For intervals in the upper half the reflected interval (-hi, -lo) is used,
    so that the tail probabilities keep their precision.
    """
    up = lo > 0
    dl = ndtr(np.where(up, -lo, lo))
    dh = ndtr(np.where(up, -hi, hi))
    width = np.maximum(np.where(up, dl - dh, dh - dl), 0.0)
    u = np.where(up, dl - w * (dl - dh), dl + w * (dh - dl))
    t = ndtri(np.clip(u, 1e-300, 1.0 - 1e-16))
    return width, np.clip(np.where(up, -t, t), lo, hi)


def _genz_terms(w, a, b, chol):
    """Integrands of the stopping probabilities of all looks at the points w."""
    nn = w.shape[0]
    k_max = a.size
    y = np.zeros((nn, k_max))
    f = np.ones(nn)
    pu = np.zeros((nn, k_max))
    pl = np.zeros((nn, k_max))
    for k in range(k_max):
        s = y[:, :k] @ chol[k, :k] if k > 0 else np.zeros(nn)
        lo = (a[k] - s) / chol[k, k]
        hi = (b[k] - s) / chol[k, k]
        pu[:, k] = f * ndtr(-hi)
        pl[:, k] = f * ndtr(lo)
        if k < k_max - 1:
            width, yk = _trunc(lo, hi, w[:, k])
            f = f * width
            y[:, k] = yk
    return pu, pl


def cv_stop(lower, upper, mean, cov, log2n=18, n_rand=8, seed=1, chunk=1 << 16, control=None):
    """Stopping probabilities with a Markov control variate.

    The control is the Markov completion of cov unless another Markov
    covariance is given. Returns (prob_upper, prob_lower) and their standard
    errors.
    """
    lower = np.asarray(lower, float)
    upper = np.asarray(upper, float)
    mean = np.asarray(mean, float)
    cov = np.asarray(cov, float)
    k_max = upper.size
    covm = markov_completion(cov) if control is None else np.asarray(control, float)
    ch = np.linalg.cholesky(cov)
    chm = np.linalg.cholesky(covm)
    a = lower - mean
    b = upper - mean
    diff = np.zeros((n_rand, 2, k_max))
    npts = 1 << log2n
    for r in range(n_rand):
        sob = qmc.Sobol(d=max(k_max - 1, 1), scramble=True, rng=np.random.default_rng([seed, r]))
        pts = sob.random_base2(log2n)
        for c0 in range(0, npts, chunk):
            w = pts[c0:c0 + chunk]
            pu, pl = _genz_terms(w, a, b, ch)
            qu, ql = _genz_terms(w, a, b, chm)
            diff[r, 0] += (pu - qu).sum(axis=0)
            diff[r, 1] += (pl - ql).sum(axis=0)
    diff /= npts
    mu, ml = markov_stop(lower, upper, mean, covm)
    d = diff.mean(axis=0)
    se = diff.std(axis=0, ddof=1) / np.sqrt(n_rand)
    return (mu + d[0], ml + d[1]), (se[0], se[1])


def reference(lower, upper, mean, cov, se_target=2.5e-8, log2n=18, log2n_max=22):
    """Reference stopping probabilities, standard errors and the method."""
    mk, _ = is_markov(cov)
    if mk:
        pu, pl = markov_stop(lower, upper, mean, cov)
        z = np.zeros(len(upper))
        return (pu, pl), (z, z.copy()), "markov"
    while True:
        (pu, pl), (su, sl) = cv_stop(lower, upper, mean, cov, log2n=log2n)
        if max(su.max(), sl.max()) <= se_target or log2n >= log2n_max:
            break
        log2n += 2
    return (pu, pl), (su, sl), "markov_cv 2^%d" % log2n


def gl_stop(lower, upper, mean, cov, nodes=160, width=9.0):
    """Stopping probabilities by tensor-product Gauss-Legendre quadrature (K <= 4).

    Integrates over the Cholesky coordinates y of the looks 1 to K - 1
    (X = mean + L y with y standard normal; each y truncated to its interval
    given the earlier coordinates and to +/- 9), with the normal probability
    of the last coordinate of each rectangle in closed form.
    """
    lower = np.asarray(lower, float)
    upper = np.asarray(upper, float)
    mean = np.asarray(mean, float)
    cov = np.asarray(cov, float)
    k_max = upper.size
    if k_max > 4:
        raise ValueError("gl_stop() is for at most 4 looks")
    chol = np.linalg.cholesky(cov)
    a = lower - mean
    b = upper - mean
    pu = np.zeros(k_max)
    pl = np.zeros(k_max)
    x0, w0 = np.polynomial.legendre.leggauss(nodes)
    y = np.zeros((1, 0))
    w = np.ones(1)
    for k in range(k_max):
        s = y @ chol[k, :k] if k > 0 else np.zeros(len(w))
        lo = (a[k] - s) / chol[k, k]
        hi = (b[k] - s) / chol[k, k]
        pu[k] = np.sum(w * ndtr(-hi))
        pl[k] = np.sum(w * ndtr(lo))
        if k == k_max - 1:
            break
        lo_c = np.maximum(lo, -width)
        hi_c = np.minimum(hi, width)
        half = np.where(hi_c > lo_c, (hi_c - lo_c) / 2, 0.0)
        mid = (hi_c + lo_c) / 2
        yk = mid[:, None] + half[:, None] * x0[None, :]
        wk = w[:, None] * half[:, None] * w0[None, :] * np.exp(-0.5 * yk ** 2) / np.sqrt(2 * np.pi)
        y = np.concatenate([np.repeat(y, nodes, axis=0), yk.reshape(-1, 1)], axis=1)
        w = wk.ravel()
        keep = w > 0
        y = y[keep]
        w = w[keep]
    return pu, pl


def bounds(info_frac, mean, sigma_alt, alpha=0.025, sf_eff="OF", param_eff=None, sf_fut="HSD",
           param_fut=-2.0, beta=0.1, binding=False, log2n=20):
    """Efficacy bounds by alpha spending under the null correlation
    sqrt(t_j / t_k) (Markov, exact recursion) and futility bounds by beta
    spending under (mean, sigma_alt), with the rules of gs_bounds(): a
    futility bound above the efficacy bound is set to the efficacy bound, and
    at the last look the futility bound is the efficacy bound. For a
    covariance that is not Markov, the probabilities of the futility bounds
    are those of cv_stop() with 2^log2n points per scrambling (the same
    points at every evaluation, so that the function of the bound is
    smooth). Roots by brentq (xtol 1e-12). Returns the futility and efficacy
    bounds and the largest standard error of the beta spent at the
    futility bounds."""
    t = np.asarray(info_frac, float)
    mean = np.asarray(mean, float)
    sigma_alt = np.asarray(sigma_alt, float)
    k_max = t.size
    corr = np.sqrt(np.minimum.outer(t, t) / np.maximum.outer(t, t))
    a_inc = np.diff(np.concatenate([[0.0], ref_gsd.spend(t, alpha, sf_eff, param_eff)]))
    b_inc = np.diff(np.concatenate([[0.0], ref_gsd.spend(t, beta, sf_fut, param_fut)]))
    lower = np.full(k_max, -np.inf)
    upper = np.full(k_max, np.inf)
    se_max = 0.0
    for k in range(k_max):
        lp = lower[:k].copy() if binding else np.full(k, -np.inf)

        def fe(c):
            lo = np.append(lp, -np.inf)
            up = np.append(upper[:k], c)
            return markov_stop(lo, up, np.zeros(k + 1), corr[:k + 1, :k + 1])[0][k] - a_inc[k]

        upper[k] = optimize.brentq(fe, 0.0, 12.0, xtol=1e-12)
        if k == k_max - 1:
            break
        sd_k = np.sqrt(sigma_alt[k, k])
        if k == 0:
            lower[0] = min(mean[0] + sd_k * special.ndtri(b_inc[0]), upper[0])
            continue
        sub = sigma_alt[:k + 1, :k + 1]
        markov = is_markov(sub)[0]

        def ff(x):
            lo = np.append(lower[:k], x)
            if markov:
                pl = markov_stop(lo, upper[:k + 1], mean[:k + 1], sub)[1]
                sl = np.zeros(k + 1)
            else:
                (_, pl), (_, sl) = cv_stop(lo, upper[:k + 1], mean[:k + 1], sub, log2n=log2n)
            ff.se = max(ff.se, sl[k])
            return pl[k] - b_inc[k]

        ff.se = 0.0
        if ff(upper[k]) <= 0:
            lower[k] = upper[k]
        else:
            lower[k] = optimize.brentq(ff, mean[k] - 15 * sd_k, upper[k], xtol=1e-12)
        se_max = max(se_max, ff.se)
    lower[k_max - 1] = upper[k_max - 1]
    return lower, upper, se_max
