"""Reference computations for the rejection probabilities of two endpoints under the Bonferroni-based graphical
procedure with recycling of the level (Maurer and Bretz 2013, Algorithm 1), for gs_prob_mtp() with upper_recycled.

Written independently of the R code, which sums rectangle probabilities over the outcomes of the procedure at the
analyses. Hypothesis k is tested with the bounds u_k at its looks until the other hypothesis is rejected and with the
bounds r_k (r_k <= u_k) afterwards; within an analysis, a hypothesis rejected with its initial bound makes the other
hypothesis be tested with its bound r at the same analysis. The trial stops when both hypotheses are rejected or when
the statistic of the first endpoint is below its futility bound l_1 at an interim look at which its hypothesis is
tested and not rejected. The functions are:

* recursion(): numerical integration over the analyses of the joint density of the two statistics, with the
  decisions of the procedure applied on the sample path (Gauss-Legendre rules on the pieces between the bounds). It
  requires the canonical joint distribution of ref_mtp.canonical(): the score processes of the two endpoints have
  independent bivariate normal increments between analyses.
* simulate(): Monte Carlo simulation of Algorithm 1 for two hypotheses (reject any hypothesis whose statistic is at
  or above the bound of its current level, update the levels, repeat within the analysis), with the stopping rules
  above.

A configuration has the fields of ref_mtp.py (tau1, tau2, rho, theta1, theta2, a1, a2, u1, u2, l1) and r1, r2.
"""
import numpy as np
from scipy import stats

import ref_mtp

NODES = 40
WIDTH = 10.0


def _gl(a, b, n):
    t, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * t + 0.5 * (a + b), 0.5 * (b - a) * w


def _grid(center, cuts):
    """Nodes, weights and the midpoint of the piece of each node on [center - WIDTH, center + WIDTH], split at the
    finite cuts."""
    lo, hi = center - WIDTH, center + WIDTH
    pts = [lo] + sorted(set(c for c in cuts if np.isfinite(c) and lo < c < hi)) + [hi]
    xs, ws, ms = [], [], []
    for a, b in zip(pts[:-1], pts[1:]):
        x, w = _gl(a, b, NODES)
        xs.append(x); ws.append(w); ms.append(np.full(NODES, 0.5 * (a + b)))
    return np.concatenate(xs), np.concatenate(ws), np.concatenate(ms)


def recursion(cfg):
    """Rejection and futility probabilities by recursive numerical integration over the analyses."""
    t1 = np.asarray(cfg["tau1"], float)
    t2 = np.asarray(cfg["tau2"], float)
    rho, th1, th2 = cfg["rho"], cfg["theta1"], cfg["theta2"]
    a1 = [int(a) - 1 for a in cfg["a1"]]
    a2 = [int(a) - 1 for a in cfg["a2"]]
    k1, k2 = len(a1), len(a2)
    u1, u2 = np.asarray(cfg["u1"], float), np.asarray(cfg["u2"], float)
    r1, r2 = np.asarray(cfg["r1"], float), np.asarray(cfg["r2"], float)
    l1 = np.full(k1, -np.inf) if cfg.get("l1") is None else np.asarray(cfg["l1"], float).copy()
    l1[k1 - 1] = -np.inf
    last1, last2 = a1[-1], a2[-1]
    p1 = np.zeros(k1); pf = np.zeros(k1); p2 = np.zeros(k2); both = 0.0
    phi = stats.norm.pdf
    # densities on the grids of the previous analysis, before its decisions were applied they were turned into:
    # A(x, y) neither hypothesis rejected; B(y) first rejected, second tested with r2; C(x) second rejected, first
    # tested with r1; D(y) first not rejected at its last look, second tested with u2
    A = B = C = D = None
    gx = gy = None
    for l in range(max(last1, last2) + 1):
        j = a1.index(l) if l in a1 else None
        i = a2.index(l) if l in a2 else None
        last_j = j is not None and j == k1 - 1
        nx = None
        if l <= last1:
            cx = [] if j is None else [u1[j], r1[j], l1[j]]
            nx = _grid(th1 * np.sqrt(t1[l]), cx)
        cy = [] if i is None else [u2[i], r2[i]]
        ny = _grid(th2 * np.sqrt(t2[l]), cy)
        # densities at analysis l, before its decisions
        if l == 0:
            X, Y = np.meshgrid(nx[0], ny[0], indexing="ij")
            Ap = stats.multivariate_normal(mean=[th1 * np.sqrt(t1[0]), th2 * np.sqrt(t2[0])],
                                           cov=[[1.0, rho], [rho, 1.0]]).pdf(np.stack([X, Y], -1))
            Bp = Cp = Dp = None
        else:
            dt1, dt2 = t1[l] - t1[l - 1], t2[l] - t2[l - 1]
            s2 = np.sqrt(dt2 / t2[l])
            mxf = lambda x: (np.sqrt(t1[l - 1]) * x + th1 * dt1) / np.sqrt(t1[l])
            myf = lambda y: (np.sqrt(t2[l - 1]) * y + th2 * dt2) / np.sqrt(t2[l])
            Ap = None
            if A is not None:
                s1 = np.sqrt(dt1 / t1[l])
                cc = rho * np.sqrt(dt1 * dt2) / np.sqrt(t1[l] * t2[l])
                kap = cc / s1 ** 2
                sc = np.sqrt(s2 ** 2 - cc ** 2 / s1 ** 2)
                g = A * gx[1][:, None] * gy[1][None, :]
                mx, my = mxf(gx[0]), myf(gy[0])
                Ap = np.zeros((nx[0].size, ny[0].size))
                for a in range(g.shape[0]):
                    if not np.any(g[a]):
                        continue
                    dx = nx[0] - mx[a]
                    kx = phi(dx / s1) / s1
                    arg = (ny[0][None, :, None] - my[None, None, :] - kap * dx[:, None, None]) / sc
                    Ap += kx[:, None] * ((phi(arg) / sc) @ g[a])

            def move(f, prev, new, mf, s):
                if f is None or new is None:
                    return None
                return (phi((new[0][:, None] - mf(prev[0])[None, :]) / s) / s) @ (f * prev[1])

            Bp = move(B, gy, ny, myf, s2)
            Dp = move(D, gy, ny, myf, s2)
            Cp = move(C, gx, nx, mxf, np.sqrt(dt1 / t1[l])) if C is not None else None
        # decisions at analysis l: flags on the pieces (False where the endpoint has no look)
        if nx is not None:
            xa = np.zeros(nx[0].size, bool) if j is None else nx[2] >= u1[j]
            xb = np.zeros(nx[0].size, bool) if j is None else nx[2] >= r1[j]
            xf = np.zeros(nx[0].size, bool) if j is None else nx[2] < l1[j]
        ya = np.zeros(ny[0].size, bool) if i is None else ny[2] >= u2[i]
        yb = np.zeros(ny[0].size, bool) if i is None else ny[2] >= r2[i]
        An = Bn = Cn = Dn = None

        def acc(old, new):
            return new if old is None else old + new

        if Ap is not None:
            W = Ap * nx[1][:, None] * ny[1][None, :]
            m_both = (xa[:, None] & yb[None, :]) | ((~xa & xb)[:, None] & ya[None, :])
            m_one = xa[:, None] & ~yb[None, :]
            m_two = (~xb)[:, None] & ya[None, :]
            m_none = (~xa)[:, None] & (~ya)[None, :]
            q = W[m_both].sum()
            if q > 0:
                p1[j] += q; p2[i] += q; both += q
            if j is not None:
                p1[j] += W[m_one].sum()
            Bn = acc(Bn, (Ap * nx[1][:, None] * m_one).sum(0))
            if i is not None:
                p2[i] += W[m_two].sum()
            if j is not None:
                pf[j] += W[m_two & xf[:, None]].sum() + W[m_none & xf[:, None]].sum()
            go2 = m_two & ~xf[:, None]
            go0 = m_none & ~xf[:, None]
            if not last_j:
                Cn = acc(Cn, (Ap * ny[1][None, :] * go2).sum(1))
                An = Ap * go0
            else:
                Dn = acc(Dn, (Ap * nx[1][:, None] * go0).sum(0))
        if Bp is not None:
            if i is not None:
                q = (Bp * ny[1])[yb].sum()
                p2[i] += q; both += q
            Bn = acc(Bn, Bp * ~yb)
        if Cp is not None:
            if j is not None:
                q = (Cp * nx[1])[xb].sum()
                p1[j] += q; both += q
                pf[j] += (Cp * nx[1])[~xb & xf].sum()
            if not last_j:
                Cn = acc(Cn, Cp * (~xb & ~xf))
        if Dp is not None:
            if i is not None:
                p2[i] += (Dp * ny[1])[ya].sum()
            Dn = acc(Dn, Dp * ~ya)
        A, B, C, D = An, Bn, Cn, Dn
        gx, gy = nx, ny
    return dict(p1=p1, pf=pf, p2=p2, first=p1.sum(), second=p2.sum(), both=both, any=p1.sum() + p2.sum() - both)


def simulate(cfg, n, seed, mean=None, sigma=None):
    """Monte Carlo simulation of Algorithm 1 of Maurer and Bretz (2013) for two hypotheses with the stopping rules.

    The statistics have the canonical joint distribution of cfg, or the given mean vector and covariance matrix.
    """
    m, s = ref_mtp.canonical(cfg) if mean is None else (np.asarray(mean, float), np.asarray(sigma, float))
    rng = np.random.default_rng(seed)
    z = rng.multivariate_normal(m, s, size=n, method="cholesky")
    an = [list(cfg["a1"]), list(cfg["a2"])]
    k1, k2 = len(an[0]), len(an[1])
    off = [0, k1]
    u = [np.asarray(cfg["u1"], float), np.asarray(cfg["u2"], float)]
    r = [np.asarray(cfg["r1"], float), np.asarray(cfg["r2"], float)]
    l1 = np.full(k1, -np.inf) if cfg.get("l1") is None else np.asarray(cfg["l1"], float)
    rej = [np.full(n, -1), np.full(n, -1)]
    fut = np.full(n, -1)
    running = np.ones(n, bool)
    for a in sorted(set(an[0]) | set(an[1])):
        look = [an[e].index(a) if a in an[e] else None for e in (0, 1)]
        changed = True
        while changed:
            changed = False
            for e in (0, 1):
                if look[e] is None:
                    continue
                bound = np.where(rej[1 - e] >= 0, r[e][look[e]], u[e][look[e]])
                new = running & (rej[e] < 0) & (z[:, off[e] + look[e]] >= bound)
                if new.any():
                    rej[e][new] = look[e]
                    changed = True
        if look[0] is not None and look[0] < k1 - 1:
            stop = running & (rej[0] < 0) & (z[:, look[0]] < l1[look[0]])
            fut[stop] = look[0]
            running &= ~stop
        running &= ~((rej[0] >= 0) & (rej[1] >= 0))
    p1 = np.array([(rej[0] == j).mean() for j in range(k1)])
    p2 = np.array([(rej[1] == i).mean() for i in range(k2)])
    pf = np.array([(fut == j).mean() for j in range(k1)])
    r1, r2 = rej[0] >= 0, rej[1] >= 0
    se = lambda p: np.sqrt(p * (1 - p) / n)
    return dict(p1=p1, pf=pf, p2=p2, first=r1.mean(), second=r2.mean(), both=(r1 & r2).mean(),
                any=(r1 | r2).mean(), se=se)
