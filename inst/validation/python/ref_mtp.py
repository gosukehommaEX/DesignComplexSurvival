"""Reference computations for the rejection probabilities of two endpoints (gs_prob_mtp).

Written independently of the R code, which sums rectangle probabilities over the outcomes of the first endpoint.
The functions are:

* recursion(): numerical integration over the analyses of the joint density of the two statistics, with the
  decisions of the procedure applied on the sample path (Gauss-Legendre rules on the pieces between the bounds).
  It requires the canonical joint distribution: the score processes of the two endpoints have independent
  bivariate normal increments between analyses.
* simulate(): Monte Carlo simulation of the procedure.
* glimm_two_look(): the formulas of Glimm et al. (2010, equations (1) and (2)) for two looks and equal
  information fractions, written as one-dimensional integrals by conditioning (Tamhane et al. 2010, Section 3.2
  and equation (8)).
* tamhane_fwer(): the familywise error rate of Tamhane et al. (2010, Section 3.2) under H2, used by
  tools/reproduce_python.py (project folder) for the reproduction of their Table 1.
"""
import numpy as np
from scipy import integrate, special, stats

QUAD = dict(epsabs=1e-13, epsrel=1e-11, limit=200)
NODES = 40
_GL400 = np.polynomial.legendre.leggauss(400)
WIDTH = 10.0


def canonical(cfg):
    """Mean vector and covariance matrix of the tested statistics (first endpoint, then second)."""
    t1 = np.asarray(cfg["tau1"], float)
    t2 = np.asarray(cfg["tau2"], float)
    d1 = np.diff(np.concatenate([[0.0], t1]))
    d2 = np.diff(np.concatenate([[0.0], t2]))
    cinc = cfg["rho"] * np.sqrt(d1 * d2)
    a1 = np.asarray(cfg["a1"]) - 1
    a2 = np.asarray(cfg["a2"]) - 1
    m = np.concatenate([cfg["theta1"] * np.sqrt(t1[a1]), cfg["theta2"] * np.sqrt(t2[a2])])
    k1, k2 = len(a1), len(a2)
    s = np.zeros((k1 + k2, k1 + k2))
    for p in range(k1):
        for q in range(k1):
            s[p, q] = np.sqrt(min(t1[a1[p]], t1[a1[q]]) / max(t1[a1[p]], t1[a1[q]]))
    for p in range(k2):
        for q in range(k2):
            s[k1 + p, k1 + q] = np.sqrt(min(t2[a2[p]], t2[a2[q]]) / max(t2[a2[p]], t2[a2[q]]))
    for p in range(k1):
        for q in range(k2):
            c = cinc[:min(a1[p], a2[q]) + 1].sum() / np.sqrt(t1[a1[p]] * t2[a2[q]])
            s[p, k1 + q] = s[k1 + q, p] = c
    return m, s


def _gl(a, b, n):
    t, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * t + 0.5 * (a + b), 0.5 * (b - a) * w


def _grid(center, cuts, labels):
    """Nodes, weights and region labels on [center - WIDTH, center + WIDTH] split at the cuts."""
    lo, hi = center - WIDTH, center + WIDTH
    pts = [lo] + [c for c in cuts if lo < c < hi] + [hi]
    xs, ws, ls = [], [], []
    for a, b in zip(pts[:-1], pts[1:]):
        x, w = _gl(a, b, NODES)
        mid = 0.5 * (a + b)
        xs.append(x); ws.append(w); ls.append(np.full(NODES, labels(mid), dtype=object))
    return np.concatenate(xs), np.concatenate(ws), np.concatenate(ls)


def recursion(cfg):
    """Rejection probabilities by recursive numerical integration over the analyses."""
    L = len(cfg["tau1"])
    t1 = np.asarray(cfg["tau1"], float)
    t2 = np.asarray(cfg["tau2"], float)
    rho, th1, th2 = cfg["rho"], cfg["theta1"], cfg["theta2"]
    a1 = list(np.asarray(cfg["a1"]) - 1)
    a2 = list(np.asarray(cfg["a2"]) - 1)
    k1, k2 = len(a1), len(a2)
    u1 = np.asarray(cfg["u1"], float)
    u2 = np.asarray(cfg["u2"], float)
    l1 = np.full(k1, -np.inf) if cfg.get("l1") is None else np.asarray(cfg["l1"], float)
    mtp, hier = cfg["mtp"], cfg.get("hierarchical", "overall")
    last1, last2 = a1[-1], a2[-1]
    p1 = np.zeros(k1); pf = np.zeros(k1); p2 = np.zeros(k2); both = 0.0
    phi = stats.norm.pdf
    A = B = C = D = None
    gx = gy = None
    for l in range(L):
        if l > max(last1, last2):
            break
        j = a1.index(l) if l in a1 else None
        i = a2.index(l) if l in a2 else None
        is_last1 = j is not None and j == k1 - 1
        # grids at analysis l
        if l <= last1:
            cuts, lab = [], (lambda m: "C")
            if j is not None:
                lo_b = l1[j] if (not is_last1 and np.isfinite(l1[j])) else -np.inf
                up_b = u1[j]
                cuts = [c for c in (lo_b, up_b) if np.isfinite(c)]
                lab = (lambda m, lo_b=lo_b, up_b=up_b: "R" if m >= up_b else ("F" if m < lo_b else "C"))
            nx = _grid(th1 * np.sqrt(t1[l]), cuts, lab)
        else:
            nx = None
        cuts, lab = [], (lambda m: "C")
        if i is not None and np.isfinite(u2[i]):
            cuts = [u2[i]]
            lab = (lambda m, b=u2[i]: "R" if m >= b else "C")
        ny = _grid(th2 * np.sqrt(t2[l]), cuts, lab)
        # densities before the decisions at analysis l
        if l == 0:
            X, Y = np.meshgrid(nx[0], ny[0], indexing="ij")
            m1, m2 = th1 * np.sqrt(t1[0]), th2 * np.sqrt(t2[0])
            Ap = stats.multivariate_normal(mean=[m1, m2], cov=[[1, rho], [rho, 1]]).pdf(
                np.stack([X, Y], -1))
            Bp = Cp = Dp = None
        else:
            dt1, dt2 = t1[l] - t1[l - 1], t2[l] - t2[l - 1]
            s1 = np.sqrt(dt1 / t1[l]) if dt1 > 0 else None
            s2 = np.sqrt(dt2 / t2[l])
            mxf = lambda x: (np.sqrt(t1[l - 1]) * x + th1 * dt1) / np.sqrt(t1[l])
            myf = lambda y: (np.sqrt(t2[l - 1]) * y + th2 * dt2) / np.sqrt(t2[l])
            Ap = None
            if A is not None and nx is not None:
                cc = rho * np.sqrt(dt1 * dt2) / np.sqrt(t1[l] * t2[l])
                kap = cc / s1 ** 2
                sc = np.sqrt(s2 ** 2 - cc ** 2 / s1 ** 2)
                g = A * gx[1][:, None] * gy[1][None, :]
                mx = mxf(gx[0]); my = myf(gy[0])
                Ap = np.zeros((nx[0].size, ny[0].size))
                for a in range(g.shape[0]):
                    if not np.any(g[a]):
                        continue
                    dx = nx[0] - mx[a]
                    kx = phi(dx / s1) / s1
                    arg = (ny[0][None, :, None] - my[None, None, :] - kap * dx[:, None, None]) / sc
                    Ap += kx[:, None] * ((phi(arg) / sc) @ g[a])
            def tr1(f, grid_prev, grid_new, mf, s):
                if f is None or grid_new is None:
                    return None
                gm = f * grid_prev[1]
                return (phi((grid_new[0][:, None] - mf(grid_prev[0])[None, :]) / s) / s) @ gm
            Bp = tr1(B, gy, ny, myf, s2)
            Dp = tr1(D, gy, ny, myf, s2)
            Cp = tr1(C, gx, nx, mxf, s1) if s1 is not None else None
        # decisions at analysis l
        yR = (ny[2] == "R") if i is not None else np.zeros(ny[0].size, bool)
        yC = ~yR
        An = Bn = Cn = Dn = None

        def add(old, new):
            return new if old is None else old + new

        if nx is not None:
            xR = nx[2] == "R"; xC = nx[2] == "C"; xF = nx[2] == "F"
        if Ap is not None:
            Wm = Ap * nx[1][:, None] * ny[1][None, :]
            if j is not None:
                p1[j] += Wm[xR].sum()
                pf[j] += Wm[xF].sum()
            if mtp == "bonferroni":
                if i is not None:
                    p2[i] += Wm[:, yR].sum()
                both += Wm[np.ix_(xR, yR)].sum()
                Bn = add(Bn, (Ap[xR] * nx[1][xR][:, None]).sum(0) * yC)
                if is_last1:
                    Dn = add(Dn, (Ap[xC] * nx[1][xC][:, None]).sum(0) * yC)
                else:
                    Cn = add(Cn, (Ap[:, yR] * ny[1][yR][None, :]).sum(1) * xC)
                    An = Ap * xC[:, None] * yC[None, :]
            else:
                q = Wm[np.ix_(xR, yR)].sum()
                if i is not None:
                    p2[i] += q
                both += q
                if hier == "overall":
                    Bn = add(Bn, (Ap[xR] * nx[1][xR][:, None]).sum(0) * yC)
                if not is_last1:
                    An = Ap * xC[:, None]
        if Bp is not None:
            if i is not None:
                q = (Bp * ny[1])[yR].sum()
                p2[i] += q; both += q
            Bn = add(Bn, Bp * yC)
        if Dp is not None:
            if i is not None:
                p2[i] += (Dp * ny[1])[yR].sum()
            Dn = add(Dn, Dp * yC)
        if Cp is not None:
            if j is not None:
                q = (Cp * nx[1])[xR].sum()
                p1[j] += q; both += q
                pf[j] += (Cp * nx[1])[xF].sum()
                if not is_last1:
                    Cn = add(Cn, Cp * xC)
            else:
                Cn = add(Cn, Cp)
        A, B, C, D = An, Bn, Cn, Dn
        gx, gy = nx, ny
    return dict(p1=p1, pf=pf, p2=p2, first=p1.sum(), second=p2.sum(), both=both,
                any=p1.sum() + p2.sum() - both)


def simulate(cfg, n, seed, mean=None, sigma=None):
    """Monte Carlo simulation of the procedure (decisions applied analysis by analysis).

    The statistics have the canonical joint distribution of cfg, or the given mean vector and covariance matrix.
    """
    m, s = canonical(cfg) if mean is None else (np.asarray(mean, float), np.asarray(sigma, float))
    rng = np.random.default_rng(seed)
    z = rng.multivariate_normal(m, s, size=n, method="cholesky")
    a1 = list(cfg["a1"]); a2 = list(cfg["a2"])
    k1, k2 = len(a1), len(a2)
    u1 = np.asarray(cfg["u1"], float); u2 = np.asarray(cfg["u2"], float)
    l1 = np.full(k1, -np.inf) if cfg.get("l1") is None else np.asarray(cfg["l1"], float)
    mtp, hier = cfg["mtp"], cfg.get("hierarchical", "overall")
    rej1 = np.full(n, -1); rej2 = np.full(n, -1); fut = np.full(n, -1)
    running = np.ones(n, bool)
    for an in sorted(set(a1) | set(a2)):
        j = a1.index(an) if an in a1 else None
        i = a2.index(an) if an in a2 else None
        prev1 = rej1 >= 0
        if j is not None:
            x = z[:, j]
            test1 = running & (rej1 < 0)
            new1 = test1 & (x >= u1[j])
            rej1[new1] = j
            stop_f = test1 & (x < l1[j]) & (j < k1 - 1)
        else:
            stop_f = np.zeros(n, bool)
        if i is not None:
            y = z[:, k1 + i]
            if mtp == "bonferroni":
                test2 = running & (rej2 < 0)
            elif hier == "overall":
                test2 = running & (rej2 < 0) & (rej1 >= 0)
            else:
                test2 = running & (rej2 < 0) & (rej1 >= 0) & ~prev1
                if j is None:
                    test2[:] = False
            rej2[test2 & (y >= u2[i])] = i
        fut[stop_f] = j if j is not None else -1
        stop = stop_f | ((rej1 >= 0) & (rej2 >= 0))
        if mtp == "fixed_sequence" and hier == "stagewise":
            stop |= rej1 >= 0
        if mtp == "fixed_sequence" and j is not None and j == k1 - 1:
            stop |= rej1 < 0
        running &= ~stop
    p1 = np.array([(rej1 == j).mean() for j in range(k1)])
    p2 = np.array([(rej2 == i).mean() for i in range(k2)])
    pf = np.array([(fut == j).mean() for j in range(k1)])
    both = ((rej1 >= 0) & (rej2 >= 0)).mean()
    se = lambda p: np.sqrt(p * (1 - p) / n)
    return dict(p1=p1, pf=pf, p2=p2, first=p1.sum(), second=p2.sum(), both=both,
                any=((rej1 >= 0) | (rej2 >= 0)).mean(), se=se)


def glimm_two_look(up, us, dp, ds, rho, t):
    """Glimm et al. (2010): R_a(H_s) (equation (1)), R_b(H_s) (equation (2)) and the primary rejection
    probabilities at the two looks, for equal information fractions t at the interim analysis.

    dp, ds: means of the final statistics (delta_p, delta_s). Conditioning on one statistic makes the other two
    conditionally independent (Tamhane et al. 2010, Section 3.2 and equation (8)).
    """
    Phi = stats.norm.cdf; phi = stats.norm.pdf
    g = np.sqrt(t)
    mp1, mp2, ms1, ms2 = g * dp, dp, g * ds, ds
    sr = np.sqrt(1 - rho ** 2); sg = np.sqrt(1 - t)
    # P(Zp1 >= up1, Zs1 >= us1): condition on U = Zp1 - mp1
    a1 = integrate.quad(lambda u: Phi((ms1 + rho * u - us[0]) / sr) * phi(u), up[0] - mp1, np.inf, **QUAD)[0]
    # P(Zp1 < up1, Zp2 >= up2, Zs2 >= us2): condition on U = Zp2 - mp2
    a2 = integrate.quad(lambda u: Phi((up[0] - mp1 - g * u) / sg) * Phi((ms2 + rho * u - us[1]) / sr) * phi(u),
                        up[1] - mp2, np.inf, **QUAD)[0]
    # P(Zp1 >= up1, Zs1 < us1, Zs2 >= us2): condition on V = Zs1 - ms1
    b3 = integrate.quad(lambda v: stats.norm.sf((up[0] - mp1 - rho * v) / sr)
                        * stats.norm.sf((us[1] - ms2 - g * v) / sg) * phi(v), -np.inf, us[0] - ms1, **QUAD)[0]
    # primary rejection at the two looks
    q1 = stats.norm.sf(up[0] - mp1)
    q2 = integrate.quad(lambda u: Phi((up[0] - mp1 - g * u) / sg) * phi(u), up[1] - mp2, np.inf, **QUAD)[0]
    return dict(a_look1=a1, a_look2=a2, a=a1 + a2, b_extra=b3, b=a1 + a2 + b3, p_look1=q1, p_look2=q2)


def tamhane_fwer(D1, rho, c, d, gamma=np.sqrt(0.5)):
    """FWER under H2 of Tamhane et al. (2010, Section 3.2) for the stagewise procedure with two looks.

    D1 is the mean of the first primary statistic (Delta_1 = Delta_11), a number or an array; c = (c1, c2),
    d = (d1, d2). For rho = 1 the secondary statistics equal the primary statistics minus their means. The integrals
    over u are computed with a Gauss-Legendre rule of 400 nodes on [lower limit, 12] (the normal density is below
    1e-31 beyond 12).
    """
    Phi = special.ndtr
    phi = lambda u: np.exp(-0.5 * u * u) / np.sqrt(2 * np.pi)
    D1 = np.atleast_1d(np.asarray(D1, float))
    D11, D12 = D1, D1 / gamma
    sg = np.sqrt(1 - gamma ** 2)
    t, w = _GL400

    def integral(lo, fun):
        lo = np.minimum(lo, 12.0)
        h = 0.5 * (12.0 - lo)
        u = h[:, None] * (t[None, :] + 1) + lo[:, None]
        return (h[:, None] * w[None, :] * fun(u) * phi(u)).sum(1)

    if rho >= 1:
        t1 = stats.norm.sf(np.maximum(c[0], d[0] + D11) - D11)
        lo = np.maximum(c[1], d[1] + D12) - D12
        out = t1 + integral(lo, lambda u: Phi((c[0] - D11[:, None] - gamma * u) / sg))
    else:
        sr = np.sqrt(1 - rho ** 2)
        f1 = integral(c[0] - D11, lambda u: Phi((-d[0] + rho * u) / sr))
        f2 = integral(c[1] - D12, lambda u: Phi((c[0] - D11[:, None] - gamma * u) / sg) * Phi((-d[1] + rho * u) / sr))
        out = f1 + f2
    return out if out.size > 1 else float(out[0])
