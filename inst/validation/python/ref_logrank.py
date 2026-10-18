"""Reference implementation of expected events and log-rank moments.

Written independently of the R code. The accrual is applied through the
enrolled fraction Abar(t - s) inside each integrand, and all integrals are
evaluated by adaptive quadrature. The variance and covariance of the
log-rank score under the alternative follow the formulas of Luo et al.
(2019, Statistics in Medicine 38:192-209, Section 3.1) as printed, and an
independent derivation through per-patient influence functions is
implemented in `influence_cov()` for cross-checking.
"""

import numpy as np
from scipy import integrate, optimize

import ref_model as rm

QUAD_OPTS = dict(epsabs=1e-12, epsrel=1e-11, limit=400)


class Accrual:
    """Piecewise uniform accrual given by rates or by proportions."""

    def __init__(self, time=(0.0,), rate=None, prop=None, duration=None):
        self.time = np.asarray(time, dtype=float)
        self.rate = None if rate is None else np.asarray(rate, dtype=float)
        self.prop = None if prop is None else np.asarray(prop, dtype=float)
        self.duration = duration

    def setup(self, n):
        """Return breaks alpha, densities rho (fraction of n per unit time), A."""
        if self.rate is not None:
            # find A with cumulative rate equal to n
            def cum(x):
                return rm.pw_cum(x, self.time, self.rate)
            hi = 1.0
            while cum(hi) < n:
                hi *= 2.0
            A = optimize.brentq(lambda x: cum(x) - n, 0.0, hi, xtol=1e-14)
            alpha = [x for x in self.time if x < A] + [A]
            rho = [self.rate[j] / n for j in range(len(alpha) - 1)]
        else:
            A = float(self.duration)
            alpha = list(self.time) + [A]
            rho = [self.prop[j] / (alpha[j + 1] - alpha[j]) for j in range(len(alpha) - 1)]
        return np.asarray(alpha), np.asarray(rho), A


def abar(alpha, rho, x):
    """Enrolled fraction by calendar time x."""
    tot = 0.0
    for j in range(rho.size):
        tot += rho[j] * max(min(x, alpha[j + 1]) - alpha[j], 0.0)
    return tot


class Scenario:
    def __init__(self, control, treatment, allocation, accrual):
        self.control = control
        self.treatment = treatment
        self.p1 = allocation / (1.0 + allocation)
        self.p0 = 1.0 - self.p1
        self.accrual = accrual

    def breaks(self):
        return sorted(set(rm.model_breaks(self.control) + rm.model_breaks(self.treatment)))


def _points(scn, alpha, t):
    pts = [x for x in scn.breaks() if 0 < x < t]
    pts += [t - a for a in alpha if 0 < t - a < t]
    return sorted(set(pts))


def events(scn, n, t, endpoint):
    """Expected events and dropouts by group at calendar time t."""
    alpha, rho, A = scn.accrual.setup(n)
    pts = _points(scn, alpha, t)
    out = {}
    for g, m, p in (("control", scn.control, scn.p0), ("treatment", scn.treatment, scn.p1)):
        fe = lambda s, m=m: abar(alpha, rho, t - s) * rm.dens(m, s, endpoint) * rm.dropout_surv(m, s)
        fd = lambda s, m=m: (abar(alpha, rho, t - s) * rm.surv(m, s, endpoint)
                             * rm.dropout_surv(m, s) * rm.dropout_haz(m, s))
        out["events_" + g] = n * p * integrate.quad(fe, 0, t, points=pts or None, **QUAD_OPTS)[0]
        out["dropouts_" + g] = n * p * integrate.quad(fd, 0, t, points=pts or None, **QUAD_OPTS)[0]
    out["events"] = out["events_control"] + out["events_treatment"]
    out["enrolled"] = n * abar(alpha, rho, t)
    return out


def event_time(scn, n, target, endpoint):
    f = lambda t: events(scn, n, t, endpoint)["events"] - target
    hi = 1.0
    while f(hi) < 0:
        hi *= 2.0
    return optimize.brentq(f, 0.0, hi, xtol=1e-11)


def _group_funcs(scn, s, endpoint):
    """Per-group S, f, G, dropout hazard at follow-up time s."""
    res = []
    for m, p in ((scn.control, scn.p0), (scn.treatment, scn.p1)):
        res.append((p, rm.surv(m, s, endpoint), rm.dens(m, s, endpoint),
                    rm.dropout_surv(m, s), rm.dropout_haz(m, s)))
    return res


def moments(scn, n, t, endpoint):
    """Score mean, information and Z means (direct and Schoenfeld) at time t."""
    alpha, rho, A = scn.accrual.setup(n)
    pts = _points(scn, alpha, t)

    def parts(s):
        ab = abar(alpha, rho, t - s)
        (p0, S0, f0, G0, _), (p1, S1, f1, G1, _) = _group_funcs(scn, s, endpoint)
        y0 = p0 * S0 * G0 * ab
        y1 = p1 * S1 * G1 * ab
        if y0 + y1 <= 0:
            return 0.0, 0.0, 0.0, 0.0
        l0 = f0 / S0
        l1 = f1 / S1
        u = y0 * y1 / (y0 + y1) * (l1 - l0)
        v = y0 * y1 / (y0 + y1) ** 2 * (y0 * l0 + y1 * l1)
        dd = ab * (p0 * f0 * G0 + p1 * f1 * G1)
        lw = dd * np.log(l1 / l0) if dd > 0 else 0.0
        return u, v, dd, lw

    vals = []
    for k in range(4):
        vals.append(integrate.quad(lambda s: parts(s)[k], 0, t, points=pts or None, **QUAD_OPTS)[0])
    score = n * vals[0]
    info = n * vals[1]
    d = n * vals[2]
    lw = n * vals[3]
    return dict(score_mean=score, info=info, events=d,
                z_direct=-score / np.sqrt(info),
                z_schoenfeld=-np.sqrt(scn.p0 * scn.p1) * lw / np.sqrt(d))


# ---------------------------------------------------------------------------
# Luo et al. (2019): variance and covariance of the score under the alternative
# ---------------------------------------------------------------------------

class LuoFuncs:
    """Functions r_k, q_k, zeta_k, zeta, a(s) of Luo et al. (2019), per patient."""

    def __init__(self, scn, endpoint):
        self.scn = scn
        self.endpoint = endpoint

    def all(self, s):
        (p0, S0, f0, G0, _), (p1, S1, f1, G1, _) = _group_funcs(self.scn, s, self.endpoint)
        r0 = p0 * S0 * G0
        r1 = p1 * S1 * G1
        r = r0 + r1
        z0 = f0 * G0
        z1 = f1 * G1
        q1 = r1 / r
        q0 = 1.0 - q1
        zeta = p1 * z1 + p0 * z0
        a = p1 * q0 * z1 - p0 * q1 * z0
        return dict(p0=p0, p1=p1, r=r, q0=q0, q1=q1, z0=z0, z1=z1, zeta=zeta, a=a)


def luo_var(scn, n, t, endpoint):
    """Per-patient U(t), I(t), V(t) = I + V1 + V2 of Luo et al. (2019)."""
    alpha, rho, A = scn.accrual.setup(n)
    lf = LuoFuncs(scn, endpoint)
    pts = _points(scn, alpha, t)
    G = lambda x: abar(alpha, rho, x)
    q = lambda f, lo, hi, p: integrate.quad(f, lo, hi, points=p or None, **QUAD_OPTS)[0]

    U = q(lambda s: G(t - s) * lf.all(s)["a"], 0, t, pts)
    I = q(lambda s: G(t - s) * lf.all(s)["q1"] * lf.all(s)["q0"] * lf.all(s)["zeta"], 0, t, pts)
    V1 = q(lambda s: G(t - s) * (lf.all(s)["q0"] - lf.all(s)["q1"]) * lf.all(s)["a"], 0, t, pts) - U ** 2

    br = [x for x in scn.breaks() if x > 0]

    def inner(s, k):
        p = [x for x in br if x < s]
        key = "q1" if k == 1 else "q0"
        return q(lambda x: lf.all(x)[key] * lf.all(x)["zeta"] / lf.all(x)["r"], 0, s, p)

    def v2int(s):
        f = lf.all(s)
        return G(t - s) * (f["q1"] * inner(s, 1) - f["q0"] * inner(s, 0)) * f["a"]

    V2 = 2.0 * q(v2int, 0, t, pts)
    return dict(U=U, I=I, V=I + V1 + V2)


def luo_cov(scn, n, t, tp, endpoint):
    """Per-patient CV(t, t') = CV1 + CV2 of Luo et al. (2019) for t' <= t."""
    assert tp <= t
    alpha, rho, A = scn.accrual.setup(n)
    lf = LuoFuncs(scn, endpoint)
    G = lambda x: abar(alpha, rho, x)
    q = lambda f, lo, hi, p: integrate.quad(f, lo, hi, points=p or None, **QUAD_OPTS)[0]
    br = [x for x in scn.breaks() if x > 0]

    Ut = luo_var(scn, n, t, endpoint)["U"]
    Utp = luo_var(scn, n, tp, endpoint)["U"]

    pts_tp = _points(scn, alpha, tp)

    def cv1(s):
        f = lf.all(s)
        return G(tp - s) * (f["p1"] * f["q0"] ** 2 * f["z1"] + f["p0"] * f["q1"] ** 2 * f["z0"])

    CV1 = q(cv1, 0, tp, pts_tp) - Ut * Utp

    def b(x, k):
        f = lf.all(x)
        return (f["q1"] if k == 1 else f["q0"]) * f["zeta"] / f["r"]

    # region s' <= s: s in [0, t], s' in [0, min(s, t')]
    def outer1(s):
        f = lf.all(s)
        hi = min(s, tp)
        if hi <= 0:
            return 0.0
        p = [x for x in br if x < hi] + [tp - a for a in alpha if 0 < tp - a < hi]
        p += [s - (t - tp)] if 0 < s - (t - tp) < hi else []
        g1 = q(lambda x: G(min(t - s, tp - x)) * b(x, 1), 0, hi, sorted(set(p)))
        g0 = q(lambda x: G(min(t - s, tp - x)) * b(x, 0), 0, hi, sorted(set(p)))
        return (g1 * f["q1"] - g0 * f["q0"]) * f["a"]

    pts1 = _points(scn, alpha, t) + [tp]
    part1 = q(outer1, 0, t, sorted(set(x for x in pts1 if 0 < x < t)))

    # region s <= s': s' in [0, t'], s in [0, s']
    def outer2(sp):
        f = lf.all(sp)
        if sp <= 0:
            return 0.0
        p = [x for x in br if x < sp] + [t - a for a in alpha if 0 < t - a < sp]
        p += [sp + (t - tp)] if 0 < sp + (t - tp) < sp else []
        g1 = q(lambda x: G(min(t - x, tp - sp)) * b(x, 1), 0, sp, sorted(set(p)))
        g0 = q(lambda x: G(min(t - x, tp - sp)) * b(x, 0), 0, sp, sorted(set(p)))
        return (g1 * f["q1"] - g0 * f["q0"]) * f["a"]

    part2 = q(outer2, 0, tp, pts_tp)
    return CV1 + part1 + part2


def influence_cov(scn, n, t, tp, endpoint, grid_pts=None):
    """E[psi(t) psi(t')] - U(t) U(t') by per-patient influence functions.

    psi = delta g_Z(X) - L_Z(X) with g_1 = q0, g_0 = -q1 and
    L_z(x) = int_0^x g_z(s) lambda_bar(s) ds. The cumulative functions are
    obtained by solving an ODE system piecewise between model breaks.
    """
    from scipy.integrate import solve_ivp
    assert tp <= t
    alpha, rho, A = scn.accrual.setup(n)
    lf = LuoFuncs(scn, endpoint)
    br = [0.0] + [x for x in scn.breaks() if 0 < x < t] + [t]

    def rhs(x, y):
        f = lf.all(x)
        lbar = f["zeta"] / f["r"]
        (p0, S0, f0, G0, d0), (p1, S1, f1, G1, d1) = _group_funcs(scn, x, endpoint)
        out = []
        for z, (L, e, c) in enumerate(((y[0], f0 * G0, S0 * G0 * d0), (y[1], f1 * G1, S1 * G1 * d1))):
            g = -f["q1"] if z == 0 else f["q0"]
            out.append((g, L, e, c))
        dy = np.zeros(6)
        for z in range(2):
            g, L, e, c = out[z]
            dy[z] = g * lbar
            dy[2 + z] = (g - L) * e - L * c          # H_z
            dy[4 + z] = (g - L) ** 2 * e + L ** 2 * c  # J_z
        return dy

    sols = []
    y0 = np.zeros(6)
    for lo, hi in zip(br[:-1], br[1:]):
        sol = solve_ivp(rhs, (lo, hi), y0, method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True)
        sols.append((lo, hi, sol))
        y0 = sol.y[:, -1]

    def Y(x):
        for lo, hi, sol in sols:
            if lo <= x <= hi:
                return sol.sol(x)
        raise ValueError(x)

    def Rz(x):
        (p0, S0, f0, G0, _), (p1, S1, f1, G1, _) = _group_funcs(scn, x, endpoint)
        return S0 * G0, S1 * G1

    def K(e):
        W = t - e
        Wp = tp - e
        yW = Y(W)
        yWp = Y(Wp)
        R0, R1 = Rz(W)
        tot = 0.0
        for z, (p, Rw) in enumerate(((scn.p0, R0), (scn.p1, R1))):
            LW, LWp = yW[z], yWp[z]
            HW, HWp = yW[2 + z], yWp[2 + z]
            JWp = yWp[4 + z]
            tot += p * (JWp + (-LWp) * (HW - HWp) + LW * LWp * Rw)
        return tot

    hi = min(tp, A)
    pts = [a for a in alpha if 0 < a < hi]
    pts += [t - x for x in scn.breaks() if 0 < t - x < hi]
    pts += [tp - x for x in scn.breaks() if 0 < tp - x < hi]
    rho_at = lambda e: rho[min(int(np.searchsorted(alpha, e, side="right")) - 1, rho.size - 1)]
    E = integrate.quad(lambda e: rho_at(e) * K(e), 0, hi, points=sorted(set(pts)) or None, **QUAD_OPTS)[0]
    Ut = luo_var(scn, n, t, endpoint)["U"]
    Utp = luo_var(scn, n, tp, endpoint)["U"]
    return E - Ut * Utp
