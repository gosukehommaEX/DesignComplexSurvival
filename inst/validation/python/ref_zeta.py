"""Reference computation of the moments behind the studentized log-rank covariance.

Written independently of the R code (Planning/covariance_note, Sections 3-4). For a patient the functionals of
one endpoint are r = delta rho(X) + B(X) with B the integral of beta:
  psi (score):          rho = g,          beta = -g lbar
  phi (observed var.):  rho = q (1 - q),  beta = (1 - 2 q) g lbar
  dlt (event):          rho = 1,          beta = 0
with q the proportion of the treatment group among patients at risk, lbar the pooled hazard among patients at
risk, g = 1 - q (treatment) or -q (control).

Same endpoint (two calendar times t <= t'): the cumulative functions L, K, A_ab and M_a are obtained by solving an
ODE system (scipy solve_ivp, DOP853) and the expectation over the entry time by adaptive quadrature; the R code
uses Gauss-Legendre rules.

PFS and OS (Fleischer model only: per group, progression A ~ Exp(a) and death D ~ Exp(b) independent, PFS =
min(A, D), OS = D, dropout C ~ Exp(c)): the expectation is computed by conditioning on the death time D. Given
the censoring times tau_P = min(C, t_p - e) and tau_O = min(C, t_o - e),
  E[h_P | D = y] = M1(y) + exp(-a y) r1_P(y) for y <= tau_P and K_P = M1(tau_P) + exp(-a tau_P) r0_P(tau_P) else,
  M1(y) = int_0^y a exp(-a x) r1_P(x) dx,
and the expectation over D, C and the entry time is computed with ODE-based cumulative functions and nested
adaptive quadrature (scipy quad_vec). The R code instead splits the product of counting-process integrals into
four terms (covariance note, equations (6)-(9)) and integrates them over the PFS and OS times.

The standardized covariance combines these moments with the derivatives of the per-patient score mean,
information and expected events with respect to calendar time, obtained here by central differences of
ref_logrank.moments() and ref_logrank.events() (the R code integrates against the accrual density).
"""
import numpy as np
from scipy import integrate

import ref_logrank as rl
import ref_model as rm

KEYS = ["psi", "phi", "dlt"]
ODE = dict(method="DOP853", rtol=1e-12, atol=1e-15)
QV = dict(epsabs=1e-13, epsrel=1e-11, limit=400)


def _rho_at(alpha, rho, e):
    j = int(np.searchsorted(alpha, e, side="right")) - 1
    return rho[min(max(j, 0), rho.size - 1)]


# ---------------------------------------------------------------------------
# Same endpoint
# ---------------------------------------------------------------------------
class SameFuncs:
    """ODE-based cumulative functions of one endpoint for both groups."""

    def __init__(self, scn, endpoint, tmax):
        self.scn = scn
        self.endpoint = endpoint
        br = [0.0] + [x for x in scn.breaks() if 0 < x < tmax] + [tmax]
        self.sols = []
        y0 = np.zeros(2 * 14)
        for lo, hi in zip(br[:-1], br[1:]):
            sol = integrate.solve_ivp(self._rhs, (lo, hi), y0, dense_output=True, **ODE)
            self.sols.append((lo, hi, sol))
            y0 = sol.y[:, -1]

    def _groups(self, s):
        out = []
        for m, p in ((self.scn.control, self.scn.p0), (self.scn.treatment, self.scn.p1)):
            S = rm.surv(m, s, self.endpoint)
            f = rm.dens(m, s, self.endpoint)
            G = rm.dropout_surv(m, s)
            out.append((p, S * G, f * G, S * G * rm.dropout_haz(m, s)))
        return out

    def _rhs(self, s, y):
        (p0, R0, e0, c0), (p1, R1, e1, c1) = self._groups(s)
        r = p0 * R0 + p1 * R1
        q = p1 * R1 / r
        lbar = (p0 * e0 + p1 * e1) / r
        dy = np.zeros_like(y)
        for z, (e, c) in enumerate(((e0, c0), (e1, c1))):
            off = 14 * z
            g = 1 - q if z == 1 else -q
            L, K = y[off], y[off + 1]
            dy[off] = g * lbar
            dy[off + 1] = (1 - 2 * q) * g * lbar
            B = np.array([-L, K, 0.0])
            r1 = np.array([g, q * (1 - q), 1.0]) + B
            r0 = B
            for a in range(3):
                for b in range(3):
                    dy[off + 2 + 3 * a + b] = r1[a] * r1[b] * e + r0[a] * r0[b] * c
                dy[off + 11 + a] = r1[a] * e + r0[a] * c
        return dy

    def y(self, x):
        for lo, hi, sol in self.sols:
            if lo - 1e-12 <= x <= hi + 1e-12:
                return sol.sol(min(max(x, lo), hi))
        raise ValueError(x)

    def parts(self, x, z):
        yy = self.y(x)
        off = 14 * z
        B = np.array([-yy[off], yy[off + 1], 0.0])
        A = yy[off + 2:off + 11].reshape(3, 3)
        M = yy[off + 11:off + 14]
        (p0, R0, _, _), (p1, R1, _, _) = self._groups(x)
        return B, A, M, (R0, R1)[z]


def same_moments(scn, n, times, endpoint):
    """Per-patient means E[r_a(t_i)] and products E[r_a(t_i) r_b(t_j)] (both groups)."""
    alpha, rho, Aend = scn.accrual.setup(n)
    times = np.asarray(times, float)
    k = times.size
    sf = SameFuncs(scn, endpoint, float(times.max()))
    mean = np.zeros((k, 3))
    prod = np.zeros((k, k, 3, 3))
    pz = (scn.p0, scn.p1)
    for i in range(k):
        for j in range(k):
            t, tp = times[i], times[j]
            tl, th = min(t, tp), max(t, tp)
            hi = min(tl, Aend)
            pts = [x for x in alpha if 0 < x < hi] + [tl - b for b in scn.breaks() if 0 < tl - b < hi]
            pts += [th - b for b in scn.breaks() if 0 < th - b < hi]

            def integrand(e):
                W, Wp = tl - e, th - e
                out = np.zeros(9)
                for z in (0, 1):
                    B, A, M, _ = sf.parts(W, z)
                    Bp, _, Mp, Rp = sf.parts(Wp, z)
                    val = A + np.outer(B, Mp - M + Bp * Rp)
                    out += pz[z] * val.ravel()
                return _rho_at(alpha, rho, e) * out
            res = integrate.quad_vec(integrand, 0, hi, points=sorted(set(pts)) or None, **QV)[0].reshape(3, 3)
            # res[a, b] = E[r_a(earlier) r_b(later)]
            prod[i, j] = res if t <= tp else res.T
        t = times[i]
        hi = min(t, Aend)
        pts = [x for x in alpha if 0 < x < hi] + [t - b for b in scn.breaks() if 0 < t - b < hi]

        def mint(e):
            out = np.zeros(3)
            for z in (0, 1):
                B, A, M, R = sf.parts(t - e, z)
                out += pz[z] * (M + B * R)
            return _rho_at(alpha, rho, e) * out
        mean[i] = integrate.quad_vec(mint, 0, hi, points=sorted(set(pts)) or None, **QV)[0]
    return mean, prod


# ---------------------------------------------------------------------------
# PFS and OS, Fleischer model
# ---------------------------------------------------------------------------
class FleischerGroup:
    """Exponential progression a, death b (same before and after progression), dropout c."""

    def __init__(self, a, b, c):
        self.a, self.b, self.c = float(a), float(b), float(c)


class CrossFleischer:
    """Cumulative functions for E[r_P(t_p) r_O(t_o)] in the Fleischer model."""

    def __init__(self, groups, p, tmax):
        self.g = groups
        self.p = p
        self.tmax = tmax
        y0 = np.zeros(2 * 4 + 2 * 21)
        self.sol = integrate.solve_ivp(self._rhs, (0.0, tmax), y0, dense_output=True, **ODE)

    def _qlbar(self, s, ep):
        R, ev = [], []
        for gr, pg in zip(self.g, self.p):
            lam = gr.a + gr.b if ep == "pfs" else gr.b
            Rg = np.exp(-(lam + gr.c) * s)
            R.append(pg * Rg)
            ev.append(pg * lam * Rg)
        r = R[0] + R[1]
        return R[1] / r, (ev[0] + ev[1]) / r

    def _rhs(self, s, y):
        # per group z: L_P, K_P, L_O, K_O (first 8 entries), then from 8 + 21 z: M1_a (3), Q1_ab (9),
        # Q2_b (3), Q3_a (3) and 3 unused entries
        dy = np.zeros_like(y)
        qP, lP = self._qlbar(s, "pfs")
        qO, lO = self._qlbar(s, "os")
        for z, gr in enumerate(self.g):
            gP = 1 - qP if z == 1 else -qP
            gO = 1 - qO if z == 1 else -qO
            dy[4 * z] = gP * lP
            dy[4 * z + 1] = (1 - 2 * qP) * gP * lP
            dy[4 * z + 2] = gO * lO
            dy[4 * z + 3] = (1 - 2 * qO) * gO * lO
            r1P, _ = self._r(s, y, z, "pfs")
            r1O, _ = self._r(s, y, z, "os")
            off = 8 + 21 * z
            M1 = y[off:off + 3]
            inner = M1 + np.exp(-gr.a * s) * r1P            # E[h_P | D = s] for s <= tau_P
            dens_d = gr.b * np.exp(-gr.b * s)
            dy[off:off + 3] = gr.a * np.exp(-gr.a * s) * r1P
            dy[off + 3:off + 12] = (dens_d * np.outer(inner, r1O)).ravel()
            dy[off + 12:off + 15] = dens_d * r1O
            dy[off + 15:off + 18] = dens_d * inner
        return dy

    def _r(self, s, y, z, ep):
        qP, _ = self._qlbar(s, ep)
        g = 1 - qP if z == 1 else -qP
        L, K = (y[4 * z], y[4 * z + 1]) if ep == "pfs" else (y[4 * z + 2], y[4 * z + 3])
        B = np.array([-L, K, 0.0])
        return np.array([g, qP * (1 - qP), 1.0]) + B, B

    def F(self, z, tauP, tauO):
        """E[r_P r_O | tau_P, tau_O] for a patient of group z (3 x 3)."""
        gr = self.g[z]
        off = 8 + 21 * z
        yP = self.sol.sol(tauP)
        yO = self.sol.sol(tauO)
        _, r0P = self._r(tauP, yP, z, "pfs")
        _, r0O = self._r(tauO, yO, z, "os")
        M1P = yP[off:off + 3]
        KP = M1P + np.exp(-gr.a * tauP) * r0P
        Q1 = lambda yy: yy[off + 3:off + 12].reshape(3, 3)
        Q2 = lambda yy: yy[off + 12:off + 15]
        Q3 = lambda yy: yy[off + 15:off + 18]
        if tauO <= tauP:
            return Q1(yO) + np.outer(Q3(yP) - Q3(yO) + KP * np.exp(-gr.b * tauP), r0O)
        return Q1(yP) + np.outer(KP, Q2(yO) - Q2(yP)) + np.outer(KP, r0O) * np.exp(-gr.b * tauO)


def cross_moments_fleischer(groups, p, alpha, rho, t_p, t_o):
    """Per-patient E[r_P(t_p) r_O(t_o)] (rows PFS, columns OS) for piecewise uniform accrual."""
    tmax = max(t_p, t_o)
    cf = CrossFleischer(groups, p, tmax)
    out = np.zeros((3, 3))
    A_end = alpha[-1]
    for z, gr in enumerate(groups):
        def over_c(e):
            WP, WO = max(t_p - e, 0.0), max(t_o - e, 0.0)
            Wmax = max(WP, WO)
            tail = np.exp(-gr.c * Wmax) * cf.F(z, WP, WO).ravel()
            if Wmax <= 0:
                return tail
            f = lambda g: gr.c * np.exp(-gr.c * g) * cf.F(z, min(g, WP), min(g, WO)).ravel()
            pts = [x for x in (WP, WO) if 0 < x < Wmax]
            body = integrate.quad_vec(f, 0.0, Wmax, points=pts or None, **QV)[0]
            return body + tail
        hi = min(A_end, max(t_p, t_o))
        pts = [x for x in list(alpha) + [t_p, t_o] if 0 < x < hi]
        val = integrate.quad_vec(lambda e: _rho_at(alpha, rho, e) * over_c(e), 0.0, hi,
                                 points=sorted(set(pts)) or None, **QV)[0]
        out += p[z] * val.reshape(3, 3)
    return out


# ---------------------------------------------------------------------------
# Standardized covariance
# ---------------------------------------------------------------------------
def per_patient_core(scn, n, t, ep, h=1e-4):
    """mu, v, their time derivatives, d and its derivative (central differences of the reference moments)."""
    def mvd(x):
        m = rl.moments(scn, n, x, ep)
        return np.array([m["score_mean"] / n, m["info"] / n, m["events"] / n])
    c = mvd(t)
    der = (mvd(t + h) - mvd(t - h)) / (2 * h)
    mu, v, d = c
    mdot = -der[0] / np.sqrt(v) + mu * der[1] / (2 * v ** 1.5)
    return dict(mu=mu, v=v, d=d, mdot=mdot, ddot=der[2])


def zeta_cov(items, mom_same, mom_cross, core):
    """Covariance of zeta for items (endpoint, time, trigger).

    mom_same[ep] = (times, mean, prod); mom_cross[(t_p, t_o)] = E[r_P r_O] (3 x 3);
    core[(ep, t)] = per_patient_core(...).
    """
    def cov_el(x, y):
        (e1, a, t1), (e2, b, t2) = x, y
        if e1 == e2:
            times, mean, prod = mom_same[e1]
            i, j = list(times).index(t1), list(times).index(t2)
            return prod[i, j, a, b] - mean[i, a] * mean[j, b]
        if e1 == "os":
            return cov_el(y, x)
        tP, mP, _ = mom_same["pfs"]
        tO, mO, _ = mom_same["os"]
        i, j = list(tP).index(t1), list(tO).index(t2)
        return mom_cross[(t1, t2)][a, b] - mP[i, a] * mO[j, b]

    def comps(item):
        ep, t, trig = item
        c = core[(ep, t)]
        out = [((ep, 0, t), -1 / np.sqrt(c["v"])), ((ep, 1, t), c["mu"] / (2 * c["v"] ** 1.5))]
        if trig is not None:
            out.append(((trig, 2, t), -c["mdot"] / core[(trig, t)]["ddot"]))
        return out
    k = len(items)
    S = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            S[i, j] = sum(ci * cj * cov_el(xi, xj) for xi, ci in comps(items[i]) for xj, cj in comps(items[j]))
    return S
