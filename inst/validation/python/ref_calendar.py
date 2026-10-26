"""Reference computations for switching that starts at a calendar time (switch_start of ms_scenario(), stage 3).

Written independently of the R code in DesignComplexSurvival. A patient enrolled at e can switch from the
follow-up time o(e) = max(switch_start - e, 0) on (Planning/covariance_note/calendar_switch_note.pdf in the
project folder). The R code uses closed forms of the post-progression transitions, splits the entry range at
switch_start - s and switch_start, swaps the integrals over the entry and the follow-up time in the products of the
per-patient functionals, and differentiates the moments analytically. This module uses
  - for the OS given o: adaptive quadrature (scipy.integrate.quad) over the progression time with the
    post-progression quantities of ref_model5.post() (adaptive quadrature over the switch time) [os_quad], and,
    for Markov models (clock "forward", switch clock "same"), products of matrix exponentials with the
    generator without switching before o and with switching from o on [os_markov];
  - for the group functions, expected events, score mean and information: Gauss-Legendre rules over the whole
    entry range (no split) with os_markov;
  - for the covariance of the per-patient functionals: equation (5) of the covariance note given the entry time,
    integrated over the entry time, with the functionals of each look represented by interpolating polynomials
    on the pieces of its rule;
  - for the time derivatives (looks driven by OS events): central differences.
"""
import numpy as np
from scipy import integrate, linalg

import ref_model5 as R5

QUAD = dict(epsabs=1e-12, epsrel=1e-10, limit=400)


# OS given the availability o ----------------------------------------------------------------------------------
def os_quad(m, s, o):
    """OS survival and density at s for switching possible from o on (scalar; any clock and switch model of
    ref_model5.Model5 without switching before progression)."""
    tm, tp = m.time, m.time_post
    fwd = m.clock == "forward"

    def part(u, which):
        t_u = u if fwd else 0.0
        t_s = s if fwd else s - u
        if u >= o:
            q = R5.post(m, t_u, t_s)
            p = R5.val(u, tm, m.switch_prob)
            if which == 0:
                return (1 - p) * (q["r22"] + q["q23"]) + p * q["s33"]
            return (1 - p) * (q["f22"] + q["f23"]) + p * q["f33"]
        t_a = min(o if fwd else o - u, t_s)
        pre = np.exp(-(R5.cum(t_a, tp, m.death_prog) - R5.cum(t_u, tp, m.death_prog)))
        q = R5.post(m, t_a, t_s)
        return pre * ((q["r22"] + q["q23"]) if which == 0 else (q["f22"] + q["f23"]))

    def k1(u):
        return R5.p1(m, u) * R5.val(u, tm, m.prog)
    cands = list(tm) + [o] + list(tp) + [s - x for x in tp] + [o - x for x in tp]
    if m.mode() == "switch":
        cands += [s - x - y for x in tp for y in tp]
    if m.mode() == "aft" and m.switch_af != 1:
        phi = m.switch_af
        cands += [s - (x - (1 - phi) * y) / phi for x in tp for y in tp]
        cands += [o * (1 - phi) + phi * s - x for x in tp]
        cands += [(x - phi * s) / (1 - phi) for x in tp]
    pts = R5._pts(cands, 0, s)
    p11 = R5.p1(m, s)
    if s <= 0 or np.all(m.prog == 0):
        return p11, R5.val(s, tm, m.death) * p11
    sv = integrate.quad(lambda u: k1(u) * part(u, 0), 0, s, points=pts, **QUAD)[0]
    dv = integrate.quad(lambda u: k1(u) * part(u, 1), 0, s, points=pts, **QUAD)[0]
    return p11 + sv, R5.val(s, tm, m.death) * p11 + dv


def _gen(m, x, switching):
    """Generators (N x 4 x 4) of states (1, 2, 3, 4) at the times x (forward clock)."""
    tm, tp = m.time, m.time_post
    idx = lambda brk: np.maximum(np.searchsorted(brk, x, side="right") - 1, 0)
    a = m.prog[idx(tm)]; b = m.death[idx(tm)]
    p = m.switch_prob[idx(tm)] if switching else np.zeros(x.size)
    r = m.switch_rate[idx(tp)] if switching else np.zeros(x.size)
    c = m.death_prog[idx(tp)]; d = m.death_switch[idx(tp)]
    q = np.zeros((x.size, 4, 4))
    q[:, 0, 0] = -(a + b); q[:, 0, 1] = a * (1 - p); q[:, 0, 2] = a * p; q[:, 0, 3] = b
    q[:, 1, 1] = -(r + c); q[:, 1, 2] = r; q[:, 1, 3] = c
    q[:, 2, 2] = -d; q[:, 2, 3] = d
    return q


def _transfer(m, x, switching):
    """P(0, x) (N x 4 x 4) for the piecewise constant generator."""
    brk = np.unique(np.concatenate([[0.0], m.time, m.time_post]))
    fb = [np.eye(4)]
    for lo, hi in zip(brk[:-1], brk[1:]):
        fb.append(fb[-1] @ linalg.expm(_gen(m, np.array([lo]), switching)[0] * (hi - lo)))
    fb = np.array(fb)
    x = np.asarray(x, float)
    k = np.maximum(np.searchsorted(brk, x, side="right") - 1, 0)
    q = _gen(m, brk[k], switching)
    return np.einsum("nij,njk->nik", fb[k], linalg.expm(q * (x - brk[k])[:, None, None]))


def os_markov(m, s, o):
    """OS survival and density at s for switching possible from o on, Markov models (vectorized pairs)."""
    s = np.atleast_1d(np.asarray(s, float)); o = np.broadcast_to(np.asarray(o, float), s.shape)
    mm = np.minimum(o, s)
    row = _transfer(m, mm, False)[:, 0, :]
    fs_m = _transfer(m, mm, True)
    fs_s = _transfer(m, s, True)
    prob = np.einsum("ni,nij->nj", row, np.linalg.solve(fs_m, fs_s))
    tm, tp = m.time, m.time_post
    b = m.death[np.maximum(np.searchsorted(tm, s, side="right") - 1, 0)]
    c = m.death_prog[np.maximum(np.searchsorted(tp, s, side="right") - 1, 0)]
    d = m.death_switch[np.maximum(np.searchsorted(tp, s, side="right") - 1, 0)]
    return prob[:, :3].sum(axis=1), prob[:, 0] * b + prob[:, 1] * c + prob[:, 2] * d


def has_switch(m):
    return bool(np.any(m.prog > 0) and (np.any(m.switch_prob > 0) or np.any(m.switch_rate > 0)))


# Gauss-Legendre rules and piecewise polynomials -----------------------------------------------------------------
def rule(a, b, breaks, h, n):
    if not b > a:
        return dict(x=np.zeros(0), w=np.zeros(0), lo=np.zeros(0), hi=np.zeros(0))
    pts = np.unique(np.concatenate([[a], [x for x in breaks if a < x < b], [b]]))
    sub = [np.linspace(lo, hi, int(np.ceil((hi - lo) / h - 1e-9)) + 1)[:-1] for lo, hi in zip(pts[:-1], pts[1:])]
    pts = np.concatenate(sub + [[b]])
    xr, wr = np.polynomial.legendre.leggauss(n)
    lo, hi = pts[:-1], pts[1:]
    half = (hi - lo) / 2; mid = (hi + lo) / 2
    return dict(x=(mid[:, None] + half[:, None] * xr[None, :]).ravel(), w=(half[:, None] * wr[None, :]).ravel(),
                lo=lo, hi=hi, n=n)


class Poly:
    """Interpolating polynomials on the pieces of a rule (Legendre basis): value and integral from the start."""

    def __init__(self, g, vals):
        vals = np.asarray(vals, float)
        self.vals2 = vals.ndim == 2
        vals = vals if self.vals2 else vals[:, None]
        n = g["n"]; npc = g["lo"].size
        self.lo, self.hi, self.n = g["lo"], g["hi"], n
        xr, _ = np.polynomial.legendre.leggauss(n)
        vinv = np.linalg.inv(np.polynomial.legendre.legvander(xr, n - 1))
        self.c = np.einsum("ij,pjf->pif", vinv, vals.reshape(npc, n, -1))
        tot = (self.hi - self.lo)[:, None] * self.c[:, 0, :]
        self.cum = np.vstack([np.zeros((1, tot.shape[1])), np.cumsum(tot, axis=0)])

    def _loc(self, x):
        x = np.clip(np.asarray(x, float), self.lo[0], self.hi[-1])
        k = np.clip(np.searchsorted(self.lo, x, side="right") - 1, 0, self.lo.size - 1)
        xi = np.clip((2 * x - self.lo[k] - self.hi[k]) / (self.hi[k] - self.lo[k]), -1, 1)
        return k, xi

    def value(self, x):
        k, xi = self._loc(x)
        out = np.einsum("ni,nif->nf", np.polynomial.legendre.legvander(xi, self.n - 1), self.c[k])
        return out if self.vals2 else out[:, 0]

    def integral(self, x):
        k, xi = self._loc(x)
        out = np.zeros((xi.size, self.c.shape[2]))
        for j in range(self.n):
            e = np.zeros(j + 1); e[j] = 1
            ij = np.polynomial.legendre.legint(e, lbnd=-1)
            out += np.polynomial.legendre.legval(xi, ij)[:, None] * self.c[k, j, :]
        out = self.cum[k] + (self.hi[k] - self.lo[k])[:, None] / 2 * out
        return out if self.vals2 else out[:, 0]


# Calendar-time functions (Markov models) -------------------------------------------------------------------------
class Design:
    """Two groups (Model5, forward clock), allocation probability p1 of the treatment group, accrual given by
    interval starts, proportions and duration, calendar time tc of the start of switching."""

    def __init__(self, m0, m1, p1, acc_time, acc_prop, duration, tc, h=0.5, n=10):
        self.m = (m0, m1); self.p = (1 - p1, p1); self.tc = float(tc)
        self.alpha = np.concatenate([acc_time, [duration]]); self.dur = duration
        self.rho = np.asarray(acc_prop, float) / np.diff(self.alpha)
        self.h = h; self.n = n

    def dens_e(self, e):
        k = np.clip(np.searchsorted(self.alpha, e, side="right") - 1, 0, self.rho.size - 1)
        return np.where((e >= 0) & (e < self.dur), self.rho[k], 0.0)

    def os(self, j, s, o):
        m = self.m[j]
        if has_switch(m):
            return os_markov(m, s, o)
        return os_markov(m, s, np.zeros(np.size(s)))

    def G(self, j, s):
        m = self.m[j]
        return np.array([np.exp(-R5.cum(x, m.time, m.dropout)) for x in np.atleast_1d(s)])

    def kap(self, j, s):
        m = self.m[j]
        return np.array([R5.val(x, m.time, m.dropout) for x in np.atleast_1d(s)])

    def breaks_s(self, t):
        b = []
        for m in self.m:
            b += list(m.time) + list(m.time_post)
        return np.array(b + list(t - self.alpha) + [t - self.tc, self.tc] + list(self.tc - self.alpha))

    def group(self, j, t, s):
        """Y and D (without dropout and allocation) at follow-up times s, by rules over the whole entry range."""
        m = self.m[j]
        ys = np.zeros(s.size); ds = np.zeros(s.size)
        for i, si in enumerate(s):
            hi = min(t - si, self.dur)
            if hi <= 0:
                continue
            br = list(self.alpha) + [self.tc, self.tc - si] + list(self.tc - m.time) + list(self.tc - m.time_post) + \
                list(self.tc - si + m.time_post)
            g = rule(0, hi, br, self.h, self.n)
            S, f = self.os(j, np.full(g["x"].size, si), np.maximum(self.tc - g["x"], 0))
            w = g["w"] * self.dens_e(g["x"])
            ys[i] = np.sum(w * S); ds[i] = np.sum(w * f)
        return ys, ds

    def look(self, t):
        g = rule(0, t, self.breaks_s(t), self.h, self.n)
        s = g["x"]
        out = dict(t=t, g=g, Y=[], D=[], G=[], kap=[])
        for j in range(2):
            Y, D = self.group(j, t, s)
            out["Y"].append(Y); out["D"].append(D); out["G"].append(self.G(j, s)); out["kap"].append(self.kap(j, s))
        y = [self.p[j] * out["G"][j] * out["Y"][j] for j in range(2)]
        d = [self.p[j] * out["G"][j] * out["D"][j] for j in range(2)]
        ysum = y[0] + y[1]; ok = ysum > 0
        q1 = np.full(s.size, self.p[1]); lb = np.zeros(s.size); u = np.zeros(s.size)
        q1[ok] = y[1][ok] / ysum[ok]; lb[ok] = (d[0][ok] + d[1][ok]) / ysum[ok]
        u[ok] = (y[0][ok] * d[1][ok] - y[1][ok] * d[0][ok]) / ysum[ok]
        out.update(q1=q1, lbar=lb, mu=np.sum(g["w"] * u), v=np.sum(g["w"] * q1 * (1 - q1) * (d[0] + d[1])),
                   ev=[np.sum(g["w"] * d[j]) for j in range(2)],
                   do=[np.sum(g["w"] * y[j] * out["kap"][j]) for j in range(2)])
        funs = []
        for j in range(2):
            gj = j - q1
            rho = np.column_stack([gj, q1 * (1 - q1), np.ones(s.size)])
            beta = np.column_stack([-gj * lb, (1 - 2 * q1) * gj * lb, np.zeros(s.size)])
            funs.append(dict(rho=rho, beta=beta, pr=Poly(g, rho), pb=Poly(g, beta)))
        out["funs"] = funs
        return out

    def mean_funs(self, L):
        avg = np.zeros(3)
        for j in range(2):
            f = L["funs"][j]; wt = L["g"]["w"] * L["G"][j]
            avg += self.p[j] * ((f["rho"] * (wt * L["D"][j])[:, None]).sum(0) +
                                (f["beta"] * (wt * L["Y"][j])[:, None]).sum(0))
        return avg

    def products(self, L1, L2):
        """Sum over the groups of p_j E_j[r(t1) r'(t2)] by equation (5) of the covariance note given the entry."""
        t1, t2 = L1["t"], L2["t"]
        out = np.zeros((3, 3))
        g_e = rule(0, min(t1, self.dur), list(self.alpha) + [self.tc] + list(t1 - self.breaks_s(t1)) +
                   list(t2 - self.breaks_s(t2)), self.h, self.n)
        for j in range(2):
            m = self.m[j]; fa = L1["funs"][j]; fb = L2["funs"][j]
            for e, we in zip(g_e["x"], g_e["w"] * self.dens_e(g_e["x"])):
                W1, W2 = t1 - e, t2 - e
                o = max(self.tc - e, 0)
                br = list(m.time) + list(m.time_post) + [self.tc - e, W1] + list(t1 - self.alpha) + \
                    list(t2 - self.alpha) + [t1 - self.tc, t2 - self.tc] + list(W1 + m.time_post) + \
                    list(o + m.time_post) + list(L1["g"]["lo"]) + list(L2["g"]["lo"])
                g = rule(0, W2, br, self.h, self.n)
                s = g["x"]; w = g["w"]
                S, f = self.os(j, s, np.full(s.size, o))
                Gs = self.G(j, s); ev = f * Gs; cv = S * Gs * self.kap(j, s)
                inw = s < W1
                B1 = fa["pb"].integral(np.minimum(s, W1)); r1a = fa["pr"].value(np.minimum(s, W1)) + B1
                B2 = fb["pb"].integral(s); r1b = fb["pr"].value(s) + B2
                part_in = ((r1a * (w * ev * inw)[:, None]).T @ r1b + (B1 * (w * cv * inw)[:, None]).T @ B2)
                part_out = (w * ev * ~inw) @ r1b + (w * cv * ~inw) @ B2
                SW2, _ = self.os(j, np.array([W2]), np.array([o]))
                piW2 = SW2[0] * self.G(j, np.array([W2]))[0]
                BW1 = fa["pb"].integral(np.array([W1]))[0]; BW2 = fb["pb"].integral(np.array([W2]))[0]
                out += self.p[j] * we * (part_in + np.outer(BW1, part_out + BW2 * piW2))
        return out

    def zeta(self, looks, trigger=None, h_diff=0.01):
        """Covariance of zeta (studentized), of psi divided by sqrt(v v') ("alternative") and the correlation."""
        k = len(looks)
        L = [self.look(t) for t in looks]
        mu = np.array([x["mu"] for x in L]); v = np.array([x["v"] for x in L])
        avg = [self.mean_funs(x) for x in L]
        cv = np.zeros((3 * k, 3 * k))
        for i in range(k):
            for l in range(i, k):
                blk = self.products(L[i], L[l]) - np.outer(avg[i], avg[l])
                cv[3 * i:3 * i + 3, 3 * l:3 * l + 3] = blk
                cv[3 * l:3 * l + 3, 3 * i:3 * i + 3] = blk.T
        coef = np.zeros((k, 3 * k))
        for i in range(k):
            coef[i, 3 * i] = -1 / np.sqrt(v[i]); coef[i, 3 * i + 1] = mu[i] / (2 * v[i] ** 1.5)
            if trigger is not None and trigger[i] == "os":
                lp, lm = self.look(looks[i] + h_diff), self.look(looks[i] - h_diff)
                m_of = lambda x: -x["mu"] / np.sqrt(x["v"])
                mdot = (m_of(lp) - m_of(lm)) / (2 * h_diff)
                ddot = (sum(lp["ev"]) - sum(lm["ev"])) / (2 * h_diff)
                coef[i, 3 * i + 2] = -mdot / ddot
        stud = coef @ cv @ coef.T
        stud = (stud + stud.T) / 2
        ip = 3 * np.arange(k)
        sd = np.sqrt(np.diag(stud))
        return dict(mu=mu, v=v, studentized=stud, alternative=cv[np.ix_(ip, ip)] / np.sqrt(np.outer(v, v)),
                    null=stud / np.outer(sd, sd), events=np.array([sum(x["ev"]) for x in L]),
                    looks=L)


def delta_null(des, t1, t2):
    """Null hypothesis (identical groups, equal dropout): delta = E[A(t1) A(t2)] - E[Delta(t1)] by equation (3) of
    the calendar note, with A(t) = Delta(t) - Lambda-bar_t(X(t)); returns (delta, E[Delta(t1)], E[Delta(t2)])."""
    L1, L2 = des.look(t1), des.look(t2)
    pl1, pl2 = Poly(L1["g"], L1["lbar"]), Poly(L2["g"], L2["lbar"])
    g = L1["g"]
    term1 = -np.sum(g["w"] * (L1["lbar"] - pl2.value(g["x"])) * pl1.integral(g["x"]) * L1["G"][0] * L1["Y"][0])
    g_e = rule(0, min(t1, des.dur), list(des.alpha) + [des.tc] + list(t1 - des.breaks_s(t1)), des.h, des.n)
    term2 = 0.0
    m = des.m[0]
    for e, we in zip(g_e["x"], g_e["w"] * des.dens_e(g_e["x"])):
        o = max(des.tc - e, 0)
        gu = rule(t1 - e, t2 - e, list(m.time) + list(m.time_post) + [des.tc - e] + list(t2 - des.alpha) +
                  [t2 - des.tc] + list(o + m.time_post) + list(L2["g"]["lo"]), des.h, des.n)
        S, f = des.os(0, gu["x"], np.full(gu["x"].size, o))
        term2 += we * pl1.integral(np.array([t1 - e]))[0] * np.sum(gu["w"] * (f - pl2.value(gu["x"]) * S) *
                                                                   des.G(0, gu["x"]))
    return term1 - term2, sum(L1["ev"]), sum(L2["ev"])
