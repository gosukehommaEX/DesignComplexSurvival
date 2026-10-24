"""Reference implementation of the multistate model with the options of stage 3 (one group).

Written independently of the R code in DesignComplexSurvival (stage 3, candidates 2, 4 and 5). The R code uses
closed forms on the pieces of the intensities and Gauss-Legendre rules; this module uses scalar loops, adaptive
quadrature (scipy.integrate.quad) over the switch, progression and entry times, products of matrix exponentials
(scipy.linalg.expm) for Markov models, and a patient-level simulation.

States: 1 progression-free (no switch), 5 progression-free after a switch before progression, 2 progressed
without switch, 3 switched, 4 dead.
  1 -> 5 switch_pre, 1 -> 2 prog (1 - switch_prob), 1 -> 3 prog switch_prob, 1 -> 4 death,
  5 -> 3 prog_switch, 5 -> 4 death_switch_pre (time since randomization),
  2 -> 3 switch_rate, 2 -> 4 death_prog (post-progression clock),
  3 -> 4: after a switch at the post-progression clock value sigma, at the clock value tau,
      switch_clock "same":   d(tau)                            (d = death_switch)
      switch_clock "switch": d(tau - sigma)
      switch_af = phi:       phi c(sigma + phi (tau - sigma))  (c = death_prog)
    after progression from state 5: d(tau) (phi c(tau) with switch_af).
The post-progression clock is the time since progression ("reset") or since randomization ("forward").
PFS: with pfs_switch "event" the time to leave {1, 5} by progression or death; with "censor" PFS is censored at
the switch 1 -> 5. dropout censors both endpoints, dropout_pfs censors PFS only.
"""
import numpy as np
from scipy import integrate, linalg

QUAD = dict(epsabs=1e-13, epsrel=1e-11, limit=400)
QIN = dict(epsabs=1e-13, epsrel=1e-11, limit=200)


def _expand(x, n):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    if x.size == 1:
        return np.repeat(x, n)
    if x.size != n:
        raise ValueError("length mismatch")
    return x


class Model5:
    """One group (arguments as in ms_model(), plus pfs_switch of ms_scenario())."""

    def __init__(self, death, time=(0.0,), prog=0.0, death_prog=0.0, switch_prob=0.0, switch_rate=0.0,
                 death_switch=None, dropout=0.0, time_post=(0.0,), clock="reset", switch_clock="same",
                 switch_af=None, switch_pre=0.0, prog_switch=None, death_switch_pre=None, dropout_pfs=0.0,
                 pfs_switch="event"):
        self.time = np.asarray(time, dtype=float)
        n = self.time.size
        self.prog = _expand(prog, n)
        self.death = _expand(death, n)
        self.switch_prob = _expand(switch_prob, n)
        self.dropout = _expand(dropout, n)
        self.switch_pre = _expand(switch_pre, n)
        self.prog_switch = self.prog.copy() if prog_switch is None else _expand(prog_switch, n)
        self.death_switch_pre = self.death.copy() if death_switch_pre is None else _expand(death_switch_pre, n)
        self.dropout_pfs = _expand(dropout_pfs, n)
        self.time_post = np.asarray(time_post, dtype=float)
        m = self.time_post.size
        self.death_prog = _expand(death_prog, m)
        self.switch_rate = _expand(switch_rate, m)
        self.switch_af = None if switch_af is None else float(switch_af)
        if self.switch_af is not None:
            self.death_switch = self.switch_af * self.death_prog
        elif death_switch is None:
            self.death_switch = self.death_prog.copy()
        else:
            self.death_switch = _expand(death_switch, m)
        self.clock = clock
        self.switch_clock = switch_clock
        self.pfs_switch = pfs_switch

    def mode(self):
        if self.switch_af is not None:
            return "aft"
        return self.switch_clock


# step functions (scalar loops) ----------------------------------------------------------------------------------
def val(x, brk, v):
    """Value at x of the right-continuous step function with values v on intervals starting at brk."""
    j = 0
    for i in range(brk.size):
        if x >= brk[i]:
            j = i
    return float(v[j])


def cum(x, brk, v):
    """Integral from 0 to max(x, 0) of the step function."""
    x = max(float(x), 0.0)
    tot = 0.0
    for i in range(brk.size):
        lo = brk[i]
        hi = brk[i + 1] if i + 1 < brk.size else np.inf
        if x <= lo:
            break
        tot += v[i] * (min(x, hi) - lo)
    return tot


def _pts(cands, lo, hi):
    return sorted(set(float(x) for x in cands if lo < x < hi)) or None


# after a switch ----------------------------------------------------------------------------------------------
def s3(m, sig, tau):
    """Survival from the switch at clock sigma to clock tau."""
    tp = m.time_post
    mode = m.mode()
    if mode == "same":
        return np.exp(-(cum(tau, tp, m.death_switch) - cum(sig, tp, m.death_switch)))
    if mode == "switch":
        return np.exp(-cum(tau - sig, tp, m.death_switch))
    phi = m.switch_af
    return np.exp(-(cum(sig + phi * (tau - sig), tp, m.death_prog) - cum(sig, tp, m.death_prog)))


def h3(m, sig, tau):
    """Death hazard at clock tau after a switch at clock sigma."""
    tp = m.time_post
    mode = m.mode()
    if mode == "same":
        return val(tau, tp, m.death_switch)
    if mode == "switch":
        return val(tau - sig, tp, m.death_switch)
    phi = m.switch_af
    return phi * val(sig + phi * (tau - sig), tp, m.death_prog)


def post(m, t0, t):
    """Post-progression quantities from progression at clock t0 to clock t, by adaptive quadrature over the time
    of a switch after progression: r22 (in state 2), q23 (switched later, alive), f22 and f23 (death densities),
    s33 and f33 (switch at progression)."""
    tp = m.time_post
    y = m.switch_rate + m.death_prog
    r22 = np.exp(-(cum(t, tp, y) - cum(t0, tp, y)))
    out = dict(r22=r22, f22=r22 * val(t, tp, m.death_prog), s33=s3(m, t0, t), f33=s3(m, t0, t) * h3(m, t0, t))
    if np.all(m.switch_rate == 0) or not t > t0:
        out.update(q23=0.0, f23=0.0)
        return out
    g = lambda x: np.exp(-(cum(x, tp, y) - cum(t0, tp, y))) * val(x, tp, m.switch_rate) * s3(m, x, t)
    cands = list(tp) + [t - x for x in tp]
    if m.mode() == "aft" and m.switch_af != 1:
        cands += [(x - m.switch_af * t) / (1 - m.switch_af) for x in tp]
    pts = _pts(cands, t0, t)
    out["q23"] = integrate.quad(g, t0, t, points=pts, **QIN)[0]
    out["f23"] = integrate.quad(lambda x: g(x) * h3(m, x, t), t0, t, points=pts, **QIN)[0]
    return out


def clock_pair(m, v, s):
    """Post-progression clock values at progression (time v) and at time s since randomization."""
    return (v, s) if m.clock == "forward" else (0.0, s - v)


def post5(m, v, s):
    """Survival and density at s after progression from state 5 at v (death_switch on the post clock)."""
    a, b = clock_pair(m, v, s)
    tp = m.time_post
    sv = np.exp(-(cum(b, tp, m.death_switch) - cum(a, tp, m.death_switch)))
    return sv, sv * val(b, tp, m.death_switch)


# before progression ----------------------------------------------------------------------------------------
def p1(m, s):
    return np.exp(-cum(s, m.time, m.prog + m.death + m.switch_pre))


def p5(m, s):
    if s <= 0 or np.all(m.switch_pre == 0):
        return 0.0
    h5 = m.prog_switch + m.death_switch_pre
    f = lambda w: p1(m, w) * val(w, m.time, m.switch_pre) * np.exp(-(cum(s, m.time, h5) - cum(w, m.time, h5)))
    return integrate.quad(f, 0, s, points=_pts(m.time, 0, s), **QUAD)[0]


def _conv_pts(m, s):
    tp = m.time_post
    c = list(m.time)
    if m.clock == "forward":
        c += list(tp)
        if m.mode() == "switch":
            c += [s - x for x in tp]
        if m.mode() == "aft" and m.switch_af != 1:
            c += [(x - m.switch_af * s) / (1 - m.switch_af) for x in tp]
    else:
        c += [s - x for x in tp]
        if m.mode() == "switch":
            c += [s - x - y for x in tp for y in tp]
        if m.mode() == "aft" and m.switch_af != 1:
            c += [s - (x - (1 - m.switch_af) * y) / m.switch_af for x in tp for y in tp]
    return _pts(c, 0, s)


def _conv(m, s, f):
    if s <= 0:
        return 0.0
    return integrate.quad(f, 0, s, points=_conv_pts(m, s), **QIN)[0]


def states(m, s):
    """P(state 1), P(state 5), P(state 2), P(state 3), P(state 4) at time s, without dropout."""
    a1, a5 = p1(m, s), p5(m, s)

    def k1(v):
        return p1(m, v) * val(v, m.time, m.prog)

    def q(v):
        a, b = clock_pair(m, v, s)
        return post(m, a, b)
    pp = lambda v: val(v, m.time, m.switch_prob)
    a2 = _conv(m, s, lambda v: k1(v) * (1 - pp(v)) * q(v)["r22"])
    a3 = _conv(m, s, lambda v: k1(v) * ((1 - pp(v)) * q(v)["q23"] + pp(v) * q(v)["s33"]))
    if np.any(m.switch_pre > 0) and np.any(m.prog_switch > 0):
        a3 += _conv(m, s, lambda v: p5(m, v) * val(v, m.time, m.prog_switch) * post5(m, v, s)[0])
    return np.array([a1, a5, a2, a3, 1.0 - a1 - a5 - a2 - a3])


def os_dens(m, s):
    b = val(s, m.time, m.death)
    b5 = val(s, m.time, m.death_switch_pre)
    out = p1(m, s) * b + p5(m, s) * b5

    def g(v):
        a, bb = clock_pair(m, v, s)
        q = post(m, a, bb)
        p = val(v, m.time, m.switch_prob)
        return p1(m, v) * val(v, m.time, m.prog) * ((1 - p) * (q["f22"] + q["f23"]) + p * q["f33"])
    out += _conv(m, s, g)
    if np.any(m.switch_pre > 0) and np.any(m.prog_switch > 0):
        out += _conv(m, s, lambda v: p5(m, v) * val(v, m.time, m.prog_switch) * post5(m, v, s)[1])
    return out


def os_surv(m, s):
    out = p1(m, s) + p5(m, s)

    def g(v):
        a, bb = clock_pair(m, v, s)
        q = post(m, a, bb)
        p = val(v, m.time, m.switch_prob)
        return p1(m, v) * val(v, m.time, m.prog) * ((1 - p) * (q["r22"] + q["q23"]) + p * q["s33"])
    out += _conv(m, s, g)
    if np.any(m.switch_pre > 0) and np.any(m.prog_switch > 0):
        out += _conv(m, s, lambda v: p5(m, v) * val(v, m.time, m.prog_switch) * post5(m, v, s)[0])
    return out


def pfs_surv(m, s):
    """Probability of being at risk for PFS (without censoring)."""
    return p1(m, s) if m.pfs_switch == "censor" else p1(m, s) + p5(m, s)


def pfs_dens(m, s):
    """Density of a PFS event (without censoring)."""
    e1 = p1(m, s) * val(s, m.time, m.prog + m.death)
    if m.pfs_switch == "censor":
        return e1
    return e1 + p5(m, s) * val(s, m.time, m.prog_switch + m.death_switch_pre)


def from5(m, w, u):
    """From state 5 at time w: probability of being alive at u and death density at u."""
    tm = m.time
    h5 = m.prog_switch + m.death_switch_pre
    stay = lambda x: np.exp(-(cum(x, tm, h5) - cum(w, tm, h5)))
    sv = stay(u)
    dv = stay(u) * val(u, tm, m.death_switch_pre)
    if u > w and np.any(m.prog_switch > 0):
        c = list(tm) + (list(m.time_post) if m.clock == "forward" else [u - x for x in m.time_post])
        pts = _pts(c, w, u)
        sv += integrate.quad(lambda v: stay(v) * val(v, tm, m.prog_switch) * post5(m, v, u)[0], w, u,
                             points=pts, **QUAD)[0]
        dv += integrate.quad(lambda v: stay(v) * val(v, tm, m.prog_switch) * post5(m, v, u)[1], w, u,
                             points=pts, **QUAD)[0]
    return sv, dv


def dropout_surv(m, s, endpoint="os"):
    rate = m.dropout + (m.dropout_pfs if endpoint == "pfs" else 0.0)
    return np.exp(-cum(s, m.time, rate))


def censor_haz(m, s, endpoint="os"):
    """Hazard of an observed censoring among the patients at risk (dropout, PFS-only censoring, switch)."""
    h = val(s, m.time, m.dropout)
    if endpoint == "pfs":
        h += val(s, m.time, m.dropout_pfs)
        if m.pfs_switch == "censor":
            h += val(s, m.time, m.switch_pre)
    return h


# Markov models: generator and matrix exponentials -----------------------------------------------------------
ORDER = (1, 5, 2, 3, 4)


def is_markov(m):
    return m.clock == "forward" and m.mode() == "same"


def generator(m, t):
    """Generator at time t of the states in ORDER (Markov models only)."""
    a = val(t, m.time, m.prog); b = val(t, m.time, m.death); p = val(t, m.time, m.switch_prob)
    e = val(t, m.time, m.switch_pre); a5 = val(t, m.time, m.prog_switch); b5 = val(t, m.time, m.death_switch_pre)
    r = val(t, m.time_post, m.switch_rate); c = val(t, m.time_post, m.death_prog)
    d = val(t, m.time_post, m.death_switch)
    return np.array([[-(a + b + e), e, a * (1 - p), a * p, b],
                     [0.0, -(a5 + b5), 0.0, a5, b5],
                     [0.0, 0.0, -(r + c), r, c],
                     [0.0, 0.0, 0.0, -d, d],
                     [0.0, 0.0, 0.0, 0.0, 0.0]])


def trans(m, s, u):
    """Transition matrix from time s to time u >= s (Markov models only)."""
    br = sorted(set([s, u] + [x for x in list(m.time) + list(m.time_post) if s < x < u]))
    P = np.eye(5)
    for lo, hi in zip(br[:-1], br[1:]):
        P = P @ linalg.expm(generator(m, 0.5 * (lo + hi)) * (hi - lo))
    return P


def markov_quantities(m, s):
    """States, OS survival and density, PFS at risk and event density at s by matrix exponentials."""
    pi = trans(m, 0.0, s)[0] if s > 0 else np.eye(5)[0]
    Q = generator(m, s)
    os_d = pi[0] * Q[0, 4] + pi[1] * Q[1, 4] + pi[2] * Q[2, 4] + pi[3] * Q[3, 4]
    if m.pfs_switch == "censor":
        pfs_s, pfs_d = pi[0], pi[0] * (Q[0, 2] + Q[0, 3] + Q[0, 4])
    else:
        pfs_s = pi[0] + pi[1]
        pfs_d = pi[0] * (Q[0, 2] + Q[0, 3] + Q[0, 4]) + pi[1] * (Q[1, 3] + Q[1, 4])
    return dict(states=pi, os_surv=1 - pi[4], os_dens=os_d, pfs_surv=pfs_s, pfs_dens=pfs_d)


def from5_markov(m, w, u):
    P = trans(m, w, u)[1] if u > w else np.eye(5)[1]
    Q = generator(m, u)
    return 1 - P[4], P[1] * Q[1, 4] + P[3] * Q[3, 4]


# patient-level simulation -------------------------------------------------------------------------------------
def inv_cum(target, start, brk, rate):
    """Smallest x >= start with integral from start to x of the step function equal to target (inf if none)."""
    target = np.asarray(target, float)
    x = np.array(start, float, copy=True)
    left = target.copy()
    out = np.full(x.size, np.inf)
    active = np.ones(x.size, bool)
    while np.any(active):
        j = np.searchsorted(brk, x, side="right") - 1
        j = np.maximum(j, 0)
        nxt = np.where(j + 1 < brk.size, brk[np.minimum(j + 1, brk.size - 1)], np.inf)
        rt = rate[j]
        with np.errstate(invalid="ignore"):
            room = np.where(np.isfinite(nxt), rt * (nxt - x), np.inf)
        hit = active & (left <= room) & (rt > 0)
        out[hit] = x[hit] + left[hit] / rt[hit]
        active &= ~hit
        last = active & ~np.isfinite(nxt)
        active &= ~last
        move = active
        left[move] -= room[move]
        x[move] = nxt[move]
    return out


def step(x, brk, v):
    return v[np.maximum(np.searchsorted(brk, x, side="right") - 1, 0)]


def simulate(m, n, rng):
    """PFS time, PFS event indicator and OS time of n patients (no dropout)."""
    tm = m.time
    tp = m.time_post
    t1 = inv_cum(rng.exponential(size=n), np.zeros(n), tm, m.prog + m.death + m.switch_pre)
    a = step(t1, tm, m.prog); b = step(t1, tm, m.death); e = step(t1, tm, m.switch_pre)
    uu = rng.uniform(size=n) * (a + b + e)
    is_prog = uu < a
    is_switch = uu >= a + b
    pfs = t1.copy(); ev = ~is_switch if m.pfs_switch == "censor" else np.ones(n, bool)
    os_ = t1.copy()
    # progression from state 1
    k = np.where(is_prog)[0]
    v = t1[k]
    origin = v if m.clock == "forward" else np.zeros(k.size)
    now = rng.uniform(size=k.size) < step(v, tm, m.switch_prob)
    tl = inv_cum(rng.exponential(size=k.size), origin, tp, m.switch_rate + m.death_prog)
    r = step(tl, tp, m.switch_rate); c = step(tl, tp, m.death_prog)
    later = ~now & (rng.uniform(size=k.size) * np.maximum(r + c, 1e-300) < r)
    sig = np.where(now, origin, tl)
    end = tl.copy()
    swi = np.where(now | later)[0]
    x = rng.exponential(size=swi.size)
    s0 = sig[swi]
    mode = m.mode()
    if mode == "same":
        end[swi] = inv_cum(x, s0, tp, m.death_switch)
    elif mode == "switch":
        end[swi] = s0 + inv_cum(x, np.zeros(swi.size), tp, m.death_switch)
    else:
        z = inv_cum(x, s0, tp, m.death_prog)
        end[swi] = s0 + (z - s0) / m.switch_af
    os_[k] = v + (end - origin)
    # switch before progression: state 5 from t1
    q = np.where(is_switch)[0]
    if q.size:
        w = t1[q]
        h5 = m.prog_switch + m.death_switch_pre
        t5 = inv_cum(rng.exponential(size=q.size), w, tm, h5)
        a5 = step(t5, tm, m.prog_switch); b5 = step(t5, tm, m.death_switch_pre)
        prog5 = rng.uniform(size=q.size) * np.maximum(a5 + b5, 1e-300) < a5
        if m.pfs_switch != "censor":
            pfs[q] = t5
        # after progression from state 5 at t5: death on the post-progression clock from oc
        oc = np.where(prog5, t5 if m.clock == "forward" else 0.0, 0.0)
        dd = inv_cum(rng.exponential(size=q.size), oc, tp, m.death_switch)
        os_[q] = np.where(prog5, t5 + (dd - oc), t5)
    return pfs, ev, os_


def draw_censoring(m, n, rng):
    """Dropout (both endpoints) and censoring of PFS only."""
    c = inv_cum(rng.exponential(size=n), np.zeros(n), m.time, m.dropout)
    cp = inv_cum(rng.exponential(size=n), np.zeros(n), m.time, m.dropout_pfs)
    return c, cp
