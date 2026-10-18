"""Reference implementation of the per-group multistate model.

This module is written independently of the R code in DesignComplexSurvival.
It uses scalar loops, matrix exponentials (scipy.linalg.expm) and adaptive
quadrature (scipy.integrate.quad), whereas the R code uses closed-form
expressions and Gauss-Legendre rules.

States: 1 progression-free and alive, 2 progressed without switch,
3 switched, 4 dead.
"""

import numpy as np
from scipy import integrate, linalg

QUAD_OPTS = dict(epsabs=1e-13, epsrel=1e-12, limit=400)


def _expand(x, n):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    if x.size == 1:
        return np.repeat(x, n)
    if x.size != n:
        raise ValueError("length mismatch")
    return x


class Model:
    """One treatment group.

    Pre-progression intensities (prog, death, switch_prob, dropout) are
    piecewise constant on `time` (time since randomization).
    Post-progression intensities (death_prog, switch_rate, death_switch) are
    piecewise constant on `time_post`, measured from randomization when
    clock = "forward" and from progression when clock = "reset".
    """

    def __init__(self, death, time=(0.0,), prog=0.0, death_prog=0.0,
                 switch_prob=0.0, switch_rate=0.0, death_switch=None,
                 dropout=0.0, time_post=(0.0,), clock="reset"):
        self.time = np.asarray(time, dtype=float)
        n = self.time.size
        self.prog = _expand(prog, n)
        self.death = _expand(death, n)
        self.switch_prob = _expand(switch_prob, n)
        self.dropout = _expand(dropout, n)
        self.time_post = np.asarray(time_post, dtype=float)
        m = self.time_post.size
        self.death_prog = _expand(death_prog, m)
        self.switch_rate = _expand(switch_rate, m)
        if death_switch is None:
            self.death_switch = self.death_prog.copy()
        else:
            self.death_switch = _expand(death_switch, m)
        self.clock = clock


def pw_value(x, time, val):
    """Value of a right-continuous step function at x."""
    i = int(np.searchsorted(time, x, side="right")) - 1
    return float(val[max(i, 0)])


def pw_cum(x, time, val):
    """Integral from 0 to x of a step function (loop over intervals)."""
    x = max(float(x), 0.0)
    total = 0.0
    for j in range(time.size):
        lo = time[j]
        hi = time[j + 1] if j + 1 < time.size else np.inf
        if x <= lo:
            break
        total += val[j] * (min(x, hi) - lo)
    return total


def _post_generator(r, c, d):
    return np.array([[-(r + c), r, c],
                     [0.0, -d, d],
                     [0.0, 0.0, 0.0]])


def post_R_reset(m, v):
    """3 x 3 transition matrix of (progressed, switched, dead) over a sojourn v."""
    P = np.eye(3)
    tp = m.time_post
    for j in range(tp.size):
        lo = tp[j]
        hi = tp[j + 1] if j + 1 < tp.size else np.inf
        if v <= lo:
            break
        dt = min(v, hi) - lo
        Q = _post_generator(m.switch_rate[j], m.death_prog[j], m.death_switch[j])
        P = P @ linalg.expm(Q * dt)
    return P


def post_R_forward(m, u, s):
    """3 x 3 transition matrix of the post-progression process from time u to s."""
    pts = sorted(set([u, s] + [x for x in m.time_post if u < x < s]))
    P = np.eye(3)
    for lo, hi in zip(pts[:-1], pts[1:]):
        mid = 0.5 * (lo + hi)
        Q = _post_generator(pw_value(mid, m.time_post, m.switch_rate),
                            pw_value(mid, m.time_post, m.death_prog),
                            pw_value(mid, m.time_post, m.death_switch))
        P = P @ linalg.expm(Q * (hi - lo))
    return P


def forward_states_expm(m, t):
    """State occupation probabilities for the Markov (forward) model by expm."""
    if t <= 0:
        return np.array([1.0, 0.0, 0.0, 0.0])
    br = [x for x in np.concatenate([m.time, m.time_post]) if 0 < x < t]
    pts = sorted(set([0.0, t] + br))
    P = np.eye(4)
    for lo, hi in zip(pts[:-1], pts[1:]):
        mid = 0.5 * (lo + hi)
        a = pw_value(mid, m.time, m.prog)
        b = pw_value(mid, m.time, m.death)
        p = pw_value(mid, m.time, m.switch_prob)
        r = pw_value(mid, m.time_post, m.switch_rate)
        c = pw_value(mid, m.time_post, m.death_prog)
        d = pw_value(mid, m.time_post, m.death_switch)
        Q = np.array([[-(a + b), a * (1 - p), a * p, b],
                      [0.0, -(r + c), r, c],
                      [0.0, 0.0, -d, d],
                      [0.0, 0.0, 0.0, 0.0]])
        P = P @ linalg.expm(Q * (hi - lo))
    return P[0].copy()


def _p11(m, u):
    return np.exp(-pw_cum(u, m.time, m.prog + m.death))


def _post_R(m, u, s):
    if m.clock == "forward":
        return post_R_forward(m, u, s), s
    return post_R_reset(m, s - u), s - u


def _conv_points(m, t):
    if m.clock == "forward":
        br = list(m.time) + list(m.time_post)
    else:
        br = list(m.time) + [t - x for x in m.time_post]
    return sorted(set(x for x in br if 0 < x < t))


def conv_states(m, t):
    """State occupation probabilities by convolution over the progression time."""
    if t <= 0:
        return np.array([1.0, 0.0, 0.0, 0.0])

    def f2(u):
        R, _ = _post_R(m, u, t)
        a = pw_value(u, m.time, m.prog)
        p = pw_value(u, m.time, m.switch_prob)
        return _p11(m, u) * a * (1 - p) * R[0, 0]

    def f3(u):
        R, _ = _post_R(m, u, t)
        a = pw_value(u, m.time, m.prog)
        p = pw_value(u, m.time, m.switch_prob)
        return _p11(m, u) * a * ((1 - p) * R[0, 1] + p * R[1, 1])

    pts = _conv_points(m, t)
    p2 = integrate.quad(f2, 0, t, points=pts or None, **QUAD_OPTS)[0]
    p3 = integrate.quad(f3, 0, t, points=pts or None, **QUAD_OPTS)[0]
    p1 = _p11(m, t)
    return np.array([p1, p2, p3, 1.0 - p1 - p2 - p3])


def states(m, t):
    if m.clock == "forward":
        return forward_states_expm(m, t)
    return conv_states(m, t)


def os_surv(m, t):
    st = states(m, t)
    return st[0] + st[1] + st[2]


def os_dens(m, t):
    """Density of overall survival."""
    if t <= 0:
        return pw_value(0.0, m.time, m.death)
    if m.clock == "forward":
        st = forward_states_expm(m, t)
        return (st[0] * pw_value(t, m.time, m.death)
                + st[1] * pw_value(t, m.time_post, m.death_prog)
                + st[2] * pw_value(t, m.time_post, m.death_switch))

    def fd(u):
        R, v = _post_R(m, u, t)
        a = pw_value(u, m.time, m.prog)
        p = pw_value(u, m.time, m.switch_prob)
        c = pw_value(v, m.time_post, m.death_prog)
        d = pw_value(v, m.time_post, m.death_switch)
        return _p11(m, u) * a * ((1 - p) * (R[0, 0] * c + R[0, 1] * d) + p * R[1, 1] * d)

    pts = _conv_points(m, t)
    conv = integrate.quad(fd, 0, t, points=pts or None, **QUAD_OPTS)[0]
    return _p11(m, t) * pw_value(t, m.time, m.death) + conv


def pfs_surv(m, t):
    return _p11(m, t)


def pfs_dens(m, t):
    return _p11(m, t) * pw_value(t, m.time, m.prog + m.death)


def surv(m, t, endpoint):
    return os_surv(m, t) if endpoint == "os" else pfs_surv(m, t)


def dens(m, t, endpoint):
    return os_dens(m, t) if endpoint == "os" else pfs_dens(m, t)


def dropout_surv(m, t):
    return np.exp(-pw_cum(t, m.time, m.dropout))


def dropout_haz(m, t):
    return pw_value(t, m.time, m.dropout)


def model_breaks(m):
    return sorted(set(list(m.time) + list(m.time_post)))
