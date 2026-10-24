"""Reference computation of the products of the PFS and OS functionals for the model of ref_model5.py.

Written independently of the R code (.cross_group(), addendum of the covariance note). For a patient of group z
the functionals are r = int rho dN + int beta Y ds over the observed follow-up, with rho and beta of the score
(psi), the observed variance (phi) and the event indicator (dlt) as in ref_zeta.py. With the entry time E, the
dropout C and the censoring of PFS only C_P independent of the event process,
  E[r_P(t_p) r_O(t_o)] = int int omega(s, u) E[(rho_P dN_P(s) + beta_P Y_P(s) ds) (rho_O dN_O(u) + beta_O Y_O(u) du)],
  omega(s, u) = Abar(min(t_p - s, t_o - u)) G(max(s, u)) G_P(s).
For a Markov model (clock "forward", switch_clock "same"), the expectations of the products of the at-risk
indicators and the counting-process increments are written with the state probabilities pi(s), the generator
Q(s) and the transition matrices P(s, u) (matrix exponentials). The double integral is computed by nested
adaptive quadrature (scipy quad_vec) over u and s; the deaths without progression, which are PFS and OS events at
the same time, give a one-dimensional term. The R code instead uses closed forms of the joint law split into the
terms C1 to C4 of the covariance note and Gauss-Legendre rules.
"""
import numpy as np
from scipy import integrate

import ref_logrank as rl
import ref_model5 as r5

QV = dict(epsabs=1e-12, epsrel=1e-9, limit=200)
I1, I5, I2, I3, I4 = range(5)


def _pfs_sets(m):
    """States at risk for PFS and, for each, the target states of a PFS event."""
    if m.pfs_switch == "censor":
        return {I1: (I2, I3, I4)}
    return {I1: (I2, I3, I4), I5: (I3, I4)}


def group_funcs(m, s, ep):
    """At-risk probability R and observed event density e of one group (with censoring)."""
    if r5.is_markov(m):
        q = r5.markov_quantities(m, s)
        S, f = (q["pfs_surv"], q["pfs_dens"]) if ep == "pfs" else (q["os_surv"], q["os_dens"])
    elif ep == "pfs":
        S, f = r5.pfs_surv(m, s), r5.pfs_dens(m, s)
    else:
        S, f = r5.os_surv(m, s), r5.os_dens(m, s)
    G = r5.dropout_surv(m, s, ep)
    return S * G, f * G


def fun_values(models, p, s, ep, z):
    """rho and beta of the three functionals for group z at follow-up time s."""
    (R0, e0), (R1, e1) = group_funcs(models[0], s, ep), group_funcs(models[1], s, ep)
    ys = p[0] * R0 + p[1] * R1
    q = p[1] * R1 / ys
    lbar = (p[0] * e0 + p[1] * e1) / ys
    g = 1 - q if z == 1 else -q
    return np.array([g, q * (1 - q), 1.0]), np.array([-g * lbar, (1 - 2 * q) * g * lbar, 0.0])


def cross_group(models, p, z, alpha, rho, t_p, t_o):
    """E[r_P(t_p) r_O(t_o)] (rows PFS, columns OS) for a patient of group z with a Markov model."""
    m = models[z]
    if not r5.is_markov(m):
        raise ValueError("Markov models only")
    risk = _pfs_sets(m)
    brk = sorted(set(list(m.time) + list(m.time_post)))
    G = lambda x: r5.dropout_surv(m, x, "os")
    GP = lambda x: np.exp(-r5.cum(x, m.time, m.dropout_pfs))
    omega = lambda s, u: rl.abar(alpha, rho, min(t_p - s, t_o - u)) * G(max(s, u)) * GP(s)
    alive = np.array([1.0, 1.0, 1.0, 1.0, 0.0])

    def inner(s, u, rho_o, beta_o):
        rho_p, beta_p = fun_values(models, p, s, "pfs", z)
        pi = r5.trans(m, 0.0, s)[0] if s > 0 else np.eye(5)[0]
        Q = r5.generator(m, s)
        w = omega(s, u)
        if w == 0.0:
            return np.zeros(9)
        if s < u:
            P = r5.trans(m, s, u)
            Qu = r5.generator(m, u)
            death_u = P @ Qu[:, I4]               # density of death at u from each state at s
            surv_u = P @ alive                    # probability of being alive at u from each state at s
            k_bb = sum(pi[i] * surv_u[i] for i in risk)
            k_bd = sum(pi[i] * death_u[i] for i in risk)
            k_db = sum(pi[i] * sum(Q[i, k] * surv_u[k] for k in risk[i]) for i in risk)
            k_dd = sum(pi[i] * sum(Q[i, k] * death_u[k] for k in risk[i] if k != I4) for i in risk)
            out = (k_dd * np.outer(rho_p, rho_o) + k_db * np.outer(rho_p, beta_o)
                   + k_bd * np.outer(beta_p, rho_o) + k_bb * np.outer(beta_p, beta_o))
        else:
            f_p = sum(pi[i] * sum(Q[i, k] for k in risk[i]) for i in risk)
            s_p = sum(pi[i] for i in risk)
            out = f_p * np.outer(rho_p, beta_o) + s_p * np.outer(beta_p, beta_o)
        return w * out.ravel()

    def outer(u):
        rho_o, beta_o = fun_values(models, p, u, "os", z)
        hi = t_p
        cands = [u, t_p - t_o + u] + [t_p - a for a in alpha] + brk
        pts = sorted(set(x for x in cands if 0 < x < hi)) or None
        return integrate.quad_vec(lambda s: inner(s, u, rho_o, beta_o), 0.0, hi, points=pts, **QV)[0]
    cands = [t_p] + [t_o - a for a in alpha] + [t_p - a for a in alpha] + brk + [x + t_o - t_p for x in brk]
    pts = sorted(set(x for x in cands if 0 < x < t_o)) or None
    res = integrate.quad_vec(outer, 0.0, t_o, points=pts, **QV)[0].reshape(3, 3)

    def diag(x):
        pi = r5.trans(m, 0.0, x)[0] if x > 0 else np.eye(5)[0]
        Q = r5.generator(m, x)
        dens = sum(pi[i] * Q[i, I4] for i in risk)
        rp, _ = fun_values(models, p, x, "pfs", z)
        ro, _ = fun_values(models, p, x, "os", z)
        return omega(x, x) * dens * np.outer(rp, ro).ravel()
    tmin = min(t_p, t_o)
    pts = sorted(set(x for x in brk + [tmin - a for a in alpha] if 0 < x < tmin)) or None
    res += integrate.quad_vec(diag, 0.0, tmin, points=pts, **QV)[0].reshape(3, 3)
    return res


def cross_moments(models, p, alpha, rho, t_p, t_o):
    """Per-patient E[r_P(t_p) r_O(t_o)] over both groups."""
    return sum(p[z] * cross_group(models, p, z, alpha, rho, t_p, t_o) for z in (0, 1))
