"""Expected values for the tests of ms_size() with the options of stage 3 (tests/testthat/test-ms_size.R).

Writes results/size_stage3.csv (columns design, quantity, value) for three designs with one look:
  recycling        scenario F of expected_zeta.py (n = 400); PFS and OS tested at one look driven by OS events;
                   Bonferroni procedure with weights 0.8 and 0.2 and full recycling of the level
                   (transition = c(1, 1)); variance "null". The time at which the smaller of the powers of PFS and
                   OS is 0.8 (target = "each"), the expected OS events at that time (the solution of
                   ms_size(solve_for = "events")) and, for these events rounded up, the time of the look and the
                   powers of PFS and OS and the probabilities of rejecting both and at least one hypothesis. As in
                   ms_size(), the search is over the time of the look, which is fixed, while the design with the
                   rounded events has a look driven by OS events, whose random time enters the correlation of the
                   statistics of PFS and OS (trigger "os" in ref_zeta.zeta_cov())
  calendar_n       scenario C with switching of the control group possible from month 15 (scn_cal() of
                   helper-scenarios.R, forward clock); one look at month 36; variance "null". The sample size with
                   power 0.8 and, for this sample size rounded up, the events and the power
  calendar_events  the same scenario with n = 1000 and one look driven by OS events. The time with power 0.8, the
                   expected OS events at that time and, for these events rounded up, the time and the power

Each quantity is computed independently of the R code:
  moments     expected_design.joint() (ref_logrank and ref_zeta: adaptive quadrature and an ODE solver), as for
              results/design2.csv, and ref_calendar.Design.look() (matrix exponentials and Gauss-Legendre rules
              over the whole entry range), as for results/calendar_looks.csv
  bounds      one look: the upper 1 - level quantile of the normal distribution at the levels 0.8 alpha, 0.2 alpha
              and alpha
  power       recycling: the decisions of Algorithm 1 of Maurer and Bretz (2013) at one look written as
              one-dimensional integrals of the bivariate normal density (adaptive quadrature); one endpoint: the
              normal distribution function
  solutions   root finding (brentq) over the time of the look; closed form for the sample size
Checks printed at the end: the powers of the recycling design by ref_graph.recursion(), the relation
any = PFS + OS - both, the moments of the calendar design at month 36 against results/calendar_looks.csv, and the
power of OS without recycling at the time of the solution.

Usage (from this folder; about 2 minutes): python3 expected_size_stage3.py
"""
import csv
import math
import os
import sys
import time

import numpy as np
from scipy import integrate, optimize, stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import expected_calendar as EC  # noqa: E402
import expected_design as ED  # noqa: E402
import ref_graph  # noqa: E402
import ref_logrank as rl  # noqa: E402

OUT = os.path.join(HERE, "results")
ALPHA = 0.025
POWER = 0.8
QUAD = dict(epsabs=1e-14, epsrel=1e-12, limit=200)

# design "recycling": levels 0.8 alpha (PFS), 0.2 alpha (OS) and alpha after recycling
C0, C1, C2 = stats.norm.isf(ALPHA), stats.norm.isf(0.8 * ALPHA), stats.norm.isf(0.2 * ALPHA)
D_REC = dict(endpoints=["pfs", "os"], trigger=["os"], fraction=[1.0], time=[None], test=None, events=None,
             followup=0.0)
N_REC = ED.N
N_CAL_EVENTS = 1000.0


def upper_box(a, b, c, m1, m2, r):
    """P(a < X <= b, Y > c) for unit variances, means m1 and m2 and correlation r."""
    s = math.sqrt(1.0 - r * r)
    hi = min(b, m1 + 12.0)
    lo = max(a, m1 - 12.0)
    if hi <= lo:
        return 0.0
    f = lambda x: stats.norm.pdf(x - m1) * stats.norm.sf((c - m2 - r * (x - m1)) / s)
    return integrate.quad(f, lo, hi, **QUAD)[0]


def rec_joint(t, fixed=True):
    """Moments at one look at time t: fixed time (fixed = True) or driven by OS events."""
    mean, sig, mom = ED.joint(ED.RATES["F"], D_REC, np.array([t]), {"pfs": [0], "os": [0]}, fixed, "null")
    return mean, sig, mom


def rec_probs(mean, r):
    """Rejection probabilities at one look with full recycling (Maurer and Bretz 2013, Algorithm 1): PFS is
    rejected if Z_P > c1, or if Z_O > c2 and Z_P > c0; OS likewise with the roles exchanged."""
    mp, mo = float(mean[0]), float(mean[1])
    pfs = stats.norm.sf(C1 - mp) + upper_box(C0, C1, C2, mp, mo, r)
    os_ = stats.norm.sf(C2 - mo) + upper_box(C0, C2, C1, mo, mp, r)
    both = upper_box(C1, np.inf, C0, mp, mo, r) + upper_box(C0, C1, C2, mp, mo, r)
    p00 = stats.norm.cdf(C1 - mp) - upper_box(-np.inf, C1, C2, mp, mo, r)
    return dict(pfs=pfs, os=os_, both=both, any=1.0 - p00)


def rec_at(t, fixed=True):
    mean, sig, mom = rec_joint(t, fixed)
    pr = rec_probs(mean, sig[0, 1])
    return pr, mean, sig[0, 1], mom


def recycling(rows, log):
    t0 = time.time()
    def f(t):
        pr = rec_at(t)[0]
        return min(pr["pfs"], pr["os"]) - POWER
    lo, hi = 30.0, 36.0
    flo, fhi = f(lo), f(hi)
    log.append("recycling: target - power at months %.0f and %.0f: %.6f, %.6f" % (lo, hi, flo, fhi))
    t_sol = optimize.brentq(f, lo, hi, xtol=1e-12, rtol=1e-15)
    pr, mean, r, mom = rec_at(t_sol)
    d_sol = mom["os"][0]["events"]
    scn = ED.scen(ED.RATES["F"])
    d_up = math.ceil(d_sol - 1e-8)
    t_up = rl.event_time(scn, N_REC, d_up, "os")
    pu, mean_u, r_u, mom_u = rec_at(t_up, fixed=False)
    pf, _, r_f, _ = rec_at(t_up)
    rows += [("recycling", "time_solution", t_sol), ("recycling", "solution", d_sol),
             ("recycling", "events_rounded", float(d_up)), ("recycling", "time_rounded", t_up),
             ("recycling", "power_pfs", pu["pfs"]), ("recycling", "power_os", pu["os"]),
             ("recycling", "power_both", pu["both"]), ("recycling", "power_any", pu["any"])]
    log.append("recycling: solution at month %.10f, OS events %.10f (PFS %.6f, OS %.10f); rounded %d at month %.10f"
               % (t_sol, d_sol, pr["pfs"], pr["os"], d_up, t_up))
    log.append("  rounded design: power PFS %.10f, OS %.10f, both %.10f, any %.10f; Z means %s, correlation %.10f"
               % (pu["pfs"], pu["os"], pu["both"], pu["any"], np.array2string(mean_u, precision=10), r_u))
    log.append("  (look at the same time but fixed: correlation %.10f, power PFS %.10f, OS %.10f)"
               % (r_f, pf["pfs"], pf["os"]))
    log.append("  check: events at the rounded time %.12f (target %d)" % (mom_u["os"][0]["events"], d_up))
    log.append("  check: any - (PFS + OS - both) = %.1e" % (pu["any"] - (pu["pfs"] + pu["os"] - pu["both"])))
    cfg = dict(tau1=[1.0], tau2=[1.0], rho=r_u, theta1=float(mean_u[0]), theta2=float(mean_u[1]), a1=[1], a2=[1],
               u1=[C1], u2=[C2], r1=[C0], r2=[C0], l1=None)
    g = ref_graph.recursion(cfg)
    dif = max(abs(g["first"] - pu["pfs"]), abs(g["second"] - pu["os"]), abs(g["both"] - pu["both"]),
              abs(g["any"] - pu["any"]))
    log.append("  check: ref_graph.recursion() against the integrals, max abs difference %.1e" % dif)
    bonf_os = stats.norm.sf(C2 - mean[1])
    log.append("  without recycling at the solution: power of OS %.6f (with recycling %.6f)" % (bonf_os, pr["os"]))
    log.append("  (%.0f s)" % (time.time() - t0))


def calendar(rows, log):
    t0 = time.time()
    des = EC.design()
    z = stats.norm.isf(ALPHA)
    zp = stats.norm.ppf(POWER)
    look = {}

    def at(t):
        if t not in look:
            L = des.look(t)
            look[t] = (-L["mu"] / math.sqrt(L["v"]), sum(L["ev"]), L["mu"], L["v"])
        return look[t]

    # calendar_n: one look at month 36
    th, ev, mu, v = at(36.0)
    with open(os.path.join(OUT, "calendar_looks.csv"), newline="") as fh:
        ref = [r for r in csv.DictReader(fh) if r["design"] == "alt" and float(r["t"]) == 36.0][0]
    log.append("calendar: month 36, relative differences from calendar_looks.csv: mu %.1e, v %.1e, events %.1e"
               % (mu / float(ref["mu"]) - 1, v / float(ref["v"]) - 1,
                  ev / (float(ref["ev0"]) + float(ref["ev1"])) - 1))
    n_sol = ((z + zp) / th) ** 2
    n_up = math.ceil(n_sol - 1e-8)
    p_up = stats.norm.sf(z - math.sqrt(n_up) * th)
    rows += [("calendar_n", "solution", n_sol), ("calendar_n", "n_rounded", float(n_up)),
             ("calendar_n", "events_rounded", n_up * ev), ("calendar_n", "power", p_up)]
    log.append("calendar_n: n %.10f, rounded %d, events %.10f, power %.12f" % (n_sol, n_up, n_up * ev, p_up))

    # calendar_events: n = 1000, one look driven by OS events
    n = N_CAL_EVENTS
    f = lambda t: stats.norm.sf(z - math.sqrt(n) * at(t)[0]) - POWER
    lo, hi = 18.0, 24.0
    log.append("calendar_events: power - target at months %.0f and %.0f: %.6f, %.6f" % (lo, hi, f(lo), f(hi)))
    t_sol = optimize.brentq(f, lo, hi, xtol=1e-12, rtol=1e-15)
    d_sol = n * at(t_sol)[1]
    d_up = math.ceil(d_sol - 1e-8)
    t_up = optimize.brentq(lambda t: n * at(t)[1] - d_up, t_sol, hi, xtol=1e-12, rtol=1e-15)
    th_u, ev_u, _, _ = at(t_up)
    p_up = stats.norm.sf(z - math.sqrt(n) * th_u)
    rows += [("calendar_events", "time_solution", t_sol), ("calendar_events", "solution", d_sol),
             ("calendar_events", "events_rounded", float(d_up)), ("calendar_events", "time_rounded", t_up),
             ("calendar_events", "power", p_up)]
    log.append("calendar_events: solution at month %.10f, OS events %.10f; rounded %d at month %.10f, power %.12f"
               % (t_sol, d_sol, d_up, t_up, p_up))
    log.append("  check: events at the rounded time %.12f (target %d)" % (n * ev_u, d_up))
    log.append("  (%d evaluations of Design.look(), %.0f s)" % (len(look), time.time() - t0))


def main():
    t_all = time.time()
    rows, log = [], []
    recycling(rows, log)
    print("\n".join(log), flush=True)
    log2 = []
    calendar(rows, log2)
    print("\n".join(log2), flush=True)
    path = os.path.join(OUT, "size_stage3.csv")
    with open(path, "w", newline="\n") as fh:
        fh.write("design,quantity,value\n")
        for d, q, v in rows:
            fh.write("%s,%s,%s\n" % (d, q, repr(float(v))))
    print("wrote %s (%d rows); total %.0f s" % (path, len(rows), time.time() - t_all))


if __name__ == "__main__":
    main()
