"""Expected values for the tests of gs_prob_mtp() with upper_recycled (graphical procedure with recycling of the
level, Maurer and Bretz 2013), in tests/testthat/test-gs_prob_mtp.R and test-mtp_graph.R.

Writes
  results/gs_prob_graph_configs.csv  the configurations (canonical joint distribution, bounds at the initial level
                                     u1, u2 and after recycling r1, r2, futility bounds l1, analyses)
  results/gs_prob_graph.csv          rejection and futility probabilities by ref_graph.recursion()
  results/gs_prob_graph_mc.csv       Monte Carlo check of the recursion (ref_graph.simulate(), 2 million trials)

Checks of the reference itself, printed at the end: with r = u the recursion equals ref_mtp.recursion() for the
Bonferroni procedure, and with u2 = Inf and r1 = u1 it equals ref_mtp.recursion() for the fixed sequence 'overall'.

Usage: python3 expected_graph.py
"""
import csv
import os
import time

import numpy as np

import ref_graph
import ref_mtp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
INF = float("inf")

CONFIGS = [
    # three looks of both endpoints at the same analyses, futility bounds for the first endpoint (the last value
    # is not used)
    dict(config="three_looks", tau1=[0.4, 0.7, 1.0], tau2=[0.4, 0.7, 1.0], rho=0.6, theta1=3.0, theta2=2.5,
         a1=[1, 2, 3], a2=[1, 2, 3], u1=[3.0, 2.4, 2.0], u2=[2.6, 2.3, 2.1], r1=[2.7, 2.1, 1.8],
         r2=[2.3, 2.0, 1.9], l1=[0.0, 1.0, 2.0]),
    # the first endpoint at analyses 1 and 2, the second at analyses 1, 2 and 3 (tau1 at analysis 3 is not used)
    dict(config="unequal", tau1=[0.6, 1.0, 1.0], tau2=[0.3, 0.6, 1.0], rho=0.7, theta1=2.8, theta2=2.2,
         a1=[1, 2], a2=[1, 2, 3], u1=[2.5, 2.0], u2=[3.2, 2.6, 2.0], r1=[2.2, 1.9], r2=[2.9, 2.3, 1.8],
         l1=[0.5, 2.0]),
    # the second endpoint is not tested at analysis 1; no futility bounds
    dict(config="late_second", tau1=[0.5, 0.75, 1.0], tau2=[0.25, 0.55, 1.0], rho=0.5, theta1=3.2, theta2=2.0,
         a1=[1, 2, 3], a2=[2, 3], u1=[2.9, 2.3, 2.0], u2=[2.4, 2.1], r1=[2.6, 2.1, 1.9], r2=[2.2, 1.95],
         l1=None),
    # a futility bound above the bound after recycling at look 2 of the first endpoint (2.2 > 2.1), and no test of
    # the second hypothesis at its initial level at look 1 (infinite bound)
    dict(config="lower_above", tau1=[0.4, 0.7, 1.0], tau2=[0.4, 0.7, 1.0], rho=0.6, theta1=2.6, theta2=2.5,
         a1=[1, 2, 3], a2=[1, 2, 3], u1=[3.0, 2.4, 2.0], u2=[INF, 2.3, 2.1], r1=[2.7, 2.1, 1.8],
         r2=[2.5, 2.0, 1.9], l1=[0.5, 2.2, 2.0]),
    # the first endpoint is tested from analysis 2 on, with a futility bound there
    dict(config="first_late", tau1=[0.3, 0.6, 1.0], tau2=[0.4, 0.7, 1.0], rho=0.5, theta1=2.6, theta2=2.8,
         a1=[2, 3], a2=[1, 2, 3], u1=[2.6, 2.0], u2=[3.0, 2.5, 2.1], r1=[2.3, 1.9], r2=[2.6, 2.2, 1.9],
         l1=[0.8, 2.0]),
]
HYPS = [("alternative", 1.0, 1.0), ("null_first", 0.0, 1.0), ("null_second", 1.0, 0.0), ("null_both", 0.0, 0.0)]


def vec(x):
    """';'-separated values; infinite values as 'Inf' (read by as.numeric() in R)."""
    f = lambda v: repr(float(v)) if np.isfinite(v) else ("Inf" if v > 0 else "-Inf")
    return "" if x is None else ";".join(f(v) for v in x)


def quantities(res, k1, k2):
    q = [("first", res["first"]), ("second", res["second"]), ("both", res["both"]), ("any", res["any"])]
    q += [("p1_%d" % (j + 1), res["p1"][j]) for j in range(k1)]
    q += [("pf_%d" % (j + 1), res["pf"][j]) for j in range(k1)]
    q += [("p2_%d" % (i + 1), res["p2"][i]) for i in range(k2)]
    return q


def self_checks():
    """Reductions of the reference to ref_mtp.recursion() (maximum absolute differences)."""
    out = []
    for c in CONFIGS:
        k1 = len(c["a1"])
        lo = None if c["l1"] is None else c["l1"]
        for lab, th1, th2 in (("alternative", c["theta1"], c["theta2"]), ("null_both", 0.0, 0.0)):
            base = dict(c, theta1=th1, theta2=th2, l1=lo)
            u2 = [v if np.isfinite(v) else r for v, r in zip(c["u2"], c["r2"])]
            a = ref_mtp.recursion(dict(base, u2=u2, mtp="bonferroni"))
            b = ref_graph.recursion(dict(base, u2=u2, r1=c["u1"], r2=u2))
            f = ref_mtp.recursion(dict(base, u2=c["r2"], mtp="fixed_sequence", hierarchical="overall"))
            g = ref_graph.recursion(dict(base, u2=[INF] * len(c["a2"]), r1=c["u1"], r2=c["r2"]))
            d_b = max(np.abs(a[k] - b[k]).max() for k in ("p1", "pf", "p2", "both"))
            d_f = max(np.abs(f[k] - g[k]).max() for k in ("p1", "pf", "p2", "both"))
            out.append((c["config"], lab, d_b, d_f))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "gs_prob_graph_configs.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["config", "tau1", "tau2", "rho", "theta1", "theta2", "a1", "a2", "u1", "u2", "r1", "r2", "l1"])
        for c in CONFIGS:
            w.writerow([c["config"], vec(c["tau1"]), vec(c["tau2"]), repr(float(c["rho"])),
                        repr(float(c["theta1"])), repr(float(c["theta2"])), vec(c["a1"]), vec(c["a2"]),
                        vec(c["u1"]), vec(c["u2"]), vec(c["r1"]), vec(c["r2"]), vec(c["l1"])])
    rows, mc_rows = [], []
    for c in CONFIGS:
        k1, k2 = len(c["a1"]), len(c["a2"])
        for hyp, f1, f2 in HYPS:
            th1, th2 = f1 * c["theta1"], f2 * c["theta2"]
            cfg = dict(c, theta1=th1, theta2=th2)
            t0 = time.time()
            res = ref_graph.recursion(cfg)
            for name, val in quantities(res, k1, k2):
                rows.append([c["config"], hyp, repr(float(th1)), repr(float(th2)), "recursion", name,
                             repr(float(val))])
            sim = ref_graph.simulate(cfg, 2_000_000, 20261024)
            for (name, val), (_, sv) in zip(quantities(res, k1, k2), quantities(sim, k1, k2)):
                se = sim["se"](sv)
                mc_rows.append([c["config"], hyp, name, repr(float(val)), repr(float(sv)), repr(float(se)),
                                repr(float((sv - val) / se)) if se > 0 else ""])
            print("%-12s %-12s %.1fs" % (c["config"], hyp, time.time() - t0), flush=True)
    with open(os.path.join(OUT, "gs_prob_graph.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["config", "hypothesis", "theta1", "theta2", "method", "quantity", "value"])
        w.writerows(rows)
    with open(os.path.join(OUT, "gs_prob_graph_mc.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["config", "hypothesis", "quantity", "recursion", "simulation", "se", "z"])
        w.writerows(mc_rows)
    z = np.array([float(r[-1]) for r in mc_rows if r[-1] != ""])
    print("written: gs_prob_graph_configs.csv, gs_prob_graph.csv (%d rows), gs_prob_graph_mc.csv (%d rows, "
          "max |z| %.2f, |z| > 4: %d)" % (len(rows), len(mc_rows), np.abs(z).max(), int(np.sum(np.abs(z) > 4))))
    print("self-checks of the reference (maximum absolute difference to ref_mtp.recursion()):")
    for cf, lab, d_b, d_f in self_checks():
        print("  %-12s %-12s r = u vs Bonferroni %.1e; u2 = Inf, r1 = u1 vs fixed sequence (overall) %.1e"
              % (cf, lab, d_b, d_f))


if __name__ == "__main__":
    main()
