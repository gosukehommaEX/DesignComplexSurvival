"""Expected values for tests/testthat/test-gs_prob_mtp.R.

Writes
  results/gs_prob_mtp_configs.csv  the configurations (canonical joint distribution, bounds, analyses)
  results/gs_prob_mtp.csv          rejection probabilities by ref_mtp.recursion() and, for two looks,
                                   by the formulas of Glimm et al. (2010) (ref_mtp.glimm_two_look())
  results/gs_prob_mtp_mc.csv       Monte Carlo check of the recursion (ref_mtp.simulate(), 2 million trials)

Usage: python3 expected_mtp.py
"""
import csv
import os
import time

import numpy as np
from scipy import stats

import ref_mtp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")

DP = stats.norm.ppf(0.975) + stats.norm.ppf(0.92)
DS = stats.norm.ppf(0.975) + stats.norm.ppf(0.86)

CONFIGS = [
    # three looks of both endpoints at the same analyses, futility bounds for the first endpoint
    # (the last value is not used)
    dict(config="three_looks", tau1=[0.4, 0.7, 1.0], tau2=[0.4, 0.7, 1.0], rho=0.6, theta1=3.0,
         theta2=2.5, a1=[1, 2, 3], a2=[1, 2, 3], u1=[3.0, 2.4, 2.0], u2=[2.6, 2.3, 2.1], l1=[0.0, 1.0, 2.0]),
    # the first endpoint at analyses 1 and 2, the second at analyses 1, 2 and 3 with other information
    # fractions (tau1 at analysis 3 is not used)
    dict(config="unequal", tau1=[0.6, 1.0, 1.0], tau2=[0.3, 0.6, 1.0], rho=0.7, theta1=2.8, theta2=2.2,
         a1=[1, 2], a2=[1, 2, 3], u1=[2.5, 2.0], u2=[3.2, 2.6, 2.0], l1=[0.5, 2.0]),
    # the second endpoint is not tested at analysis 1; no futility bounds
    dict(config="late_second", tau1=[0.5, 0.75, 1.0], tau2=[0.25, 0.55, 1.0], rho=0.5, theta1=3.2,
         theta2=2.0, a1=[1, 2, 3], a2=[2, 3], u1=[2.9, 2.3, 2.0], u2=[2.4, 2.1], l1=None),
    # Glimm et al. (2010), Section 5: classical O'Brien-Fleming bounds for the primary and Pocock bounds
    # for the secondary endpoint (one-sided 0.025, interim at half of the information), drifts from
    # powers of 92 and 86 per cent; correlations 0 and 0.8
    dict(config="glimm_rho0", tau1=[0.5, 1.0], tau2=[0.5, 1.0], rho=0.0, theta1=DP, theta2=DS,
         a1=[1, 2], a2=[1, 2], u1=[2.7965, 1.9774], u2=[2.1783, 2.1783], l1=None),
    dict(config="glimm_rho08", tau1=[0.5, 1.0], tau2=[0.5, 1.0], rho=0.8, theta1=DP, theta2=DS,
         a1=[1, 2], a2=[1, 2], u1=[2.7965, 1.9774], u2=[2.1783, 2.1783], l1=None),
]
PROCS = [("fixed_sequence", "overall"), ("fixed_sequence", "stagewise"), ("bonferroni", "overall")]


def vec(x):
    return "" if x is None else ";".join(repr(float(v)) for v in x)


def quantities(res, k1, k2):
    q = [("first", res["first"]), ("second", res["second"]), ("both", res["both"]), ("any", res["any"])]
    q += [("p1_%d" % (j + 1), res["p1"][j]) for j in range(k1)]
    q += [("pf_%d" % (j + 1), res["pf"][j]) for j in range(k1)]
    q += [("p2_%d" % (i + 1), res["p2"][i]) for i in range(k2)]
    return q


def main():
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "gs_prob_mtp_configs.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["config", "tau1", "tau2", "rho", "theta1", "theta2", "a1", "a2", "u1", "u2", "l1"])
        for c in CONFIGS:
            w.writerow([c["config"], vec(c["tau1"]), vec(c["tau2"]), repr(float(c["rho"])),
                        repr(float(c["theta1"])), repr(float(c["theta2"])), vec(c["a1"]), vec(c["a2"]),
                        vec(c["u1"]), vec(c["u2"]), vec(c["l1"])])
    rows, mc_rows = [], []
    for c in CONFIGS:
        k1, k2 = len(c["a1"]), len(c["a2"])
        glimm = c["config"].startswith("glimm")
        hyps = [("alternative", c["theta1"], c["theta2"])]
        if not glimm:
            hyps += [("null_second", c["theta1"], 0.0), ("null_both", 0.0, 0.0)]
        for mtp, hier in PROCS:
            for hyp, th1, th2 in hyps:
                cfg = dict(c, theta1=th1, theta2=th2, mtp=mtp, hierarchical=hier)
                t0 = time.time()
                res = ref_mtp.recursion(cfg)
                for name, val in quantities(res, k1, k2):
                    rows.append([c["config"], mtp, hier, hyp, repr(float(th1)), repr(float(th2)), "recursion",
                                 name, repr(float(val))])
                if glimm and mtp == "fixed_sequence":
                    g = ref_mtp.glimm_two_look(c["u1"], c["u2"], th1, th2, c["rho"], c["tau1"][0])
                    gq = [("second", g["a"] if hier == "stagewise" else g["b"]),
                          ("p1_1", g["p_look1"]), ("p1_2", g["p_look2"])]
                    if hier == "stagewise":
                        gq += [("p2_1", g["a_look1"]), ("p2_2", g["a_look2"])]
                    for name, val in gq:
                        rows.append([c["config"], mtp, hier, hyp, repr(float(th1)), repr(float(th2)), "glimm",
                                     name, repr(float(val))])
                sim = ref_mtp.simulate(cfg, 2_000_000, 20261020)
                for (name, val), (_, sv) in zip(quantities(res, k1, k2), quantities(sim, k1, k2)):
                    se = sim["se"](sv)
                    mc_rows.append([c["config"], mtp, hier, hyp, name, repr(float(val)), repr(float(sv)),
                                    repr(float(se)), repr(float((sv - val) / se)) if se > 0 else ""])
                print("%-12s %-15s %-10s %-12s %.1fs" % (c["config"], mtp, hier, hyp, time.time() - t0), flush=True)
    with open(os.path.join(OUT, "gs_prob_mtp.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["config", "mtp", "hierarchical", "hypothesis", "theta1", "theta2", "method", "quantity", "value"])
        w.writerows(rows)
    with open(os.path.join(OUT, "gs_prob_mtp_mc.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["config", "mtp", "hierarchical", "hypothesis", "quantity", "recursion", "simulation", "se", "z"])
        w.writerows(mc_rows)
    z = np.array([float(r[-1]) for r in mc_rows if r[-1] != ""])
    print("written: gs_prob_mtp_configs.csv, gs_prob_mtp.csv (%d rows), gs_prob_mtp_mc.csv (%d rows, max |z| %.2f)"
          % (len(rows), len(mc_rows), np.abs(z).max()))


if __name__ == "__main__":
    main()
