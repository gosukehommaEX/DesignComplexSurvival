"""Type I error of the validation scenarios V2 and V7 at the looks driven by OS events, by larger simulations under
the null hypothesis than those of validate_fastsurvival.R and simulate_ms.py (see ../README.md, "Type I error at the
looks driven by OS events").

The patients are generated with the functions of simulate_ms.py (simulate_group(), draw_entry(), logrank_z(),
kth_time()); the analysis at the d-th OS death and the log-rank Z (positive when the treatment group has fewer deaths
than expected) are as in simulate_ms.run(). Bounds: ref_gsd.bounds() for the information fractions d / d_K. New
seeds (SeedSequence with 20261029).

Variants
  V7          scenario C, forward clock (simulate_ms.scenario_c("forward")), control model in both groups, simple
              randomization 1:1, n = 500, looks at 150 and 250 OS deaths; as V7, with 10 times the trials
  B_1:2       scenario B (helper-scenarios.R scn_b()), control model (death 0.05 then 0.04 from month 6, dropout
              0.01) in both groups, n = 300 with fixed group sizes 100 (control) and 200 (treatment) and fixed numbers
              of patients in the accrual intervals (30 in months 0 to 3, 270 in months 3 to 12), as V2 in
              FastSurvival; looks at 75 and 150 OS deaths
  B_2:1       the same with 200 (control) and 100 (treatment)
  B_1:1       the same with 150 and 150
For each variant: the rejection rate in the direction of benefit (Z_k >= b_k at a look), the rate in the opposite
direction (Z_k <= -b_k), the rate at look 1, and the mean, standard deviation and skewness of Z at each look.

Usage (numpy and scipy; about 11 minutes):
  python3 inst/validation/python/simulate_type1.py [--reps-v7 200000] [--reps-b 400000] [--out <csv>]
Writes ../results/sim_python_type1.csv by default (columns variant, quantity, value, se, n).
"""
import argparse
import csv
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ref_gsd  # noqa: E402
import simulate_ms as sm  # noqa: E402

BASE_SEED = 20261029
CHUNK_PATIENTS = 2_000_000


def z_at_event_looks(rng, nr, n, entry, trt, ctl_model, events):
    """Log-rank Z of OS at the looks driven by OS deaths (as simulate_ms.run(), control model in both groups)."""
    rep = np.repeat(np.arange(nr), n)
    pfs = np.empty(nr * n); os_ = np.empty(nr * n); sw = np.empty(nr * n); drop = np.empty(nr * n)
    for flag in (False, True):
        ix = np.where(trt == flag)[0]
        pfs[ix], os_[ix], sw[ix], drop[ix] = sm.simulate_group(rng, ctl_model, ix.size)
    death_cal = entry + os_
    seen = os_ <= drop
    out = []
    for d in events:
        cut = sm.kth_time(rep[seen], death_cal[seen], d, nr)
        cut = np.where(np.isnan(cut), np.inf, cut)
        f = cut[rep] - entry
        enr = f >= 0
        obs_ev = enr & (os_ <= drop) & (entry + os_ <= cut[rep])
        ti = np.where(obs_ev, os_, np.minimum(np.minimum(os_, drop), np.maximum(f, 0)))
        out.append(sm.logrank_z(rep[enr], ti[enr], obs_ev[enr], trt[enr], nr))
    return out


def run_v7(reps):
    sc = sm.scenario_c("forward")
    n = sc["n"]
    brk, mass = sm.accrual_setup(sc["accrual"], n)
    per = max(1, CHUNK_PATIENTS // n)
    ss = np.random.SeedSequence([BASE_SEED, 7])
    zs = [[], []]
    for ci, child in enumerate(ss.spawn(int(np.ceil(reps / per)))):
        rng = np.random.default_rng(child)
        nr = min(per, reps - ci * per)
        entry = sm.draw_entry(rng, brk, mass, nr * n)
        trt = rng.random(nr * n) < 0.5
        z = z_at_event_looks(rng, nr, n, entry, trt, sc["control"], sc["events"])
        for k in range(2):
            zs[k].append(z[k])
    return [np.concatenate(z) for z in zs], sc["events"]


def run_b(reps, n_ctl, n_trt, tag):
    ctl = sm.Group([0.05, 0.04], time=[0, 6], dropout=0.01)
    n = n_ctl + n_trt
    events = [75, 150]
    per = max(1, CHUNK_PATIENTS // n)
    ss = np.random.SeedSequence([BASE_SEED, 2, tag])
    zs = [[], []]
    for ci, child in enumerate(ss.spawn(int(np.ceil(reps / per)))):
        rng = np.random.default_rng(child)
        nr = min(per, reps - ci * per)
        # fixed numbers of patients in the accrual intervals and fixed group sizes (FastSurvival)
        e1 = rng.random((nr, 30)) * 3.0
        e2 = 3.0 + rng.random((nr, 270)) * 9.0
        entry = np.concatenate([e1, e2], axis=1).ravel()
        lab = np.concatenate([np.zeros(n_ctl, bool), np.ones(n_trt, bool)])
        trt = np.array([rng.permutation(lab) for _ in range(nr)]).ravel()
        z = z_at_event_looks(rng, nr, n, entry, trt, ctl, events)
        for k in range(2):
            zs[k].append(z[k])
    return [np.concatenate(z) for z in zs], events


def summarize(name, z, events, rows):
    frac = np.asarray(events, float) / events[-1]
    _, b = ref_gsd.bounds(frac, alpha=0.025, sf_eff="OF")
    z1, z2 = z
    ok = np.isfinite(z1) & np.isfinite(z2)
    z1, z2 = z1[ok], z2[ok]
    m = z1.size

    def prop(x):
        p = x.mean()
        return p, np.sqrt(p * (1 - p) / m)
    up1 = z1 >= b[0]; upa = up1 | (z2 >= b[1])
    lo1 = z1 <= -b[0]; loa = lo1 | (z2 <= -b[1])
    for q, x in (("reject_upper_look1", up1), ("reject_upper_any", upa), ("reject_lower_look1", lo1),
                 ("reject_lower_any", loa)):
        p, se = prop(x)
        rows.append((name, q, p, se, m))
    for k, zk in enumerate((z1, z2), 1):
        mu, sd = zk.mean(), zk.std(ddof=1)
        sk = np.mean(((zk - mu) / sd) ** 3)
        rows.append((name, "z_mean_look%d" % k, mu, sd / np.sqrt(m), m))
        rows.append((name, "z_sd_look%d" % k, sd, sd / np.sqrt(2 * (m - 1)), m))
        rows.append((name, "z_skewness_look%d" % k, sk, np.sqrt(6.0 / m), m))
    r = np.corrcoef(z1, z2)[0, 1]
    rows.append((name, "z_corr", r, (1 - r ** 2) / np.sqrt(m), m))
    rows.append((name, "bound_look1", b[0], 0.0, 0))
    rows.append((name, "bound_look2", b[1], 0.0, 0))
    rows.append((name, "not_reached", float((~ok).sum()), 0.0, ok.size))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps-v7", type=int, default=200000)
    ap.add_argument("--reps-b", type=int, default=400000)
    ap.add_argument("--out", default=os.path.join(HERE, "..", "results", "sim_python_type1.csv"))
    a = ap.parse_args()
    rows = []
    t0 = time.time()
    for name, n_ctl, n_trt, tag in (("B_1:2", 100, 200, 12), ("B_2:1", 200, 100, 21), ("B_1:1", 150, 150, 11)):
        z, ev = run_b(a.reps_b, n_ctl, n_trt, tag)
        summarize(name, z, ev, rows)
        print("%s done (%.0f s)" % (name, time.time() - t0), flush=True)
    z, ev = run_v7(a.reps_v7)
    summarize("V7", z, ev, rows)
    print("V7 done (%.0f s)" % (time.time() - t0), flush=True)
    path = os.path.abspath(a.out)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="\n") as f:
        f.write("variant,quantity,value,se,n\n")
        for r in rows:
            f.write("%s,%s,%s,%s,%d\n" % (r[0], r[1], repr(float(r[2])), repr(float(r[3])), r[4]))
    for r in rows:
        print("%-6s %-20s %12.6f  se %.6f  n %d" % r)
    print("written %s; total %.0f s" % (path, time.time() - t0))


if __name__ == "__main__":
    main()
