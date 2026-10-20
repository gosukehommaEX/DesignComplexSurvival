"""Covariance estimator of Danzer et al. (2025) applied to large simulated trials under the global null hypothesis.

Written independently of the R code. For scenarios V3 and V6 of simulate_two.py under the global null hypothesis (the
control model in both groups), 10 trials of n = 1,000,000 patients are simulated with simulate_ms.simulate_group().
At the calendar looks of the scenario the covariances of the log-rank statistics of PFS and OS are estimated as in
Danzer et al. (2025, Section 2), with mu^{Z=g}(t, s) = 1 - Y^{Z=g}(t, s) / Y(t, s) and
psi^{Z=g}(t, s) = int_0^s mu^{Z=g}(t, u) Lambda-hat(t, du) (Lambda-hat the Nelson-Aalen estimate):
  variance of U_E(t):            (1/n) sum_i mu^{Z=0} mu^{Z=1} Delta_E,i  (at X_E,i(t))
  covariance of U_PFS(t1), U_OS(t2):
      (1/n) sum_i (mu^{Z_i}_PFS Delta_PFS,i - psi^{Z_i}_PFS) (mu^{Z_i}_OS Delta_OS,i - psi^{Z_i}_OS)
  with the PFS quantities at calendar time t1 and the OS quantities at t2, each at the patient's observed time.
The estimated correlation of PFS at look j and OS at look k (quantity "corr_po", look "j-k") is written to
results/danzer_check.csv (columns of sim_python.csv) as the mean over the 10 trials, with its standard error.

Usage (about 2 minutes):  python3 danzer_check.py
"""
import argparse, csv, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import simulate_ms as sm  # noqa: E402
import simulate_two as st  # noqa: E402

SEEDS = {"V3": 33, "V6": 63}


def functionals(entry, tt, drop, trt, t):
    """Per-patient factor mu Delta - psi and the variance term mu0 mu1 Delta of one endpoint at calendar time t."""
    n = entry.size
    enr = entry <= t
    obs = enr & (tt <= drop) & (entry + tt <= t)
    x = np.where(obs, tt, np.minimum(np.minimum(tt, drop), np.maximum(t - entry, 0.0)))
    idx = np.where(enr)[0]
    o = idx[np.argsort(x[idx], kind="stable")]
    m = o.size
    y = (m - np.arange(m)).astype(float)                          # at risk at the sorted times (no ties)
    g = trt[o].astype(float)
    y1 = np.cumsum(g[::-1])[::-1]                                 # treatment patients at risk
    mu1 = 1 - y1 / y                                              # mu^{Z=1} = Y^{Z=0} / Y
    mu0 = y1 / y                                                  # mu^{Z=0} = Y^{Z=1} / Y
    d = obs[o].astype(float)
    psi1 = np.cumsum(d * mu1 / y)                                 # int_0^X mu^{Z=1} dLambda-hat
    psi0 = np.cumsum(d * mu0 / y)
    h = np.zeros(n); v = np.zeros(n)
    h[o] = np.where(g == 1, mu1 * d - psi1, mu0 * d - psi0)
    v[o] = mu0 * mu1 * d
    return h, v


def run(name, reps=10, n=1_000_000, seed=20261022):
    sc = st.scenario(name)
    cal = sc["cal"]; nc = len(cal)
    brk, mass = sm.accrual_setup(sc["accrual"], n)
    p1 = sc["alloc"] / (1 + sc["alloc"])
    ss = np.random.SeedSequence([seed, SEEDS[name]])
    est = np.zeros((reps, nc, nc))
    for r, child in enumerate(ss.spawn(reps)):
        rng = np.random.default_rng(child)
        entry = sm.draw_entry(rng, brk, mass, n)
        trt = rng.random(n) < p1
        lat = {"pfs": np.empty(n), "os": np.empty(n)}
        drop = np.empty(n)
        for flag in (False, True):
            ix = np.where(trt == flag)[0]
            lat["pfs"][ix], lat["os"][ix], _, drop[ix] = sm.simulate_group(rng, sc["control"], ix.size)
        f = {(ep, j): functionals(entry, lat[ep], drop, trt, c) for ep in ("pfs", "os") for j, c in enumerate(cal)}
        for j in range(nc):
            for k in range(nc):
                hp, vp = f[("pfs", j)]
                ho, vo = f[("os", k)]
                est[r, j, k] = np.mean(hp * ho) / np.sqrt(np.mean(vp) * np.mean(vo))
        print("  %s trial %d/%d" % (name, r + 1, reps), flush=True)
    rows = []
    for j in range(nc):
        for k in range(nc):
            x = est[:, j, k]
            rows.append(dict(source="danzer", scenario=name, hypothesis="null", look_type="calendar",
                             look="%d-%d" % (j + 1, k + 1), time="", target_events="", quantity="corr_po",
                             endpoint="pfs-os", group="", mean=repr(float(x.mean())), sd=repr(float(x.std(ddof=1))),
                             se=repr(float(x.std(ddof=1) / np.sqrt(reps))), n=reps))
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenarios", default="V3,V6")
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--n", type=int, default=1_000_000)
    ap.add_argument("--out", default=os.path.join(HERE, "..", "results", "danzer_check.csv"))
    a = ap.parse_args()
    rows = []
    for nm in a.scenarios.split(","):
        rows += run(nm, a.reps, a.n)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=st.COLS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("written", a.out, len(rows), "rows")
