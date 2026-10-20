"""Expected values for the tests of ms_power() and ms_size() with two endpoints or with looks driven by the events of
different endpoints (stage 2, step 3).

Run from this folder with:  python3 expected_design.py
Requires numpy and scipy. Results are written to results/design2.csv (exact values) and results/design2_mc.csv
(Monte Carlo check of the rejection probabilities) and read by tests/testthat (helper-expected.R).

The designs use scenario F of expected_zeta.py (Fleischer-type model, n = 400) and are defined identically in
tests/testthat/helper-designs.R:
  D1  PFS, then OS (fixed sequence, overall); looks driven by PFS events at 60 and 100 per cent of 250 and by 220
      OS events; PFS tested at looks 1 and 2, OS at looks 2 and 3; O'Brien-Fleming type spending for both;
      variance "studentized"; null scenario with the same OS in both groups ("os")
  D2  PFS and OS (Bonferroni, weights 0.8 and 0.2); look 1 at month 15, look 2 driven by OS events 12 months after
      the end of accrual; both endpoints at both looks; O'Brien-Fleming type spending for PFS and Pocock type for
      OS; non-binding futility for PFS by beta spending (HSD, -2, beta 0.1); variance "null"; null scenario with
      the same PFS in both groups ("pfs")
  D3  PFS, then OS (fixed sequence, stagewise); looks driven by OS events at 50 and 100 per cent of 220; both
      endpoints at both looks; O'Brien-Fleming type spending; variance "alternative"
  D4  OS only; look 1 driven by 250 PFS events, look 2 by 220 OS events; O'Brien-Fleming type spending;
      variance "studentized"

Each quantity is computed independently of the R code:
  analysis times      root finding on ref_logrank.events() (adaptive quadrature)
  events, Z means     ref_logrank.moments()
  covariance          ref_zeta.zeta_cov() with the moments of ref_zeta.same_moments() and
                      ref_zeta.cross_moments_fleischer(), then the choice of 'variance' (written here)
  bounds              ref_gsd.bounds() (recursive conditioning with adaptive quadrature)
  rejection probs     the procedure written here as a sum over the disjoint sample-path outcomes of rectangle
                      probabilities computed by ref_gsd.mvn_box() (recursive conditioning with adaptive
                      quadrature); for one endpoint ref_gsd.stop_probs()
  Monte Carlo         ref_mtp.simulate() with the mean vector and covariance matrix above (4 million trials)
"""
import csv
import os
import sys
import time

import numpy as np
from scipy import optimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ref_gsd as rg  # noqa: E402
import ref_logrank as rl  # noqa: E402
import ref_model as rm  # noqa: E402
import ref_mtp as rmtp  # noqa: E402
import ref_zeta as rz  # noqa: E402
from expected_zeta import F_RATES, scenario_f, write  # noqa: E402

LN2 = np.log(2.0)
N = 400
MC_TRIALS = 4_000_000


def fleischer_groups(rates):
    """(a, b, c) of each group for ref_zeta.cross_moments_fleischer."""
    return [rz.FleischerGroup(a, b, c) for a, b, c in rates]


def scen(rates):
    (a0, b0, c0), (a1, b1, c1) = rates
    return rl.Scenario(rm.Model(b0, prog=a0, death_prog=b0, dropout=c0),
                       rm.Model(b1, prog=a1, death_prog=b1, dropout=c1), 2.0,
                       rl.Accrual(time=[0, 4], rate=[10, 25]))


B0, B1 = F_RATES["b"]
RATES = {
    "F": ((LN2 / 3 - B0, B0, 0.004), (LN2 / 4.3 - B1, B1, 0.006)),
    # same OS in both groups (death hazard of the control group), PFS as in F
    "F_os_null": ((LN2 / 3 - B0, B0, 0.004), (LN2 / 4.3 - B0, B0, 0.006)),
    # same PFS in both groups (median 3), OS as in F
    "F_pfs_null": ((LN2 / 3 - B0, B0, 0.004), (LN2 / 3 - B1, B1, 0.006)),
    "F_global_null": ((LN2 / 3 - B0, B0, 0.004), (LN2 / 3 - B0, B0, 0.004)),
}

DESIGNS = {
    "D1": dict(endpoints=["pfs", "os"], mtp="fixed_sequence", hierarchical="overall", weights=None,
               sf=["OF", "OF"], param=[None, None], futility="none", variance="studentized",
               trigger=["pfs", "pfs", "os"], fraction=[0.6, 1.0, 1.0], time=[None, None, None],
               test=[["pfs"], ["pfs", "os"], ["os"]], events={"pfs": 250.0, "os": 220.0}, followup=None,
               null=[("os", "F_os_null")]),
    "D2": dict(endpoints=["pfs", "os"], mtp="bonferroni", hierarchical="overall", weights=[0.8, 0.2],
               sf=["OF", "Pocock"], param=[None, None], futility="beta", sf_fut="HSD", param_fut=-2.0,
               beta=0.1, variance="null", trigger=["calendar", "os"], fraction=[None, 1.0], time=[15.0, None],
               test=None, events=None, followup=12.0, null=[("pfs", "F_pfs_null")]),
    "D3": dict(endpoints=["pfs", "os"], mtp="fixed_sequence", hierarchical="stagewise", weights=None,
               sf=["OF", "OF"], param=[None, None], futility="none", variance="alternative",
               trigger=["os", "os"], fraction=[0.5, 1.0], time=[None, None], test=None,
               events={"os": 220.0}, followup=None, null=[]),
    "D4": dict(endpoints=["os"], sf=["OF"], param=[None], futility="none", variance="studentized",
               trigger=["pfs", "os"], fraction=[1.0, 1.0], time=[None, None], test=None,
               events={"pfs": 250.0, "os": 220.0}, followup=None, null=[]),
}


# ---------------------------------------------------------------------------
# Analysis times
# ---------------------------------------------------------------------------
def n_events(scn, t, ep):
    return rl.events(scn, N, t, ep)["events"]


def look_times(scn, d, targets=None):
    """Times of the looks and event targets (event-driven looks)."""
    k = len(d["trigger"])
    eps = [e for e in ("pfs", "os") if e in d["trigger"]]
    A = scn.accrual.setup(N)[2]
    if targets is None:
        if d["followup"] is not None:
            t_ref = A + d["followup"]
            final = {e: n_events(scn, t_ref, e) for e in eps}
        elif d["events"] is not None:
            t_ref = None
            final = dict(d["events"])
        else:
            t_ref = d["time"][-1]
            final = {e: n_events(scn, t_ref, e) for e in eps}
        targets = [final[d["trigger"][j]] * d["fraction"][j] if d["trigger"][j] != "calendar" else None
                   for j in range(k)]
    else:
        final, t_ref = None, None
    times = []
    for j in range(k):
        if d["trigger"][j] == "calendar":
            times.append(d["time"][j])
        elif j == k - 1 and d["followup"] is not None and t_ref is not None:
            times.append(t_ref)
        else:
            times.append(rl.event_time(scn, N, targets[j], d["trigger"][j]))
    assert all(np.diff(times) > 0), times
    return np.array(times), targets, final


def look_trigger(d, j, fixed_final):
    """Endpoint whose observed events determine the time of look j (None for calendar and fixed times)."""
    if d["trigger"][j] == "calendar" or (fixed_final and j == len(d["trigger"]) - 1):
        return None
    return d["trigger"][j]


# ---------------------------------------------------------------------------
# Joint distribution
# ---------------------------------------------------------------------------
def joint(rates, d, times, tested, fixed_final, variance, mean_zero=False):
    scn = scen(rates)
    alpha, rho, _ = scn.accrual.setup(N)
    items = [(e, times[j], look_trigger(d, j, fixed_final)) for e in d["endpoints"] for j in tested[e]]
    tset = {e: sorted(set([t for (ee, t, _) in items if ee == e] + [t for (_, t, tr) in items if tr == e]))
            for e in ("pfs", "os")}
    tset = {e: v for e, v in tset.items() if v}
    mom_same = {e: (tset[e],) + rz.same_moments(scn, N, tset[e], e) for e in tset}
    mom_cross = {}
    if len(tset) == 2:
        groups = fleischer_groups(rates)
        for tp in tset["pfs"]:
            for to in tset["os"]:
                mom_cross[(tp, to)] = rz.cross_moments_fleischer(groups, (scn.p0, scn.p1), alpha, rho, tp, to)
    core = {(e, t): rz.per_patient_core(scn, N, t, e) for e in tset for t in tset[e]}
    S = rz.zeta_cov(items, mom_same, mom_cross, core)
    v = np.array([core[(e, t)]["v"] for (e, t, _) in items])
    same = np.equal.outer([e for (e, _, _) in items], [e for (e, _, _) in items])
    if variance == "studentized":
        sig = S
    elif variance == "null":
        sd = np.sqrt(np.diag(S))
        sig = S / np.outer(sd, sd)
        vv = np.sqrt(np.minimum.outer(v, v) / np.maximum.outer(v, v))
        sig[same] = vv[same]
    else:
        def cov_psi(x, y):
            (e1, t1, _), (e2, t2, _) = x, y
            if e1 == e2:
                ts, mean, prod = mom_same[e1]
                i, j = ts.index(t1), ts.index(t2)
                return prod[i, j, 0, 0] - mean[i, 0] * mean[j, 0]
            if e1 == "os":
                return cov_psi(y, x)
            iP, iO = mom_same["pfs"][0].index(t1), mom_same["os"][0].index(t2)
            return mom_cross[(t1, t2)][0, 0] - mom_same["pfs"][1][iP, 0] * mom_same["os"][1][iO, 0]
        sig = np.array([[cov_psi(x, y) for y in items] for x in items]) / np.sqrt(np.outer(v, v))
    mom = {e: [rl.moments(scn, N, times[j], e) for j in tested[e]] for e in d["endpoints"]}
    mean = np.zeros(len(items)) if mean_zero else np.concatenate(
        [[m["z_direct"] for m in mom[e]] for e in d["endpoints"]])
    return mean, (sig + sig.T) / 2, mom


# ---------------------------------------------------------------------------
# Rejection probabilities of the procedure
# ---------------------------------------------------------------------------
def procedure_probs(mean, sig, an1, an2, u1, u2, l1, mtp, hier):
    """Rejection probabilities by a sum over the disjoint outcomes of the sample path.

    Outcomes of the first hypothesis: rejection at look j, stop for futility at interim look j (finite bound),
    or no rejection. Given the outcome, the second hypothesis is tested at its looks taking place while the trial
    runs and, for the fixed sequence, after (overall) or at (stagewise) the rejection of the first; its outcome is
    the first of these looks at which its statistic reaches the bound.
    """
    k1, k2 = len(u1), len(u2)
    l1 = np.full(k1, -np.inf) if l1 is None else np.array(l1, float)
    l1[-1] = -np.inf

    def box(idx, lo, hi):
        lo, hi = np.asarray(lo, float), np.asarray(hi, float)
        keep = np.isfinite(lo) | np.isfinite(hi)
        i = np.asarray(idx)[keep]
        if i.size == 0:
            return 1.0
        return rg.mvn_box(lo[keep], hi[keep], mean[i], sig[np.ix_(i, i)])

    p1 = np.zeros(k1); pf = np.zeros(k1); p2 = np.zeros(k2); both = 0.0
    outcomes = [("reject", j) for j in range(k1)] + [("futility", j) for j in range(k1 - 1)
                                                     if np.isfinite(l1[j])] + [("none", k1 - 1)]
    for kind, j in outcomes:
        stop_an = an1[j] if kind == "futility" else None
        if mtp == "bonferroni":
            elig = [i for i in range(k2) if stop_an is None or an2[i] <= stop_an]
        elif kind != "reject":
            elig = []
        elif hier == "overall":
            elig = [i for i in range(k2) if an2[i] >= an1[j]]
        else:
            elig = [i for i in range(k2) if an2[i] == an1[j]]
        lo1 = [l1[m] for m in range(j)] + [u1[j] if kind == "reject" else -np.inf]
        hi1 = [u1[m] for m in range(j)] + [{"reject": np.inf, "futility": l1[j], "none": u1[j]}[kind]]
        idx1 = list(range(j + 1))
        pr = box(idx1, lo1, hi1)
        if kind == "reject":
            p1[j] = pr
        elif kind == "futility":
            pf[j] = pr
        for r, i in enumerate(elig):
            q = box(idx1 + [k1 + m for m in elig[:r + 1]],
                    lo1 + [-np.inf] * r + [u2[i]], hi1 + [u2[m] for m in elig[:r]] + [np.inf])
            p2[i] += q
            if kind == "reject":
                both += q
    return dict(p1=p1, pf=pf, p2=p2, first=p1.sum(), second=p2.sum(), both=both,
                any=p1.sum() + p2.sum() - both)


def bounds_of(d, mom, j_ep, sig_block, beta=None):
    e = d["endpoints"][j_ep]
    dd = np.array([m["events"] for m in mom[e]])
    tfrac = dd / dd[-1]
    corr0 = np.sqrt(np.minimum.outer(dd, dd) / np.maximum.outer(dd, dd))
    alpha = 0.025 * (d["weights"][j_ep] if d.get("mtp") == "bonferroni" else 1.0)
    if j_ep == 0 and d["futility"] == "beta":
        mu = np.array([m["z_direct"] for m in mom[e]])
        return rg.bounds(tfrac, alpha=alpha, sf_eff=d["sf"][j_ep], futility="beta", sf_fut=d["sf_fut"],
                         param_fut=d["param_fut"], beta=d["beta"], mean=mu, corr=corr0, sigma_alt=sig_block)
    return rg.bounds(tfrac, alpha=alpha, sf_eff=d["sf"][j_ep], corr=corr0)


def run(name):
    d = DESIGNS[name]
    t0 = time.time()
    rows = []
    mc_rows = []

    def add(q, ep="", look="", value=None, row="", col=""):
        rows.append(dict(config=name, quantity=q, endpoint=ep, look=look, row=row, col=col, value=float(value)))
    eps = d["endpoints"]
    k = len(d["trigger"])
    tested = {e: [j for j in range(k) if d["test"] is None or e in d["test"][j]] for e in eps}
    an = [[j + 1 for j in tested[e]] for e in eps]
    fixed_final = d["followup"] is not None
    times, targets, final = look_times(scen(RATES["F"]), d)
    for j in range(k):
        add("time", look=j + 1, value=times[j])
    for e, v in final.items():
        add("final_events", ep=e, value=v)
    mean, sig, mom = joint(RATES["F"], d, times, tested, fixed_final, d["variance"])
    for i in range(len(mean)):
        for j in range(len(mean)):
            add("sigma", row=i + 1, col=j + 1, value=sig[i, j])
    off = 0
    bnd = []
    for je, e in enumerate(eps):
        kj = len(tested[e])
        blk = sig[off:off + kj, off:off + kj]
        lo, up = bounds_of(d, mom, je, blk)
        bnd.append((lo, up))
        for m in range(kj):
            look = tested[e][m] + 1
            add("events", e, look, mom[e][m]["events"])
            add("z_mean", e, look, mom[e][m]["z_direct"])
            add("z_sd", e, look, np.sqrt(blk[m, m]))
            add("z_upper", e, look, up[m])
            add("z_lower", e, look, lo[m])
        off += kj
    if len(eps) == 1:
        lo, up = bnd[0]
        pu, pl = rg.stop_probs(lo, up, mean, sig)
        for m in range(len(up)):
            add("prob_upper", eps[0], tested[eps[0]][m] + 1, pu[m])
            add("prob_lower", eps[0], tested[eps[0]][m] + 1, pl[m])
        add("power", eps[0], value=pu.sum())
        s0 = scen(RATES["F_global_null"])
        times0, _, _ = look_times(s0, d, targets)
        for j in range(k):
            add("time0", look=j + 1, value=times0[j])
        v0 = np.array([rl.moments(s0, N, t, eps[0])["info"] for t in times0])
        c0 = np.sqrt(np.minimum.outer(v0, v0) / np.maximum.outer(v0, v0))
        pu0, _ = rg.stop_probs(np.full(len(up), -np.inf), up, np.zeros(len(up)), c0)
        pu0b, _ = rg.stop_probs(lo, up, np.zeros(len(up)), c0)
        add("type1_nonbinding", value=pu0.sum())
        add("type1_binding", value=pu0b.sum())
        print(name, "%.0fs" % (time.time() - t0), flush=True)
        return rows, mc_rows
    (l1, u1), (_, u2) = bnd
    lfut = l1 if d["futility"] != "none" else None
    pr = procedure_probs(mean, sig, an[0], an[1], u1, u2, lfut, d["mtp"], d["hierarchical"])
    for m, j in enumerate(tested[eps[0]]):
        add("prob_upper", eps[0], j + 1, pr["p1"][m])
        add("prob_lower", eps[0], j + 1, pr["pf"][m])
    for m, j in enumerate(tested[eps[1]]):
        add("prob_upper", eps[1], j + 1, pr["p2"][m])
    add("power", eps[0], value=pr["first"])
    add("power", eps[1], value=pr["second"])
    add("power", "both", value=pr["both"])
    add("power", "any", value=pr["any"])
    cfg = dict(a1=an[0], a2=an[1], u1=u1, u2=u2, l1=lfut, mtp=d["mtp"], hierarchical=d["hierarchical"])
    sim = rmtp.simulate(cfg, MC_TRIALS, 20261020, mean=mean, sigma=sig)
    for q in ("first", "second", "both", "any"):
        mc_rows.append(dict(config=name, scenario="design", quantity=q, exact=pr[q], mc=sim[q], se=sim["se"](sim[q]),
                            trials=MC_TRIALS))
    # familywise error rate under the global null scenario and the partial null scenarios
    for label, key in [("global", "F_global_null")] + list(d["null"]):
        rates = RATES[key] if label != "global" else RATES["F_global_null"]
        sq = scen(rates)
        tq, _, _ = look_times(sq, d, targets)
        if label == "global":
            for j in range(k):
                add("time0", look=j + 1, value=tq[j])
        mq, sq_sig, _ = joint(rates, d, tq, tested, False, "null" if label == "global" else d["variance"],
                              mean_zero=(label == "global"))
        nb = procedure_probs(mq, sq_sig, an[0], an[1], u1, u2, None, d["mtp"], d["hierarchical"])
        bd = nb if lfut is None else procedure_probs(mq, sq_sig, an[0], an[1], u1, u2, lfut, d["mtp"],
                                                     d["hierarchical"])
        qk = "any" if label == "global" else ("first" if label == eps[0] else "second")
        add("fwer_nonbinding", label, value=nb[qk])
        add("fwer_binding", label, value=bd[qk])
        if label == "global":
            for m, j in enumerate(tested[eps[0]]):
                add("cum_alpha", eps[0], j + 1, np.cumsum(nb["p1"])[m])
            for m, j in enumerate(tested[eps[1]]):
                add("cum_alpha", eps[1], j + 1, np.cumsum(nb["p2"])[m])
            sim0 = rmtp.simulate(dict(cfg, l1=None), MC_TRIALS, 20261021, mean=mq, sigma=sq_sig)
            mc_rows.append(dict(config=name, scenario="global", quantity="any", exact=nb["any"], mc=sim0["any"],
                                se=sim0["se"](sim0["any"]), trials=MC_TRIALS))
    print(name, "%.0fs" % (time.time() - t0), flush=True)
    return rows, mc_rows


def main():
    names = sys.argv[1:] or list(DESIGNS)
    rows, mc = [], []
    for name in names:
        r, m = run(name)
        rows += r
        mc += m
    write("design2.csv", rows)
    write("design2_mc.csv", mc)


if __name__ == "__main__":
    main()
