"""Patient-level simulation of the four-state model, written independently of the R code.

Scenarios V6-V8 of the Phase 3 plan (section 4.2), which FastSurvival cannot simulate: a hazard of switching after
progression, the forward clock, and treatment discontinuation (dilution) of Luo et al. (2019), Example 1.

Each patient is generated directly: entry time from the piecewise uniform accrual, group by simple randomization,
the exit from the progression-free state by inversion of the piecewise linear cumulative hazard (the cause drawn in
proportion to the hazards), the switch at progression with the probability at the progression time, the exit from the
progressed state (death or switch) and the death after a switch in the same way on the post-progression clock (time
since progression for clock = "reset", time since randomization for clock = "forward"), and the dropout time.

Only summaries are written (no patient-level data), to results/sim_python.csv:
  calendar looks: group-wise OS and PFS events and dropouts, log-rank Z of OS and PFS (positive when the treatment
    is better), and the correlation of Z between the looks;
  event looks (calendar time of the d-th OS death in each trial): look time, group-wise events, Z, correlation,
    rejection of the O'Brien-Fleming design (bounds from ref_gsd.py for the information fractions d / d_K);
  latent quantities (alternative hypothesis only, no dropout): state occupation probabilities and OS and PFS
    survival at times since randomization.

Usage:
  python3 simulate_ms.py                    # all scenarios, 20000 trials under each hypothesis
  python3 simulate_ms.py --scenarios V6 --reps 2000 --out /tmp/x.csv
"""
import argparse, csv, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ref_gsd  # noqa: E402


# Model -----------------------------------------------------------------------
def _expand(x, n):
    x = np.atleast_1d(np.asarray(x, float))
    return np.repeat(x, n) if x.size == 1 else x


class Group:
    """One treatment group, with the arguments of ms_model()."""

    def __init__(self, death, time=(0.0,), prog=0.0, death_prog=0.0, switch_prob=0.0, switch_rate=0.0,
                 death_switch=None, dropout=0.0, time_post=(0.0,), clock="reset"):
        self.time = np.asarray(time, float); k = self.time.size
        self.prog = _expand(prog, k); self.death = _expand(death, k)
        self.switch_prob = _expand(switch_prob, k); self.dropout = _expand(dropout, k)
        self.time_post = np.asarray(time_post, float); m = self.time_post.size
        self.death_prog = _expand(death_prog, m); self.switch_rate = _expand(switch_rate, m)
        self.death_switch = self.death_prog.copy() if death_switch is None else _expand(death_switch, m)
        self.clock = clock


def step(x, brk, val):
    """Right-continuous step function with pieces starting at brk."""
    return val[np.maximum(np.searchsorted(brk, x, side="right") - 1, 0)]


def cumhaz(x, brk, rate):
    x = np.maximum(x, 0.0)
    hb = np.concatenate([[0.0], np.cumsum(rate[:-1] * np.diff(brk))])
    i = np.maximum(np.searchsorted(brk, x, side="right") - 1, 0)
    return hb[i] + rate[i] * (x - brk[i])


def inv_cumhaz(c, brk, rate):
    """Smallest x with cumulative hazard c (infinite when the hazard is zero from there on)."""
    hb = np.concatenate([[0.0], np.cumsum(rate[:-1] * np.diff(brk))])
    i = np.searchsorted(hb, c, side="right") - 1
    i = np.clip(i, 0, brk.size - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        x = brk[i] + (c - hb[i]) / rate[i]
    return np.where(rate[i] > 0, x, np.inf)


def draw_after(rng, start, brk, rate):
    """Event time on a clock, given survival to clock value start (vector)."""
    e = rng.exponential(size=start.size)
    return inv_cumhaz(cumhaz(start, brk, rate) + e, brk, rate)


def simulate_group(rng, g, n):
    """Latent times since randomization: pfs, os, switch time (inf if none), dropout."""
    tot = g.prog + g.death
    t1 = draw_after(rng, np.zeros(n), g.time, tot)
    with np.errstate(invalid="ignore", divide="ignore"):
        pp = np.where(np.isfinite(t1), step(t1, g.time, g.prog) / step(t1, g.time, tot), 0.0)
    prog = rng.random(n) < pp
    os_ = t1.copy()
    sw = np.full(n, np.inf)
    ip = np.where(prog)[0]
    if ip.size:
        u = t1[ip]
        tau0 = np.zeros(ip.size) if g.clock == "reset" else u.copy()
        sw0 = rng.random(ip.size) < step(u, g.time, g.switch_prob)
        y = g.switch_rate + g.death_prog
        tau2 = draw_after(rng, tau0, g.time_post, y)
        with np.errstate(invalid="ignore", divide="ignore"):
            ps = np.where(np.isfinite(tau2), step(tau2, g.time_post, g.switch_rate) / step(tau2, g.time_post, y), 0.0)
        later = (~sw0) & (rng.random(ip.size) < ps)
        tau_s = np.where(sw0, tau0, np.where(later, tau2, np.inf))
        switched = sw0 | later
        taud = np.where(switched, draw_after(rng, np.where(switched, tau_s, 0.0), g.time_post, g.death_switch), tau2)
        os_[ip] = u + (taud - tau0)
        sw[ip] = np.where(switched, u + (tau_s - tau0), np.inf)
    drop = draw_after(rng, np.zeros(n), g.time, g.dropout)
    return t1, os_, sw, drop


def accrual_setup(acc, n):
    if acc["type"] == "rate":
        tm = np.asarray(acc["time"], float); rate = np.asarray(acc["rate"], float); k = tm.size
        cs = np.concatenate([[0.0], np.cumsum(rate[:-1] * np.diff(tm))]); ce = np.concatenate([cs[1:], [np.inf]])
        i = int(np.where(ce >= n)[0][0]); dur = tm[i] + (n - cs[i]) / rate[i]
        brk = np.concatenate([tm[tm < dur], [dur]]); mass = rate[:brk.size - 1] * np.diff(brk) / n
    else:
        brk = np.concatenate([np.asarray(acc["time"], float), [float(acc["duration"])]])
        mass = np.asarray(acc["prop"], float)
    return brk, mass / mass.sum()


def draw_entry(rng, brk, mass, size):
    j = rng.choice(mass.size, size=size, p=mass)
    return brk[j] + rng.random(size) * (brk[j + 1] - brk[j])


# Log-rank statistic per trial -----------------------------------------------------
def logrank_z(rep, t, ev, trt, nrep):
    """Z = -(O - E) / sqrt(V) of the treatment group for each trial (positive when the treatment is better)."""
    o = np.lexsort((t, rep))
    r = rep[o]; e = ev[o]; g = trt[o].astype(float)
    nt = np.bincount(r, minlength=nrep)
    start = np.concatenate([[0], np.cumsum(nt)[:-1]])
    pos = np.arange(r.size) - start[r]
    y = (nt[r] - pos).astype(float)
    cg = np.concatenate([np.cumsum(g[::-1])[::-1], [0.0]])
    yt = cg[:-1] - cg[(start + nt)[r]]
    m = e
    u = np.bincount(r[m], weights=g[m] - yt[m] / y[m], minlength=nrep)
    v = np.bincount(r[m], weights=yt[m] * (y[m] - yt[m]) / y[m] ** 2, minlength=nrep)
    with np.errstate(invalid="ignore", divide="ignore"):
        return -u / np.sqrt(v)


def kth_time(rep, x, k, nrep):
    """k-th smallest value of x within each trial (nan if fewer than k)."""
    o = np.lexsort((x, rep)); r = rep[o]; xs = x[o]
    nt = np.bincount(r, minlength=nrep); start = np.concatenate([[0], np.cumsum(nt)[:-1]])
    out = np.full(nrep, np.nan)
    ok = nt >= k
    out[ok] = xs[start[ok] + k - 1]
    return out


# Scenarios ---------------------------------------------------------------------
def scenario_c(clock):
    ctl = Group(0.01, time=[0, 4], prog=[0.10, 0.08], switch_prob=0.3, time_post=[0, 3], death_prog=[0.05, 0.06],
                switch_rate=0.04, death_switch=[0.03, 0.035], dropout=0.005, clock=clock)
    trt = Group(0.008, time=[0, 4], prog=[0.06, 0.05], switch_prob=0.0, time_post=[0, 3], death_prog=[0.045, 0.05],
                switch_rate=0.0, dropout=0.005, clock=clock)
    return dict(control=ctl, treatment=trt, alloc=1.0, n=500,
                accrual=dict(type="prop", time=[0.0, 6.0], prop=[0.3, 0.7], duration=18.0),
                cal=[20.0, 30.0], events=[150, 250], latent=[2.0, 5.0, 10.0, 20.0, 30.0])


def scenario_luo1():
    tm = [0, 6, 12, 24]
    lam0 = np.array([0.09, 0.09, 0.07, 0.05]) / 12
    drop = np.array([0, 0, 0.015, 0.015]) / 12
    disc = np.array([0.04, 0.04, 0.03, 0.03]) / 12
    hr = np.array([1, 0.75, 0.75, 0.75])
    ctl = Group(lam0, time=tm, dropout=drop)
    trt = Group(hr * lam0, time=tm, dropout=drop, prog=disc, death_prog=lam0, time_post=tm, clock="forward")
    return dict(control=ctl, treatment=trt, alloc=1.0, n=4600,
                accrual=dict(type="rate", time=list(range(14)), rate=[100, 200, 300, 300] + [400] * 8 + [300, 200]),
                cal=[30.0, 40.0], events=[361, 722], latent=[6.0, 12.0, 24.0, 36.0])


SCENARIOS = {"V6": lambda: scenario_c("reset"), "V7": lambda: scenario_c("forward"), "V8": scenario_luo1}
SEEDS = {"V6": 61, "V7": 71, "V8": 81}


# Simulation -----------------------------------------------------------------------
def run(name, reps, chunk_patients=2_000_000, base_seed=20261019):
    sc = SCENARIOS[name]()
    n = sc["n"]; p1 = sc["alloc"] / (1 + sc["alloc"])
    brk, mass = accrual_setup(sc["accrual"], n)
    frac = np.asarray(sc["events"], float) / sc["events"][-1]
    _, b = ref_gsd.bounds(frac, alpha=0.025, sf_eff="OF")
    rows = []
    for hyp in ("alternative", "null"):
        trt_model = sc["treatment"] if hyp == "alternative" else sc["control"]
        per = max(1, chunk_patients // n)
        ss = np.random.SeedSequence([base_seed, SEEDS[name], 0 if hyp == "alternative" else 1])
        nchunk = int(np.ceil(reps / per))
        acc = {}  # per-trial statistics, concatenated over chunks
        lat = {}  # latent counts
        for ci, child in enumerate(ss.spawn(nchunk)):
            rng = np.random.default_rng(child)
            nr = min(per, reps - ci * per)
            rep = np.repeat(np.arange(nr), n)
            entry = draw_entry(rng, brk, mass, nr * n)
            trt = rng.random(nr * n) < p1
            pfs = np.empty(nr * n); os_ = np.empty(nr * n); sw = np.empty(nr * n); drop = np.empty(nr * n)
            for flag, g in ((False, sc["control"]), (True, trt_model)):
                ix = np.where(trt == flag)[0]
                pfs[ix], os_[ix], sw[ix], drop[ix] = simulate_group(rng, g, ix.size)
            # latent quantities (no dropout), alternative only
            if hyp == "alternative":
                for gname, flag in (("control", False), ("treatment", True)):
                    gi = trt == flag
                    for t in sc["latent"]:
                        s1 = pfs[gi] > t; dead = os_[gi] <= t; swt = (sw[gi] <= t) & ~dead & ~s1
                        cnt = np.array([s1.sum(), (~s1 & ~dead & ~swt).sum(), swt.sum(), dead.sum(),
                                        (os_[gi] > t).sum(), (pfs[gi] > t).sum(), gi.sum()], float)
                        key = (gname, t)
                        lat[key] = lat.get(key, 0) + cnt

            def at_cutoff(cut, tag):
                f = cut[rep] - entry
                enr = f >= 0
                out = {}
                for ep, tt in (("os", os_), ("pfs", pfs)):
                    obs_ev = enr & (tt <= drop) & (tt <= f)
                    obs_do = enr & (drop < tt) & (drop <= f)
                    for gname, flag in (("control", False), ("treatment", True)):
                        out["events_%s_%s" % (ep, gname)] = np.bincount(rep, weights=obs_ev & (trt == flag), minlength=nr)
                        out["dropouts_%s_%s" % (ep, gname)] = np.bincount(rep, weights=obs_do & (trt == flag), minlength=nr)
                    ti = np.minimum(np.minimum(tt, drop), np.maximum(f, 0))
                    out["z_%s" % ep] = logrank_z(rep[enr], ti[enr], obs_ev[enr], trt[enr], nr)
                for k, v in out.items():
                    acc.setdefault((tag, k), []).append(v)

            for li, c in enumerate(sc["cal"], 1):
                at_cutoff(np.full(nr, c), ("calendar", li, c))
            death_cal = entry + os_
            seen = os_ <= drop
            for li, d in enumerate(sc["events"], 1):
                cut = kth_time(rep[seen], death_cal[seen], d, nr)
                acc.setdefault((("event", li, d), "time"), []).append(cut)
                cut_f = np.where(np.isnan(cut), np.inf, cut)
                at_cutoff(cut_f, ("event", li, d))
            print("  %s %s chunk %d/%d (%d trials)" % (name, hyp, ci + 1, nchunk, nr), flush=True)
        st = {k: np.concatenate(v) for k, v in acc.items()}
        nrep_tot = reps

        def add(look_type, look, tval, quantity, endpoint, group, x, target=""):
            x = np.asarray(x, float); ok = np.isfinite(x)
            m = x[ok].mean(); s = x[ok].std(ddof=1) if ok.sum() > 1 else 0.0
            rows.append(dict(source="python", scenario=name, hypothesis=hyp, look_type=look_type, look=look, time=tval,
                             target_events=target, quantity=quantity, endpoint=endpoint, group=group, mean=m, sd=s,
                             se=s / np.sqrt(ok.sum()), n=int(ok.sum())))

        tags = [("calendar", i + 1, c) for i, c in enumerate(sc["cal"])] + [("event", i + 1, d) for i, d in enumerate(sc["events"])]
        for tag in tags:
            lt, li, val = tag
            tval = val if lt == "calendar" else np.nan
            tgt = val if lt == "event" else ""
            if lt == "event":
                add(lt, li, np.nan, "time", "os", "", st[(tag, "time")], tgt)
            for ep in ("os", "pfs"):
                for gname in ("control", "treatment"):
                    add(lt, li, tval, "events", ep, gname, st[(tag, "events_%s_%s" % (ep, gname))], tgt)
                    add(lt, li, tval, "dropouts", ep, gname, st[(tag, "dropouts_%s_%s" % (ep, gname))], tgt)
                add(lt, li, tval, "z", ep, "", st[(tag, "z_%s" % ep)], tgt)
        for lt in ("calendar", "event"):
            tg = [t for t in tags if t[0] == lt]
            for ep in ("os", "pfs"):
                z1 = st[(tg[0], "z_%s" % ep)]; z2 = st[(tg[1], "z_%s" % ep)]
                ok = np.isfinite(z1) & np.isfinite(z2)
                r = np.corrcoef(z1[ok], z2[ok])[0, 1]
                rows.append(dict(source="python", scenario=name, hypothesis=hyp, look_type=lt, look="1-2", time=np.nan,
                                 quantity="z_corr", endpoint=ep, group="", mean=r, sd=np.nan,
                                 se=(1 - r ** 2) / np.sqrt(ok.sum()), n=int(ok.sum())))
        # rejection of the O'Brien-Fleming design at the event looks (OS)
        tg = [t for t in tags if t[0] == "event"]
        z1 = st[(tg[0], "z_os")]; z2 = st[(tg[1], "z_os")]
        r1 = z1 >= b[0]; ra = r1 | (z2 >= b[1])
        for q, x in (("reject_look1", r1), ("reject_any", ra)):
            p = x.mean()
            rows.append(dict(source="python", scenario=name, hypothesis=hyp, look_type="event", look="", time=np.nan,
                             quantity=q, endpoint="os", group="", mean=p, sd=np.nan,
                             se=np.sqrt(p * (1 - p) / x.size), n=x.size))
        for li, bv in enumerate(b, 1):
            rows.append(dict(source="python", scenario=name, hypothesis=hyp, look_type="event", look=li, time=np.nan,
                             quantity="bound", endpoint="os", group="", mean=bv, sd=0.0, se=0.0, n=0))
        rows.append(dict(source="python", scenario=name, hypothesis=hyp, look_type="event", look="", time=np.nan,
                         quantity="not_reached", endpoint="os", group="",
                         mean=float(np.isnan(st[(tg[-1], "time")]).sum()), sd=np.nan, se=np.nan, n=nrep_tot))
        if hyp == "alternative":
            names = ["p_state1", "p_state2", "p_state3", "p_state4", "surv_os", "surv_pfs"]
            for (gname, t), cnt in sorted(lat.items()):
                ng = cnt[-1]
                for j, q in enumerate(names):
                    p = cnt[j] / ng
                    rows.append(dict(source="python", scenario=name, hypothesis=hyp, look_type="latent", look="",
                                     time=t, quantity=q, endpoint="os" if q == "surv_os" else ("pfs" if q == "surv_pfs" else ""),
                                     group=gname, mean=p, sd=np.nan, se=np.sqrt(p * (1 - p) / ng), n=int(ng)))
    return rows


COLS = ["source", "scenario", "hypothesis", "look_type", "look", "time", "target_events", "quantity", "endpoint",
        "group", "mean", "sd", "se", "n"]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenarios", default="V6,V7,V8")
    ap.add_argument("--reps", type=int, default=20000)
    ap.add_argument("--out", default=os.path.join(HERE, "..", "results", "sim_python.csv"))
    a = ap.parse_args()
    allrows = []
    for nm in a.scenarios.split(","):
        t0 = time.time()
        allrows += run(nm, a.reps)
        print("%s done in %.0f s" % (nm, time.time() - t0), flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, lineterminator="\n")
        w.writeheader()
        for r in allrows:
            r = dict(r)
            r.setdefault("target_events", "")
            for k in ("mean", "sd", "se", "time"):
                if isinstance(r[k], float):
                    r[k] = "" if np.isnan(r[k]) else repr(float(r[k]))
            w.writerow(r)
    print("written", a.out, len(allrows), "rows")
