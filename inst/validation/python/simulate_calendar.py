"""Patient-level simulation of switching from a calendar time (switch_start of ms_scenario(); stage 3, candidate 3).

Written independently of the R code. Scenarios V13-V15 are scenario C of the tests (switching at and after
progression in the control group, n = 500, accrual of 30% in months 0 to 6 and 70% in months 6 to 18) with switching
possible from calendar month 15 (switch_start = 15):
  V13  clock forward
  V14  clock reset
  V15  clock forward, accelerated failure time effect of the switch (switch_af = 0.6, no death_switch)

A patient enrolled at e can switch from the follow-up time o = max(15 - e, 0) on. The patient is generated as in
simulate_ms.py (entry time, group by simple randomization, exit from the progression-free state by inversion of the
cumulative hazard), with these changes for a patient who progresses at u: the switch at progression (probability
switch_prob) happens only when u >= o; a patient who does not switch at progression has only the death hazard
death_prog until the post-progression clock reaches the value at follow-up time o (o for the forward clock, o - u
for the reset clock), and from then on the competing hazards switch_rate and death_prog; after a switch the death
hazard is death_switch on the post-progression clock (with switch_af = phi, the time from the switch to death is the
time to death without switch divided by phi). Dropout censors OS.

Only summaries are written, to results/sim_python_calendar.csv (columns of sim_python.csv):
  per-patient quantities (look_type "perpatient", calendar months 24 and 36): OS events by group and in total
    ("events"), the log-rank score U = sum over deaths of (treatment indicator - Y1 / Y) ("score") and its variance
    estimate sum of Y0 Y1 / Y^2 ("info"), divided by the number of patients, from 100 trials of 20000 patients
    (mean and standard error over the trials);
  calendar looks (look_type "calendar", n = 500, 20000 trials): log-rank Z of OS at months 24 and 36 (positive when
    the treatment is better) and its correlation between the looks;
  event looks (look_type "event", n = 500, the same trials): calendar time of the 150th and 250th OS death, Z at
    these times, its correlation, and the rejection of the O'Brien-Fleming design (bounds of ref_gsd.py for the
    information fractions 150 / 250 and 1); looks not reached are counted.
Each quantity is simulated under the alternative and under the null hypothesis (the control model in both groups).

Usage (about 3 minutes):
  python3 simulate_calendar.py
  python3 simulate_calendar.py --scenarios V13 --reps 2000 --big-reps 20 --out /tmp/x.csv
"""
import argparse, csv, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ref_gsd  # noqa: E402
import simulate_ms as sm  # noqa: E402

TC = 15.0
SEEDS = {"V13": 131, "V14": 141, "V15": 151}


def scenario(name):
    clock = "reset" if name == "V14" else "forward"
    af = 0.6 if name == "V15" else None
    ctl = sm.Group(0.01, time=[0, 4], prog=[0.10, 0.08], switch_prob=0.3, time_post=[0, 3], death_prog=[0.05, 0.06],
                   switch_rate=0.04, death_switch=None if af is not None else [0.03, 0.035], dropout=0.005,
                   clock=clock, switch_af=af)
    trt = sm.Group(0.008, time=[0, 4], prog=[0.06, 0.05], switch_prob=0.0, time_post=[0, 3], death_prog=[0.045, 0.05],
                   switch_rate=0.0, dropout=0.005, clock=clock)
    return dict(control=ctl, treatment=trt, alloc=1.0, n=500,
                accrual=dict(type="prop", time=[0.0, 6.0], prop=[0.3, 0.7], duration=18.0),
                cal=[24.0, 36.0], events=[150, 250])


def death_after_switch(rng, g, sig):
    """Post-progression clock value at death after a switch at the clock value sig."""
    if g.switch_af is not None:
        t0 = sm.draw_after(rng, sig, g.time_post, g.death_prog)
        return sig + (t0 - sig) / g.switch_af
    if g.switch_clock == "switch":
        return sig + sm.draw_after(rng, np.zeros(sig.size), g.time_post, g.death_switch)
    return sm.draw_after(rng, sig, g.time_post, g.death_switch)


def simulate_group_cal(rng, g, e, tc):
    """OS and dropout times since randomization of patients enrolled at e, switching possible from calendar time tc."""
    n = e.size
    o = np.maximum(tc - e, 0.0)
    tot = g.prog + g.death
    t1 = sm.draw_after(rng, np.zeros(n), g.time, tot)
    with np.errstate(invalid="ignore", divide="ignore"):
        pp = np.where(np.isfinite(t1), sm.step(t1, g.time, g.prog) / sm.step(t1, g.time, tot), 0.0)
    prog = rng.random(n) < pp
    os_ = t1.copy()
    ip = np.where(prog)[0]
    if ip.size:
        u = t1[ip]
        oi = o[ip]
        tau0 = np.zeros(ip.size) if g.clock == "reset" else u.copy()
        # post-progression clock value from which switching is possible
        avail = np.maximum(oi - u, 0.0) if g.clock == "reset" else np.maximum(oi, u)
        now = (u >= oi) & (rng.random(ip.size) < sm.step(u, g.time, g.switch_prob))
        # before avail only death (death_prog); from avail on switch (switch_rate) or death (death_prog)
        td_a = sm.draw_after(rng, tau0, g.time_post, g.death_prog)
        dead_a = ~now & (td_a <= avail)
        y = g.switch_rate + g.death_prog
        tau2 = sm.draw_after(rng, avail, g.time_post, y)
        with np.errstate(invalid="ignore", divide="ignore"):
            ps = np.where(np.isfinite(tau2), sm.step(tau2, g.time_post, g.switch_rate) / sm.step(tau2, g.time_post, y),
                          0.0)
        later = ~now & ~dead_a & (rng.random(ip.size) < ps)
        switched = now | later
        sig = np.where(now, tau0, np.where(later, tau2, 0.0))
        taud_s = death_after_switch(rng, g, sig)
        taud = np.where(switched, taud_s, np.where(dead_a, td_a, tau2))
        os_[ip] = u + (taud - tau0)
    drop = sm.draw_after(rng, np.zeros(n), g.time, g.dropout)
    return os_, drop


def score_info(rep, t, ev, trt, nrep):
    """Log-rank score U (treatment group) and its variance estimate for each trial."""
    o = np.lexsort((t, rep))
    r = rep[o]; e = ev[o]; g = trt[o].astype(float)
    nt = np.bincount(r, minlength=nrep)
    start = np.concatenate([[0], np.cumsum(nt)[:-1]])
    pos = np.arange(r.size) - start[r]
    y = (nt[r] - pos).astype(float)
    cg = np.concatenate([np.cumsum(g[::-1])[::-1], [0.0]])
    yt = cg[:-1] - cg[(start + nt)[r]]
    u = np.bincount(r[e], weights=g[e] - yt[e] / y[e], minlength=nrep)
    v = np.bincount(r[e], weights=yt[e] * (y[e] - yt[e]) / y[e] ** 2, minlength=nrep)
    return u, v


def draw_trials(rng, sc, trt_model, n, nr, brk, mass, p1):
    rep = np.repeat(np.arange(nr), n)
    entry = sm.draw_entry(rng, brk, mass, nr * n)
    trt = rng.random(nr * n) < p1
    os_ = np.empty(nr * n); drop = np.empty(nr * n)
    for flag, g in ((False, sc["control"]), (True, trt_model)):
        ix = np.where(trt == flag)[0]
        os_[ix], drop[ix] = simulate_group_cal(rng, g, entry[ix], TC)
    return rep, entry, trt, os_, drop


def observed(cut, rep, entry, os_, drop):
    f = cut[rep] - entry
    enr = f >= 0
    # compared on the calendar scale, as in kth_time()
    ev = enr & (os_ <= drop) & (entry + os_ <= cut[rep])
    ti = np.where(ev, os_, np.minimum(np.minimum(os_, drop), np.maximum(f, 0)))
    return enr, ev, ti


def run(name, reps, big_reps, big_n=20000, chunk_patients=2_000_000, base_seed=20261026):
    sc = scenario(name)
    n = sc["n"]; p1 = sc["alloc"] / (1 + sc["alloc"])
    frac = np.asarray(sc["events"], float) / sc["events"][-1]
    _, b = ref_gsd.bounds(frac, alpha=0.025, sf_eff="OF")
    rows = []
    for hyp in ("alternative", "null"):
        trt_model = sc["treatment"] if hyp == "alternative" else sc["control"]
        ss = np.random.SeedSequence([base_seed, SEEDS[name], 0 if hyp == "alternative" else 1])
        s_big, s_trials = ss.spawn(2)

        def add(look_type, look, tval, quantity, group, mean, sd, se, nn, target=""):
            rows.append(dict(source="python", scenario=name, hypothesis=hyp, look_type=look_type, look=look,
                             time=tval, target_events=target, quantity=quantity, endpoint="os", group=group,
                             mean=mean, sd=sd, se=se, n=nn))

        # per-patient quantities from large trials
        brk, mass = sm.accrual_setup(sc["accrual"], big_n)
        big = {}
        for bi, child in enumerate(s_big.spawn(big_reps)):
            rng = np.random.default_rng(child)
            rep, entry, trt, os_, drop = draw_trials(rng, sc, trt_model, big_n, 1, brk, mass, p1)
            for c in sc["cal"]:
                enr, ev, ti = observed(np.full(1, c), rep, entry, os_, drop)
                u, v = score_info(rep[enr], ti[enr], ev[enr], trt[enr], 1)
                for key, val in (("events_control", np.sum(ev & ~trt)), ("events_treatment", np.sum(ev & trt)),
                                 ("events", np.sum(ev)), ("score", u[0]), ("info", v[0])):
                    big.setdefault((c, key), []).append(val / big_n)
        for c in sc["cal"]:
            for key, q, g in (("events_control", "events", "control"), ("events_treatment", "events", "treatment"),
                              ("events", "events", ""), ("score", "score", ""), ("info", "info", "")):
                x = np.asarray(big[(c, key)])
                add("perpatient", "", c, q, g, x.mean(), x.std(ddof=1), x.std(ddof=1) / np.sqrt(x.size), x.size)
        print("  %s %s per-patient done" % (name, hyp), flush=True)

        # trials of n patients
        brk, mass = sm.accrual_setup(sc["accrual"], n)
        per = max(1, chunk_patients // n)
        nchunk = int(np.ceil(reps / per))
        acc = {}
        for ci, child in enumerate(s_trials.spawn(nchunk)):
            rng = np.random.default_rng(child)
            nr = min(per, reps - ci * per)
            rep, entry, trt, os_, drop = draw_trials(rng, sc, trt_model, n, nr, brk, mass, p1)
            for li, c in enumerate(sc["cal"], 1):
                enr, ev, ti = observed(np.full(nr, c), rep, entry, os_, drop)
                acc.setdefault(("calendar", li), []).append(sm.logrank_z(rep[enr], ti[enr], ev[enr], trt[enr], nr))
            death_cal = entry + os_
            seen = os_ <= drop
            for li, d in enumerate(sc["events"], 1):
                cut = sm.kth_time(rep[seen], death_cal[seen], d, nr)
                acc.setdefault(("time", li), []).append(cut)
                cut_f = np.where(np.isnan(cut), np.inf, cut)
                enr, ev, ti = observed(cut_f, rep, entry, os_, drop)
                acc.setdefault(("event", li), []).append(sm.logrank_z(rep[enr], ti[enr], ev[enr], trt[enr], nr))
            print("  %s %s chunk %d/%d (%d trials)" % (name, hyp, ci + 1, nchunk, nr), flush=True)
        st = {k: np.concatenate(v) for k, v in acc.items()}

        def add_x(look_type, look, tval, quantity, x, target=""):
            x = np.asarray(x, float); ok = np.isfinite(x)
            s = x[ok].std(ddof=1)
            add(look_type, look, tval, quantity, "", x[ok].mean(), s, s / np.sqrt(ok.sum()), int(ok.sum()), target)

        def add_corr(look_type, x, y):
            ok = np.isfinite(x) & np.isfinite(y)
            r = np.corrcoef(x[ok], y[ok])[0, 1]
            add(look_type, "1-2", np.nan, "z_corr", "", r, np.nan, (1 - r ** 2) / np.sqrt(ok.sum()), int(ok.sum()))

        for li, c in enumerate(sc["cal"], 1):
            add_x("calendar", li, c, "z", st[("calendar", li)])
        add_corr("calendar", st[("calendar", 1)], st[("calendar", 2)])
        for li, d in enumerate(sc["events"], 1):
            add_x("event", li, np.nan, "time", st[("time", li)], d)
            add_x("event", li, np.nan, "z", st[("event", li)], d)
        z1, z2 = st[("event", 1)], st[("event", 2)]
        add_corr("event", z1, z2)
        r1 = z1 >= b[0]
        ra = r1 | (z2 >= b[1])
        for q, x in (("reject_look1", r1), ("reject_any", ra)):
            p = x.mean()
            add("event", "", np.nan, q, "", p, np.nan, np.sqrt(p * (1 - p) / x.size), x.size)
        for li, bv in enumerate(b, 1):
            add("event", li, np.nan, "bound", "", bv, 0.0, 0.0, 0)
        add("event", "", np.nan, "not_reached", "", float(np.isnan(st[("time", len(sc["events"]))]).sum()), np.nan,
            np.nan, reps)
    return rows


COLS = ["source", "scenario", "hypothesis", "look_type", "look", "time", "target_events", "quantity", "endpoint",
        "group", "mean", "sd", "se", "n"]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenarios", default="V13,V14,V15")
    ap.add_argument("--reps", type=int, default=20000)
    ap.add_argument("--big-reps", type=int, default=100)
    ap.add_argument("--out", default=os.path.join(HERE, "..", "results", "sim_python_calendar.csv"))
    a = ap.parse_args()
    allrows = []
    for nm in a.scenarios.split(","):
        t0 = time.time()
        allrows += run(nm, a.reps, a.big_reps)
        print("%s done in %.0f s" % (nm, time.time() - t0), flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, lineterminator="\n")
        w.writeheader()
        for r in allrows:
            r = dict(r)
            for k in ("mean", "sd", "se", "time"):
                if isinstance(r[k], float):
                    r[k] = "" if np.isnan(r[k]) else repr(float(r[k]))
            w.writerow(r)
    print("written", a.out, len(allrows), "rows")
