"""Patient-level simulation of the expected study duration, events and enrolled patients when the trial stops at the
bounds (version 0.4.0, Planning/stage4_plan.md, section 2.5), written independently of the R code.

The patients are generated with the functions of simulate_ms.py (model, accrual, log-rank statistic, time of the
d-th death); simulate_ms.py and its output are not changed. Each trial is analysed at its looks with the
O'Brien-Fleming type efficacy bounds (ref_gsd.py) and the fixed futility bound Z = 0 at the interim looks; the
efficacy bounds ignore the futility bound (non-binding), and the trial stops at the first look at which Z crosses
either bound, as in the expected values of the package. Z is the log-rank statistic of OS, positive when the
treatment is better.

Designs:
  E1  scenario V1 (exponential, hazard ratio 0.7, dropout 0.005, 20 patients per month, n = 400), looks at 150, 225
      and 300 OS deaths
  E2  scenario V6 (test scenario C: switching at and after progression, clock reset, n = 500), looks at 150 and 250
      OS deaths
  E3  scenario V1, looks at calendar months 24, 30 and 36; the bounds are computed at the fractions of the expected
      OS deaths at the looks under the alternative (computed here by quadrature), as in the package

Each design is simulated 20000 times under the alternative and under the null hypothesis (the treatment group with
the control model). For each trial: the look and reason of the stop, and the calendar time, the OS deaths and the
enrolled patients at the stop. Only summaries are written, to results/sim_python_expected.csv (columns of
results/sim_python.csv):
  overall: prob_efficacy_early and prob_futility_early (stopping at an interim look), study_duration, events and
    enrolled (means at the stop), and study_duration_indep, events_indep and enrolled_indep: the sums over the looks of
    the simulated probability of stopping at the look times the mean over all trials of the time, the deaths and the
    enrolled patients at the look (the expected values without the correlation between the look time and the
    decision to stop; their se is not computed);
  by look: stop_efficacy and stop_futility (probability of stopping at the look for each reason), look_time,
    events_at_look and enrolled_at_look (means over all trials reaching the look time, whether or not they stopped
    earlier), and the efficacy bound (bound) used;
  not_reached: the number of trials in which an event-driven look is never reached (then the trial stops at the
    last reached look; expected 0).

Usage (numpy and scipy; about 4 minutes):
  python3 simulate_expected.py
  python3 simulate_expected.py --designs E1 --reps 2000 --out /tmp/x.csv
"""
import argparse, csv, os, sys, time
import numpy as np
from scipy import integrate

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ref_gsd  # noqa: E402
import simulate_ms as S  # noqa: E402


def scenario_v1():
    lam = np.log(2) / 12
    return dict(control=S.Group(lam, dropout=0.005), treatment=S.Group(0.7 * lam, dropout=0.005), alloc=1.0, n=400,
                accrual=dict(type="rate", time=[0.0], rate=[20.0]))


def expected_deaths_v1(t):
    """Expected OS deaths of V1 at calendar time t (uniform accrual over 20 months, exponential death and dropout)."""
    lam = np.log(2) / 12
    gam = 0.005
    dur = 20.0

    def p_death(u, h):
        return h / (h + gam) * (1 - np.exp(-(h + gam) * u)) if u > 0 else 0.0

    tot = 0.0
    for h in (lam, 0.7 * lam):
        f = lambda e, h=h: p_death(t - e, h)
        tot += 200 * integrate.quad(f, 0, min(t, dur), epsabs=1e-12, epsrel=1e-12, limit=200)[0] / dur
    return tot


DESIGNS = {
    "E1": lambda: dict(scenario=scenario_v1(), look_type="event", looks=[150, 225, 300]),
    "E2": lambda: dict(scenario=S.SCENARIOS["V6"](), look_type="event", looks=[150, 250]),
    "E3": lambda: dict(scenario=scenario_v1(), look_type="calendar", looks=[24.0, 30.0, 36.0]),
}
SEEDS = {"E1": 11, "E2": 21, "E3": 31}


def design_bounds(d):
    if d["look_type"] == "event":
        frac = np.asarray(d["looks"], float) / d["looks"][-1]
    else:
        ev = np.array([expected_deaths_v1(t) for t in d["looks"]])
        frac = ev / ev[-1]
    _, b = ref_gsd.bounds(frac, alpha=0.025, sf_eff="OF")
    return frac, b


def run(name, reps, chunk_patients=2_000_000, base_seed=20261102):
    d = DESIGNS[name]()
    sc = d["scenario"]
    n = sc["n"]
    p1 = sc["alloc"] / (1 + sc["alloc"])
    brk, mass = S.accrual_setup(sc["accrual"], n)
    frac, b = design_bounds(d)
    kmax = len(d["looks"])
    lower = np.zeros(kmax)
    rows = []
    for hyp in ("alternative", "null"):
        trt_model = sc["treatment"] if hyp == "alternative" else sc["control"]
        per = max(1, chunk_patients // n)
        ss = np.random.SeedSequence([base_seed, SEEDS[name], 0 if hyp == "alternative" else 1])
        nchunk = int(np.ceil(reps / per))
        acc = {}
        for ci, child in enumerate(ss.spawn(nchunk)):
            rng = np.random.default_rng(child)
            nr = min(per, reps - ci * per)
            rep = np.repeat(np.arange(nr), n)
            entry = S.draw_entry(rng, brk, mass, nr * n)
            trt = rng.random(nr * n) < p1
            os_ = np.empty(nr * n)
            drop = np.empty(nr * n)
            for flag, g in ((False, sc["control"]), (True, trt_model)):
                ix = np.where(trt == flag)[0]
                _, os_[ix], _, drop[ix] = S.simulate_group(rng, g, ix.size)
            death_cal = entry + os_
            seen = os_ <= drop
            for k, lv in enumerate(d["looks"]):
                if d["look_type"] == "event":
                    cut = S.kth_time(rep[seen], death_cal[seen], lv, nr)
                else:
                    cut = np.full(nr, lv)
                cut_f = np.where(np.isnan(cut), np.inf, cut)
                f = cut_f[rep] - entry
                enr = f >= 0
                obs_ev = enr & seen & (death_cal <= cut_f[rep])
                ti = np.where(obs_ev, os_, np.minimum(np.minimum(os_, drop), np.maximum(f, 0)))
                z = S.logrank_z(rep[enr], ti[enr], obs_ev[enr], trt[enr], nr)
                for key, v in (("time", cut), ("z", z),
                               ("events", np.bincount(rep, weights=obs_ev, minlength=nr)),
                               ("enrolled", np.bincount(rep, weights=enr, minlength=nr))):
                    acc.setdefault((k, key), []).append(v)
            print("  %s %s chunk %d/%d (%d trials)" % (name, hyp, ci + 1, nchunk, nr), flush=True)
        st = {k: np.concatenate(v) for k, v in acc.items()}
        reached = np.ones(reps, bool)
        for k in range(kmax):
            reached &= np.isfinite(st[(k, "time")])
        not_reached = int((~reached).sum())
        # stopping look and reason (0-based look; reason 1 efficacy, 2 futility, 3 final)
        stop = np.full(reps, -1)
        reason = np.zeros(reps, int)
        for k in range(kmax):
            z = st[(k, "z")]
            act = stop < 0
            eff = act & (z >= b[k])
            fut = act & ~eff & (z < lower[k]) & (k < kmax - 1)
            stop[eff] = k
            reason[eff] = 1
            stop[fut] = k
            reason[fut] = 2
            # final look, or the next look is never reached in the trial
            last = np.ones(reps, bool) if k == kmax - 1 else ~np.isfinite(st[(k + 1, "time")])
            rest = (stop < 0) & last
            stop[rest] = k
            reason[rest] = 3
        pick = lambda key: np.choose(stop, [st[(k, key)] for k in range(kmax)])

        def add(look, quantity, x, se=None, mean=None):
            x = np.asarray(x, float)
            ok = np.isfinite(x)
            m = x[ok].mean() if mean is None else mean
            s = x[ok].std(ddof=1) if ok.sum() > 1 and mean is None else np.nan
            rows.append(dict(source="python", scenario=name, hypothesis=hyp, look_type=d["look_type"], look=look,
                             time=d["looks"][look - 1] if (d["look_type"] == "calendar" and look != "") else np.nan,
                             target_events=d["looks"][look - 1] if (d["look_type"] == "event" and look != "") else "",
                             quantity=quantity, endpoint="os", group="", mean=m, sd=s,
                             se=(s / np.sqrt(ok.sum()) if se is None else se), n=int(ok.sum())))

        interim = stop < kmax - 1
        add("", "prob_efficacy_early", interim & (reason == 1))
        add("", "prob_futility_early", interim & (reason == 2))
        add("", "study_duration", pick("time"))
        add("", "events", pick("events"))
        add("", "enrolled", pick("enrolled"))
        p_stop = np.array([(stop == k).mean() for k in range(kmax)])
        for key, q in (("time", "study_duration_indep"), ("events", "events_indep"), ("enrolled", "enrolled_indep")):
            means = np.array([np.nanmean(st[(k, key)]) for k in range(kmax)])
            add("", q, [np.nan], se=np.nan, mean=float(p_stop @ means))
        for k in range(kmax):
            add(k + 1, "stop_efficacy", (stop == k) & (reason == 1))
            if k < kmax - 1:
                add(k + 1, "stop_futility", (stop == k) & (reason == 2))
            add(k + 1, "look_time", st[(k, "time")])
            add(k + 1, "events_at_look", st[(k, "events")])
            add(k + 1, "enrolled_at_look", st[(k, "enrolled")])
            rows.append(dict(source="python", scenario=name, hypothesis=hyp, look_type=d["look_type"], look=k + 1,
                             time=d["looks"][k] if d["look_type"] == "calendar" else np.nan,
                             target_events=d["looks"][k] if d["look_type"] == "event" else "",
                             quantity="bound", endpoint="os", group="", mean=float(b[k]), sd=0.0, se=0.0, n=0))
        rows.append(dict(source="python", scenario=name, hypothesis=hyp, look_type=d["look_type"], look="",
                         time=np.nan, target_events="", quantity="not_reached", endpoint="os", group="",
                         mean=float(not_reached), sd=np.nan, se=np.nan, n=reps))
    return rows


COLS = ["source", "scenario", "hypothesis", "look_type", "look", "time", "target_events", "quantity", "endpoint",
        "group", "mean", "sd", "se", "n"]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--designs", default="E1,E2,E3")
    ap.add_argument("--reps", type=int, default=20000)
    ap.add_argument("--out", default=os.path.join(HERE, "..", "results", "sim_python_expected.csv"))
    a = ap.parse_args()
    allrows = []
    for nm in a.designs.split(","):
        t0 = time.time()
        allrows += run(nm, a.reps)
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
