"""Expected values for the tests of version 0.4.0 (Planning/stage4_plan.md in the project folder): the expected
study duration, events and enrolled patients when the trial stops at the bounds (step 1).

Run from this folder with:  python3 expected_stage4.py
Requires numpy and scipy. Results are written to results/expected_stage4.csv and read by tests/testthat
(helper-expected.R, stage4_expected()).

Designs (defined identically in tests/testthat/test-ms_power.R):
  B      scenario B of expected_values.py, n = 300, OS looks at 50, 75, 100 per cent of 120 events,
         O'Brien-Fleming type efficacy bounds, non-binding futility bounds by beta spending (HSD, -2, beta 0.1),
         variance "null" (the design of results/power_B.csv)
  A_cal  scenario A, n = 400, calendar looks at months 24 and 36, binding futility bound hazard ratio 1
         (Z = 0) at look 1, variance "null"
  D4     design D4 of expected_design.py (scenario F, n = 400, OS tested at a look driven by 250 PFS events and at
         220 OS events, variance "studentized"); the analysis times, events, bounds and stopping probabilities
         under the design scenario are those of results/design2.csv
  R_binding, R_nonbinding
         the design of the rpact vignette "Designing Group Sequential Trials with Two Groups and a Survival
         Endpoint with rpact" (piecewise exponential death hazards 0.025, 0.04, 0.015, 0.01, 0.007 from months
         0, 6, 9, 15, 21 in the control group, hazard ratio 0.7, exponential dropout of 5 per cent in 12 months,
         accrual of 42 per month, n = 1000, OS looks at 50, 75, 100 per cent of 387 events, O'Brien-Fleming type
         efficacy bounds) with the fixed futility bounds Z = 0 and 0.5, binding and non-binding, mean
         "schoenfeld" and variance "null" (correlations sqrt(V_j / V_k) of the information under the design
         scenario)

Each design is evaluated under the design scenario and under the null scenario (the control model in both groups),
where event-driven looks take place when the expected events of the trigger endpoint reach the same numbers and
calendar-driven looks keep their times:
  analysis times      root finding on ref_logrank.events() (adaptive quadrature)
  enrolled patients   ref_logrank.events()["enrolled"]
  bounds              ref_gsd.bounds() (recursive conditioning with adaptive quadrature)
  stopping probs      ref_gsd.stop_probs() with the efficacy and futility bounds; under the null scenario with the
                      correlations sqrt(V_j / V_k) of the information under the null scenario
  expected values     sum over the looks of the stopping probability times the value at the look, with the
                      probability of stopping at the last look equal to one minus the probability of stopping
                      earlier
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ref_gsd as rg  # noqa: E402
import ref_logrank as rl  # noqa: E402
import ref_model as rm  # noqa: E402
from expected_values import scenario_a, scenario_b  # noqa: E402
from expected_zeta import scenario_f  # noqa: E402

OUT = os.path.join(HERE, "results")


def corr_of(v):
    v = np.asarray(v, float)
    return np.sqrt(np.minimum.outer(v, v) / np.maximum.outer(v, v))


def null_of(scn):
    return rl.Scenario(scn.control, scn.control, scn.p1 / scn.p0, scn.accrual)


def expected(pu, pl, time, events, enrolled):
    """Expected values written out look by look."""
    k = len(time)
    stop = [pu[j] + pl[j] for j in range(k - 1)]
    stop.append(1.0 - sum(stop))
    out = dict(prob_efficacy_early=float(sum(pu[:k - 1])), prob_futility_early=float(sum(pl[:k - 1])))
    for name, x in (("study_duration", time), ("events", events), ("enrolled", enrolled)):
        out[name] = float(sum(stop[j] * x[j] for j in range(k)))
    return out


def rows_of(config, ex_design, ex_null, looks_null):
    rows = []
    for hyp, ex in (("design", ex_design), ("null", ex_null)):
        for q, v in ex.items():
            rows.append(dict(config=config, hypothesis=hyp, quantity=q, look="", value=v))
    for q in ("time", "enrolled", "events", "prob_upper", "prob_lower"):
        for k, v in enumerate(looks_null[q]):
            rows.append(dict(config=config, hypothesis="null", quantity=q, look=k + 1, value=float(v)))
    return rows


def design_b():
    scn, n = scenario_b(), 300
    targets = [60.0, 90.0, 120.0]
    times = [rl.event_time(scn, n, d, "os") for d in targets]
    mos = [rl.moments(scn, n, t, "os") for t in times]
    d = np.array([m["events"] for m in mos])
    info = np.array([m["info"] for m in mos])
    mu = np.array([m["z_direct"] for m in mos])
    lo, up = rg.bounds(d / d[-1], futility="beta", sf_fut="HSD", param_fut=-2, beta=0.1, mean=mu,
                       corr=corr_of(d), sigma_alt=corr_of(info))
    pu, pl = rg.stop_probs(lo, up, mu, corr_of(info))
    enr = [rl.events(scn, n, t, "os")["enrolled"] for t in times]
    s0 = null_of(scn)
    t0 = [rl.event_time(s0, n, x, "os") for x in targets]
    e0 = [rl.events(s0, n, t, "os") for t in t0]
    v0 = [rl.moments(s0, n, t, "os")["info"] for t in t0]
    pu0, pl0 = rg.stop_probs(lo, up, np.zeros(3), corr_of(v0))
    l0 = dict(time=t0, enrolled=[e["enrolled"] for e in e0], events=[e["events"] for e in e0],
              prob_upper=pu0, prob_lower=pl0)
    inner = dict(design=dict(pu=pu, pl=pl, time=times, events=d, enrolled=enr, upper=up, lower=lo,
                             alpha_spent=rg.spend(d / d[-1], 0.025, "OF")), null=l0)
    return rows_of("B", expected(pu, pl, times, d, enr), expected(pu0, pl0, t0, l0["events"], l0["enrolled"]),
                   l0), inner


def design_a_cal():
    scn, n = scenario_a(), 400
    times = [24.0, 36.0]
    mos = [rl.moments(scn, n, t, "os") for t in times]
    d = np.array([m["events"] for m in mos])
    info = np.array([m["info"] for m in mos])
    mu = np.array([m["z_direct"] for m in mos])
    lo, up = rg.bounds(d / d[-1], futility="fixed", fut_bound=[0.0], binding=True, corr=corr_of(d))
    pu, pl = rg.stop_probs(lo, up, mu, corr_of(info))
    enr = [rl.events(scn, n, t, "os")["enrolled"] for t in times]
    s0 = null_of(scn)
    e0 = [rl.events(s0, n, t, "os") for t in times]
    v0 = [rl.moments(s0, n, t, "os")["info"] for t in times]
    pu0, pl0 = rg.stop_probs(lo, up, np.zeros(2), corr_of(v0))
    l0 = dict(time=times, enrolled=[e["enrolled"] for e in e0], events=[e["events"] for e in e0],
              prob_upper=pu0, prob_lower=pl0)
    inner = dict(design=dict(pu=pu, pl=pl, time=times, events=d, enrolled=enr, upper=up, lower=lo), null=l0)
    return rows_of("A_cal", expected(pu, pl, times, d, enr),
                   expected(pu0, pl0, times, l0["events"], l0["enrolled"]), l0), inner


def design_d4():
    with open(os.path.join(OUT, "design2.csv"), newline="") as fh:
        ex = [r for r in csv.DictReader(fh) if r["config"] == "D4"]

    def val(q, look=None):
        return [float(r["value"]) for r in ex if r["quantity"] == q and (look is None or r["look"] == str(look))]
    scn, n = scenario_f(), 400
    times = val("time")
    d = [val("events", k)[0] for k in (1, 2)]
    up = np.array([val("z_upper", k)[0] for k in (1, 2)])
    lo = np.array([val("z_lower", k)[0] for k in (1, 2)])
    pu = [val("prob_upper", k)[0] for k in (1, 2)]
    pl = [val("prob_lower", k)[0] for k in (1, 2)]
    enr = [rl.events(scn, n, t, "os")["enrolled"] for t in times]
    s0 = null_of(scn)
    t0 = [rl.event_time(s0, n, 250.0, "pfs"), rl.event_time(s0, n, 220.0, "os")]
    for a, b in zip(t0, val("time0")):
        assert abs(a - b) < 1e-8, (a, b)
    e0 = [rl.events(s0, n, t, "os") for t in t0]
    v0 = [rl.moments(s0, n, t, "os")["info"] for t in t0]
    pu0, pl0 = rg.stop_probs(lo, up, np.zeros(2), corr_of(v0))
    l0 = dict(time=t0, enrolled=[e["enrolled"] for e in e0], events=[e["events"] for e in e0],
              prob_upper=pu0, prob_lower=pl0)
    inner = dict(design=dict(pu=pu, pl=pl, time=times, events=d, enrolled=enr, upper=up, lower=lo), null=l0)
    return rows_of("D4", expected(pu, pl, times, d, enr), expected(pu0, pl0, t0, l0["events"], l0["enrolled"]),
                   l0), inner


def scenario_rpact():
    haz = [0.025, 0.04, 0.015, 0.01, 0.007]
    tm = [0, 6, 9, 15, 21]
    drop = -np.log(0.95) / 12
    return rl.Scenario(rm.Model(death=haz, time=tm, dropout=drop),
                       rm.Model(death=[0.7 * h for h in haz], time=tm, dropout=drop), 1.0,
                       rl.Accrual(time=[0], rate=[42]))


def design_rpact(binding):
    scn, n = scenario_rpact(), 1000
    targets = [193.5, 290.25, 387.0]
    times = [rl.event_time(scn, n, x, "os") for x in targets]
    mos = [rl.moments(scn, n, t, "os") for t in times]
    d = np.array([m["events"] for m in mos])
    info = np.array([m["info"] for m in mos])
    mu = np.array([m["z_schoenfeld"] for m in mos])
    lo, up = rg.bounds(d / d[-1], futility="fixed", fut_bound=[0.0, 0.5], binding=binding, corr=corr_of(d))
    pu, pl = rg.stop_probs(lo, up, mu, corr_of(info))
    enr = [rl.events(scn, n, t, "os")["enrolled"] for t in times]
    s0 = null_of(scn)
    t0 = [rl.event_time(s0, n, x, "os") for x in targets]
    e0 = [rl.events(s0, n, t, "os") for t in t0]
    v0 = [rl.moments(s0, n, t, "os")["info"] for t in t0]
    pu0, pl0 = rg.stop_probs(lo, up, np.zeros(3), corr_of(v0))
    l0 = dict(time=t0, enrolled=[e["enrolled"] for e in e0], events=[e["events"] for e in e0],
              prob_upper=pu0, prob_lower=pl0)
    name = "R_binding" if binding else "R_nonbinding"
    inner = dict(design=dict(pu=pu, pl=pl, time=times, events=d, enrolled=enr, upper=up, lower=lo, mean=mu),
                 null=l0)
    return rows_of(name, expected(pu, pl, times, d, enr), expected(pu0, pl0, t0, l0["events"], l0["enrolled"]),
                   l0), inner


def main():
    rows = (design_b()[0] + design_a_cal()[0] + design_d4()[0] + design_rpact(True)[0] +
            design_rpact(False)[0])
    path = os.path.join(OUT, "expected_stage4.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["config", "hypothesis", "quantity", "look", "value"], lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(dict(r, value=repr(float(r["value"]))))
    print("wrote", path, len(rows), "rows")


if __name__ == "__main__":
    main()
