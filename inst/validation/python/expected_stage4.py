"""Expected values for the tests of version 0.4.0 (Planning/stage4_plan.md in the project folder): the expected
study duration, events and enrolled patients when the trial stops at the bounds (step 1) and the conditional power
(step 2).

Run from this folder with:  python3 expected_stage4.py [expected | cond_power]
Requires numpy and scipy. Results are written to results/expected_stage4.csv (step 1) and
results/cond_power_stage4.csv (step 2), or only to the file of the argument, and read by tests/testthat
(helper-expected.R, stage4_expected() and stage4_cp()).

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

Conditional power (step 2), with one row per case, quantity and look i (or pair of coordinates i, j):
  cm_*   conditional mean and covariance of a normal vector with three coordinates (.cond_moments() in R)
  B_*    design B at look 1 (z = 0.5) and at look 2 (z = 1) under the design scenario, the null scenario, the
         current trend, the hazard ratio 0.75 and scenario B with the death hazards 0.05 and 0.03 in the treatment
         group (same allocation and accrual, looks at 60, 90 and 120 expected deaths); and at look 2 with z = -1.5,
         below the non-binding futility bound, under the design scenario; B_design has the bounds of the design
  S_*    scenario B, n = 300, OS looks at 40, 80 and 120 events, O'Brien-Fleming type efficacy
         bounds, variance "studentized" (covariance without independent increments, ref_zeta.zeta_cov() with the
         moments of ref_zeta.same_moments() as in expected_design.py), under the design scenario, given the
         statistic of the look or the statistics of the looks up to it; S_design has the analysis times, events,
         means, covariance and efficacy bounds of the design
  M5_*   a fixed mean vector, covariance matrix (not that of a Markov sequence) and bounds of five looks
         (.cond_probs() in R) at look 1 (z = 0.9) and at look 2 (given 1.4, or 0.9 and 1.4)
  G_*    scenario A, n = 400, OS looks at 100, 200 and 300 events, O'Brien-Fleming type efficacy bounds,
         non-binding futility bounds by beta spending (HSD, -2, beta 0.1), mean "schoenfeld", variance "null",
         under the design scenario with the correlations of the information under the design scenario (the
         package, G_*_design) and of the events (gsDesign::gsCP() with theta = -log(0.7) / 2, G_*_events), at
         look 1 (z = 0.5) and at look 2 (z = 1.5); G_design has the analysis times, events, means and bounds
For statistics with the correlations sqrt(s_j / s_k), the conditional moments are written from the independent
increments of W_k = Z_k sqrt(s_k); otherwise from the precision matrix Q of the coordinates involved (covariance
inv(Q_FF), mean mu_F - inv(Q_FF) Q_FG (z - mu_G)), whereas the R code regresses on the given coordinates. The
stopping probabilities of the later looks are ref_gsd.stop_probs() of the conditional moments, with the futility
bounds of the later interim looks (cond_power) and without them (cond_power_nonbinding).
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
import ref_zeta as rz  # noqa: E402
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


# ---------------------------------------------------------------------------
# Step 2: conditional power
# ---------------------------------------------------------------------------
def cond_markov(s, mu, look, z):
    """Conditional moments of the statistics of the looks after 'look' given Z_look = z, for statistics with the
    means mu and the correlations sqrt(s_j / s_k): W_k = Z_k sqrt(s_k) has increments W_k - W_look independent of
    W_look, with means mu_k sqrt(s_k) - mu_look sqrt(s_look) and covariances min(s_k, s_l) - s_look."""
    s = np.asarray(s, float)
    mu = np.asarray(mu, float)
    j = look - 1
    sf = s[look:]
    cm = (z * np.sqrt(s[j]) + mu[look:] * np.sqrt(sf) - mu[j] * np.sqrt(s[j])) / np.sqrt(sf)
    cc = (np.minimum.outer(sf, sf) - s[j]) / np.sqrt(np.outer(sf, sf))
    return cm, cc


def cond_general(mean, sig, given, future, z):
    """Conditional moments of the coordinates 'future' given the coordinates 'given' (0-based) from the precision
    matrix Q of the coordinates given + future: covariance inv(Q_FF), mean mu_F - inv(Q_FF) Q_FG (z - mu_G)."""
    mean = np.asarray(mean, float)
    idx = list(given) + list(future)
    q = np.linalg.inv(np.asarray(sig, float)[np.ix_(idx, idx)])
    g = len(given)
    cc = np.linalg.inv(q[g:, g:])
    cm = mean[list(future)] - cc @ q[g:, :g] @ (np.asarray(z, float) - mean[list(given)])
    return cm, (cc + cc.T) / 2


def is_markov(sig, tol=1e-8):
    sd = np.sqrt(np.diag(sig))
    prec = np.linalg.inv(sig / np.outer(sd, sd))
    k = sig.shape[0]
    off = np.abs(np.subtract.outer(np.arange(k), np.arange(k))) >= 2
    return k <= 2 or np.max(np.abs(prec[off])) <= tol * np.max(np.abs(prec))


def cp_rows(case, look, cm, cc, lo, up, time=None, events=None, extra=None):
    """Rows of a case: the conditional power with and without the futility bounds of the later interim looks and,
    at each later look, the conditional moments and stopping probabilities (and the time and events)."""
    lo = np.asarray(lo, float)
    up = np.asarray(up, float)
    pu, pl = rg.stop_probs(lo, up, cm, cc)
    lo_nb = np.append(np.full(lo.size - 1, -np.inf), lo[-1])
    pu_nb, _ = rg.stop_probs(lo_nb, up, cm, cc)
    rows = [dict(case=case, quantity="cond_power", i="", j="", value=float(np.sum(pu))),
            dict(case=case, quantity="cond_power_nonbinding", i="", j="", value=float(np.sum(pu_nb)))]
    for a in range(len(cm)):
        k = look + 1 + a
        vals = [("cond_mean", cm[a]), ("cond_sd", np.sqrt(cc[a, a])), ("prob_upper", pu[a]), ("prob_lower", pl[a])]
        if time is not None:
            vals += [("time", time[look + a]), ("events", events[look + a])]
        for q, v in vals:
            rows.append(dict(case=case, quantity=q, i=k, j="", value=float(v)))
    for q, v in (extra or {}).items():
        rows.append(dict(case=case, quantity=q, i="", j="", value=float(v)))
    return rows


def cp_moments_cases():
    mean = np.array([0.5, 1.0, 1.5])
    sig = np.array([[1.0, 0.6, 0.4], [0.6, 1.2, 0.7], [0.4, 0.7, 1.5]])
    rows = []
    for case, given, future, z in (("cm_1", [0], [1, 2], [0.2]), ("cm_12", [0, 1], [2], [0.2, 1.4]),
                                   ("cm_13", [0, 2], [1], [0.2, 2.0])):
        cm, cc = cond_general(mean, sig, given, future, z)
        for a in range(len(future)):
            rows.append(dict(case=case, quantity="mean", i=a + 1, j="", value=float(cm[a])))
            for b in range(len(future)):
                rows.append(dict(case=case, quantity="sigma", i=a + 1, j=b + 1, value=float(cc[a, b])))
    return rows


def scenario_b_trt3():
    """Scenario B with the death hazards 0.05 and 0.03 (from months 0 and 6) in the treatment group."""
    b = scenario_b()
    trt = rm.Model(death=[0.05, 0.03], time=[0, 6], dropout=0.012)
    return rl.Scenario(b.control, trt, b.p1 / b.p0, b.accrual)


def cp_design_b():
    _, inner = design_b()
    scn, n = scenario_b(), 300
    ds = inner["design"]
    lo, up = np.asarray(ds["lower"]), np.asarray(ds["upper"])
    times = ds["time"]
    mos = [rl.moments(scn, n, t, "os") for t in times]
    d = np.array([m["events"] for m in mos])
    info = np.array([m["info"] for m in mos])
    mu = np.array([m["z_direct"] for m in mos])
    s0 = null_of(scn)
    l0 = inner["null"]
    v0 = np.array([rl.moments(s0, n, t, "os")["info"] for t in l0["time"]])
    s2 = scenario_b_trt3()
    t2 = [rl.event_time(s2, n, x, "os") for x in (60.0, 90.0, 120.0)]
    m2 = [rl.moments(s2, n, t, "os") for t in t2]
    d2 = [m["events"] for m in m2]
    p01 = scn.p0 * scn.p1
    rows = []
    for k in range(3):
        for q, v in (("z_upper", up[k]), ("z_lower", lo[k])):
            rows.append(dict(case="B_design", quantity=q, i=k + 1, j="", value=float(v)))
    for look, z in ((1, 0.5), (2, 1.0)):
        fl = slice(look, 3)
        h_tr = np.exp(-z / np.sqrt(p01 * d[look - 1]))
        cases = (("design", info, mu, times, d, {}),
                 ("null", v0, np.zeros(3), l0["time"], l0["events"], {}),
                 ("trend", d, z * np.sqrt(d / d[look - 1]), times, d, dict(hr=h_tr)),
                 ("hr", d, -np.log(0.75) * np.sqrt(p01 * d), times, d, dict(hr=0.75)),
                 ("scenario", [m["info"] for m in m2], [m["z_direct"] for m in m2], t2, d2, {}))
        for name, s, m, tm, ev, extra in cases:
            cm, cc = cond_markov(s, m, look, z)
            rows += cp_rows("B_l%d_%s" % (look, name), look, cm, cc, lo[fl], up[fl], tm, ev, extra)
    cm, cc = cond_markov(info, mu, 2, -1.5)
    rows += cp_rows("B_l2low_design", 2, cm, cc, lo[2:], up[2:], times, d)
    return rows


def cp_design_s():
    scn, n = scenario_b(), 300
    targets = [40.0, 80.0, 120.0]
    times = [rl.event_time(scn, n, x, "os") for x in targets]
    mos = [rl.moments(scn, n, t, "os") for t in times]
    d = np.array([m["events"] for m in mos])
    mu = np.array([m["z_direct"] for m in mos])
    mean_, prod = rz.same_moments(scn, n, times, "os")
    core = {("os", t): rz.per_patient_core(scn, n, t, "os") for t in times}
    sig = rz.zeta_cov([("os", t, "os") for t in times], {"os": (list(times), mean_, prod)}, {}, core)
    sig = (sig + sig.T) / 2
    assert not is_markov(sig)
    lo, up = rg.bounds(d / d[-1], corr=corr_of(d))
    rows = []
    for k in range(3):
        for q, v in (("time", times[k]), ("events", d[k]), ("z_mean", mu[k]), ("z_upper", up[k])):
            rows.append(dict(case="S_design", quantity=q, i=k + 1, j="", value=float(v)))
        for m in range(3):
            rows.append(dict(case="S_design", quantity="sigma", i=k + 1, j=m + 1, value=float(sig[k, m])))
    for case, look, z in (("S_l1", 1, [0.8]), ("S_l2_hist", 2, [0.8, 1.2]), ("S_l2_cur", 2, [1.2])):
        given = list(range(look)) if len(z) > 1 else [look - 1]
        cm, cc = cond_general(mu, sig, given, list(range(look, 3)), z)
        rows += cp_rows(case, look, cm, cc, lo[look:], up[look:], times, d)
    return rows


def m5_inputs():
    t = np.arange(1.0, 6.0)
    c = 0.9 * corr_of(t) + 0.1
    sd = np.sqrt([1.0, 1.1, 0.95, 1.05, 1.0])
    sig = c * np.outer(sd, sd)
    mean = np.array([0.8, 1.3, 1.7, 2.0, 2.3])
    lo = np.array([-0.5, 0.2, 0.8, 1.4, 2.0])
    up = np.array([3.5, 2.9, 2.5, 2.2, 2.0])
    return mean, sig, lo, up


def cp_m5():
    mean, sig, lo, up = m5_inputs()
    rows = []
    for case, look, z in (("M5_l1", 1, [0.9]), ("M5_l2_hist", 2, [0.9, 1.4]), ("M5_l2_cur", 2, [1.4])):
        given = list(range(look)) if len(z) > 1 else [look - 1]
        cm, cc = cond_general(mean, sig, given, list(range(look, 5)), z)
        if look == 1:
            assert not is_markov(cc)
        rows += cp_rows(case, look, cm, cc, lo[look:], up[look:])
    return rows


def cp_design_g():
    scn, n = scenario_a(), 400
    targets = [100.0, 200.0, 300.0]
    times = [rl.event_time(scn, n, x, "os") for x in targets]
    mos = [rl.moments(scn, n, t, "os") for t in times]
    d = np.array([m["events"] for m in mos])
    info = np.array([m["info"] for m in mos])
    mu = np.array([m["z_schoenfeld"] for m in mos])
    lo, up = rg.bounds(d / d[-1], futility="beta", sf_fut="HSD", param_fut=-2, beta=0.1, mean=mu,
                       corr=corr_of(d), sigma_alt=corr_of(info))
    rows = []
    for k in range(3):
        for q, v in (("time", times[k]), ("events", d[k]), ("z_mean", mu[k]), ("z_upper", up[k]),
                     ("z_lower", lo[k])):
            rows.append(dict(case="G_design", quantity=q, i=k + 1, j="", value=float(v)))
    for look, z in ((1, 0.5), (2, 1.5)):
        for name, s in (("design", info), ("events", d)):
            cm, cc = cond_markov(s, mu, look, z)
            rows += cp_rows("G_l%d_%s" % (look, name), look, cm, cc, lo[look:], up[look:], times, d)
    return rows


def write_csv(name, fields, rows):
    path = os.path.join(OUT, name)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(dict(r, value=repr(float(r["value"]))))
    print("wrote", path, len(rows), "rows", flush=True)


def main(parts):
    if "expected" in parts:
        rows = (design_b()[0] + design_a_cal()[0] + design_d4()[0] + design_rpact(True)[0] +
                design_rpact(False)[0])
        write_csv("expected_stage4.csv", ["config", "hypothesis", "quantity", "look", "value"], rows)
    if "cond_power" in parts:
        rows = []
        for f in (cp_moments_cases, cp_design_b, cp_design_g, cp_m5, cp_design_s):
            rows += f()
            print(f.__name__, "done", flush=True)
        write_csv("cond_power_stage4.csv", ["case", "quantity", "i", "j", "value"], rows)


if __name__ == "__main__":
    main(sys.argv[1:] or ["expected", "cond_power"])
