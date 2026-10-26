"""Expected values for the tests of switching from a calendar time (switch_start, stage 3).

Uses ref_calendar.py (adaptive quadrature and matrix exponentials, written independently of the R code) and
writes
  results/calendar_os.csv          OS given the availability o for eight post-progression models (os_quad)
  results/calendar_group.csv       group functions Y and D at t = 36 (forward clock: os_markov over the whole
                                   entry range; reset clock: os_quad)
  results/calendar_looks.csv       per-patient events, dropouts, score mean and information at months 24 and 36
  results/calendar_zeta.csv        covariance matrices of the standardized statistics (studentized, alternative,
                                   null), calendar-driven looks and looks driven by OS events
  results/calendar_delta.csv       null hypothesis: correlation between the looks from the closed expression of
                                   delta (calendar note, equation (3))
The designs are scenario C of tests/testthat/helper-scenarios.R (scn_cal(): forward clock, switching of the
control group possible from month 15; the null design has the control model in both groups).

Usage (from this folder; about 15 minutes): python3 expected_calendar.py
"""
import os
import time

import numpy as np

import ref_calendar as RC
import ref_model5 as R5

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
TC = 15.0
LOOKS = [24.0, 36.0]


def models_c(clock):
    """Scenario C of the tests (models_c() of helper-scenarios.R)."""
    ctl = R5.Model5(0.01, time=[0, 4], prog=[0.10, 0.08], switch_prob=0.3, time_post=[0, 3],
                    death_prog=[0.05, 0.06], switch_rate=0.04, death_switch=[0.03, 0.035], dropout=0.005, clock=clock)
    trt = R5.Model5(0.008, time=[0, 4], prog=[0.06, 0.05], switch_prob=0.0, time_post=[0, 3],
                    death_prog=[0.045, 0.05], switch_rate=0.0, dropout=0.005, clock=clock)
    return ctl, trt


def models_post():
    """Eight post-progression models (two clocks; switch clock "same" or "switch", AFT effect 0.6 or 1.5), as
    model_post_cal() of helper-scenarios.R."""
    base = dict(time=[0, 4], prog=[0.10, 0.08], switch_prob=0.3, time_post=[0, 3, 8],
                death_prog=[0.05, 0.08, 0.04], switch_rate=[0.04, 0.06, 0.02], dropout=0.005)
    out = []
    for clock in ("reset", "forward"):
        out.append(("%s_same" % clock, R5.Model5(0.01, death_switch=[0.03, 0.05, 0.02], clock=clock, **base)))
        out.append(("%s_switch" % clock, R5.Model5(0.01, death_switch=[0.03, 0.05, 0.02], clock=clock,
                                                    switch_clock="switch", **base)))
        out.append(("%s_aft0.6" % clock, R5.Model5(0.01, switch_af=0.6, clock=clock, **base)))
        out.append(("%s_aft1.5" % clock, R5.Model5(0.01, switch_af=1.5, clock=clock, **base)))
    return out


def design(null=False, clock="forward", h=0.5):
    ctl, trt = models_c(clock)
    return RC.Design(ctl, ctl if null else trt, 0.5, np.array([0.0, 6.0]), np.array([0.3, 0.7]), 18.0, TC, h=h)


def write(name, header, rows):
    with open(os.path.join(OUT, name), "w", newline="\n") as f:
        f.write(",".join(header) + "\n")
        for r in rows:
            f.write(",".join(x if isinstance(x, str) else repr(float(x)) for x in r) + "\n")


def main():
    t_all = time.time()
    pairs = [(20.0, 7.5), (30.0, 2.0), (12.0, 11.0), (9.0, 15.0), (36.0, 0.5), (25.0, 0.0)]
    rows = []
    for name, m in models_post():
        for s, o in pairs:
            sv, dv = RC.os_quad(m, s, o)
            rows.append((name, s, o, sv, dv))
        if name == "forward_same":
            mk = RC.os_markov(m, np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs]))
            dif = max(np.max(np.abs(mk[0] - np.array([r[3] for r in rows[-len(pairs):]]))),
                      np.max(np.abs(mk[1] - np.array([r[4] for r in rows[-len(pairs):]]))))
            print("forward_same: matrix exponentials against adaptive quadrature, max abs difference %.1e" % dif)
    write("calendar_os.csv", ["model", "s", "o", "surv", "dens"], rows)
    print("calendar_os.csv: %d rows (%.0f s)" % (len(rows), time.time() - t_all))

    rows = []
    des = design()
    s = np.array([3.0, 12.0, 20.0, 30.0])
    for j in range(2):
        y, d = des.group(j, 36.0, s)
        for i in range(s.size):
            rows.append(("forward", j, 36.0, s[i], y[i], d[i]))
    ctl_r, _ = models_c("reset")
    a = np.concatenate([des.alpha])
    for si in (12.0, 30.0):
        hi = min(36.0 - si, des.dur)
        br = list(a) + [TC, TC - si] + list(TC - ctl_r.time) + list(TC - ctl_r.time_post) + list(TC - si + ctl_r.time_post)
        g = RC.rule(0, hi, br, 1.0, 10)
        vals = np.array([RC.os_quad(ctl_r, si, max(TC - e, 0)) for e in g["x"]])
        w = g["w"] * des.dens_e(g["x"])
        rows.append(("reset", 0, 36.0, si, np.sum(w * vals[:, 0]), np.sum(w * vals[:, 1])))
    write("calendar_group.csv", ["clock", "group", "t", "s", "Y", "D"], rows)
    print("calendar_group.csv: %d rows (%.0f s)" % (len(rows), time.time() - t_all))

    rows_l, rows_z = [], []
    for dname, dd in (("alt", design()), ("null", design(null=True))):
        for trig in ("none", "os"):
            z = dd.zeta(LOOKS, trigger=None if trig == "none" else ["os"] * len(LOOKS))
            if trig == "none":
                for i, L in enumerate(z["looks"]):
                    rows_l.append((dname, L["t"], L["ev"][0], L["ev"][1], L["do"][0], L["do"][1], L["mu"], L["v"]))
            for q in ("studentized", "alternative", "null"):
                for i in range(len(LOOKS)):
                    for k in range(len(LOOKS)):
                        rows_z.append((dname, trig, q, i + 1, k + 1, z[q][i, k]))
            print("%s %s done (%.0f s)" % (dname, trig, time.time() - t_all))
    write("calendar_looks.csv", ["design", "t", "ev0", "ev1", "do0", "do1", "mu", "v"], rows_l)
    write("calendar_zeta.csv", ["design", "trigger", "variance", "i", "j", "value"], rows_z)

    dl, e1, e2 = RC.delta_null(design(null=True), LOOKS[0], LOOKS[1])
    corr = (e1 + dl) / np.sqrt(e1 * e2)
    write("calendar_delta.csv", ["t1", "t2", "delta", "events1", "events2", "corr"],
          [(LOOKS[0], LOOKS[1], dl, e1, e2, corr)])
    zn = [r for r in rows_z if r[0] == "null" and r[1] == "none" and r[2] == "null" and r[3] == 1 and r[4] == 2][0][5]
    print("null correlation: delta formula %.10f, moments %.10f, difference %.1e" % (corr, zn, corr - zn))
    print("Total %.0f s" % (time.time() - t_all))


if __name__ == "__main__":
    main()
