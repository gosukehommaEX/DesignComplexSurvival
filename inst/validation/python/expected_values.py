"""Expected values for the tests of DesignComplexSurvival.

Run from this folder with:  python3 expected_values.py
Requires numpy and scipy. Results are written to results/*.csv and are
copied into tests/testthat with the source file and row given in comments.

The scenarios are defined identically in the R tests:

A  exponential, HR 0.7, uniform accrual 20/month, n = 400, 1:1
B  piecewise exponential with a delayed effect, piecewise accrual,
   n = 300, 2:1
C  progression with switching (clock "reset"), accrual by proportions,
   n = 500, 1:1;  Cf is C with clock "forward"
"""

import csv
import os
import sys
import time

import numpy as np
from scipy import optimize, stats

import ref_gsd as rg
import ref_logrank as rl
import ref_model as rm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
LN2 = np.log(2.0)


def write(name, rows):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    keys = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (repr(float(v)) if isinstance(v, (float, np.floating)) else v)
                        for k, v in r.items()})
    print("wrote", path, flush=True)


def scenario_a():
    ctl = rm.Model(death=LN2 / 12, dropout=0.005)
    trt = rm.Model(death=0.7 * LN2 / 12, dropout=0.005)
    return rl.Scenario(ctl, trt, 1.0, rl.Accrual(time=[0], rate=[20]))


def scenario_b():
    ctl = rm.Model(death=[0.05, 0.04], time=[0, 6], dropout=0.01)
    trt = rm.Model(death=[0.05, 0.025], time=[0, 6], dropout=0.012)
    return rl.Scenario(ctl, trt, 2.0, rl.Accrual(time=[0, 3], rate=[10, 30]))


def models_c(clock):
    ctl = rm.Model(time=[0, 4], prog=[0.10, 0.08], death=0.01, switch_prob=0.3,
                   time_post=[0, 3], death_prog=[0.05, 0.06], switch_rate=0.04,
                   death_switch=[0.03, 0.035], dropout=0.005, clock=clock)
    trt = rm.Model(time=[0, 4], prog=[0.06, 0.05], death=0.008, switch_prob=0.0,
                   time_post=[0, 3], death_prog=[0.045, 0.05], switch_rate=0.0,
                   death_switch=None, dropout=0.005, clock=clock)
    return ctl, trt


def scenario_c(clock="reset"):
    ctl, trt = models_c(clock)
    return rl.Scenario(ctl, trt, 1.0,
                       rl.Accrual(time=[0, 6], prop=[0.3, 0.7], duration=18.0))


def main(which):
    t_start = time.time()

    if "model" in which:
        rows = []
        ctl_r, trt_r = models_c("reset")
        ctl_f, trt_f = models_c("forward")
        for label, m in (("C_ctl_reset", ctl_r), ("C_trt_reset", trt_r),
                         ("C_ctl_forward", ctl_f)):
            for t in (2.0, 5.0, 10.0, 20.0):
                st = rm.states(m, t)
                if m.clock == "forward":
                    st2 = rm.conv_states(m, t)
                    assert np.max(np.abs(st - st2)) < 1e-10, (st, st2)
                rows.append(dict(model=label, time=t, p1=st[0], p2=st[1], p3=st[2],
                                 p4=st[3], os_surv=rm.os_surv(m, t), os_dens=rm.os_dens(m, t),
                                 pfs_surv=rm.pfs_surv(m, t), pfs_dens=rm.pfs_dens(m, t)))
        write("model_states.csv", rows)

    if "events" in which:
        rows = []
        scn = scenario_a()
        for t in (12.0, 24.0, 36.0):
            e = rl.events(scn, 400, t, "os")
            # closed form for exponential event and dropout with uniform accrual
            A = 20.0
            for g, lam in (("control", LN2 / 12), ("treatment", 0.7 * LN2 / 12)):
                h = lam + 0.005
                u = min(t, A)
                cf = 200 * (lam / h) * (u - (np.exp(-h * (t - u)) - np.exp(-h * t)) / h) / A
                assert abs(cf - e["events_" + g]) < 1e-8, (cf, e["events_" + g])
            rows.append(dict(scenario="A", endpoint="os", n=400, time=t, **e))
        scn = scenario_b()
        for t in (15.0, 25.0):
            rows.append(dict(scenario="B", endpoint="os", n=300, time=t,
                             **rl.events(scn, 300, t, "os")))
        scn = scenario_c("reset")
        for ep in ("pfs", "os"):
            rows.append(dict(scenario="C", endpoint=ep, n=500, time=20.0,
                             **rl.events(scn, 500, 20.0, ep)))
        scn = scenario_c("forward")
        rows.append(dict(scenario="Cf", endpoint="os", n=500, time=20.0,
                         **rl.events(scn, 500, 20.0, "os")))
        write("events.csv", rows)

    if "event_time" in which:
        rows = []
        for label, scn, n, target in (("A", scenario_a(), 400, 150.0),
                                      ("B", scenario_b(), 300, 60.0),
                                      ("B", scenario_b(), 300, 120.0)):
            rows.append(dict(scenario=label, endpoint="os", n=n, events=target,
                             time=rl.event_time(scn, n, target, "os")))
        write("event_time.csv", rows)

    if "moments" in which:
        rows = []
        for label, scn, n, ep, times in (("A", scenario_a(), 400, "os", (24.0, 36.0)),
                                         ("B", scenario_b(), 300, "os", (15.0, 25.0)),
                                         ("C", scenario_c("reset"), 500, "pfs", (10.0, 20.0)),
                                         ("C", scenario_c("reset"), 500, "os", (20.0,)),
                                         ("Cf", scenario_c("forward"), 500, "os", (20.0,))):
            for t in times:
                mo = rl.moments(scn, n, t, ep)
                rows.append(dict(scenario=label, endpoint=ep, n=n, time=t, **mo))
                print(label, ep, t, mo, flush=True)
        write("moments.csv", rows)

    if "alt" in which:
        rows = []
        scn = scenario_b()
        v15 = rl.luo_var(scn, 300, 15.0, "os")
        v25 = rl.luo_var(scn, 300, 25.0, "os")
        cv = rl.luo_cov(scn, 300, 25.0, 15.0, "os")
        rows.append(dict(scenario="B", endpoint="os", n=300, t1=15.0, t2=15.0, U1=v15["U"],
                         U2=v15["U"], I1=v15["I"], I2=v15["I"], cov=v15["V"]))
        rows.append(dict(scenario="B", endpoint="os", n=300, t1=25.0, t2=25.0, U1=v25["U"],
                         U2=v25["U"], I1=v25["I"], I2=v25["I"], cov=v25["V"]))
        rows.append(dict(scenario="B", endpoint="os", n=300, t1=15.0, t2=25.0, U1=v15["U"],
                         U2=v25["U"], I1=v15["I"], I2=v25["I"], cov=cv))
        scn = scenario_c("reset")
        v20 = rl.luo_var(scn, 500, 20.0, "pfs")
        rows.append(dict(scenario="C", endpoint="pfs", n=500, t1=20.0, t2=20.0, U1=v20["U"],
                         U2=v20["U"], I1=v20["I"], I2=v20["I"], cov=v20["V"]))
        write("alt_cov.csv", rows)

    if "alt_check" in which:
        # the influence-function derivation (used in the R code) against the
        # printed formulas of Luo et al. (2019)
        scn = scenario_b()
        rows = []
        for t, tp in ((25.0, 15.0), (25.0, 25.0)):
            luo = rl.luo_cov(scn, 300, t, tp, "os") if tp < t else rl.luo_var(scn, 300, t, "os")["V"]
            infl = rl.influence_cov(scn, 300, t, tp, "os")
            assert abs(luo - infl) < 1e-10, (luo, infl)
            rows.append(dict(scenario="B", t=t, t_prime=tp, luo=luo, influence=infl,
                             abs_diff=abs(luo - infl)))
        write("alt_check.csv", rows)

    if "gs" in which:
        rows = []
        cases = [
            ("a", dict(info_frac=[0.5, 0.75, 1.0])),
            ("b", dict(info_frac=[1 / 3, 2 / 3, 1.0], sf_eff="Pocock")),
            ("c", dict(info_frac=[0.5, 0.75, 1.0], futility="beta", sf_fut="HSD",
                       param_fut=-2, beta=0.1, mean=[0.8, 1.9, 2.8])),
            ("d", dict(info_frac=[0.5, 0.75, 1.0], futility="beta", sf_fut="HSD",
                       param_fut=-2, beta=0.1, mean=[0.8, 1.9, 2.8], binding=True)),
            ("e", dict(info_frac=[0.5, 0.75, 1.0], futility="fixed", fut_bound=[0.0, 0.5],
                       binding=True)),
        ]
        sd = np.array([1.05, 1.03, 1.02])
        rcor = np.array([[1.0, 0.80, 0.68], [0.80, 1.0, 0.85], [0.68, 0.85, 1.0]])
        sig_alt = np.outer(sd, sd) * rcor
        cases.append(("f", dict(info_frac=[0.5, 0.75, 1.0], futility="beta", sf_fut="OF",
                                beta=0.15, mean=[1.2, 2.0, 2.7], sigma_alt=sig_alt)))
        for label, kw in cases:
            lo, up = rg.bounds(**kw)
            for k in range(len(up)):
                rows.append(dict(case=label, look=k + 1, upper=up[k], lower=lo[k]))
            print(label, lo, up, flush=True)
        write("gs_bounds.csv", rows)
        # stopping probabilities for case c
        lo = np.array([r["lower"] for r in rows if r["case"] == "c"])
        up = np.array([r["upper"] for r in rows if r["case"] == "c"])
        t = np.array([0.5, 0.75, 1.0])
        corr = np.sqrt(np.minimum.outer(t, t) / np.maximum.outer(t, t))
        prow = []
        for label, mean in (("c_h1", [0.8, 1.9, 2.8]), ("c_h0", [0.0, 0.0, 0.0])):
            pu, pl = rg.stop_probs(lo, up, mean, corr)
            for k in range(3):
                prow.append(dict(case=label, look=k + 1, prob_upper=pu[k], prob_lower=pl[k]))
        write("gs_prob.csv", prow)
        srow = []
        for kind, par, alpha in (("OF", None, 0.025), ("Pocock", None, 0.025),
                                 ("HSD", -4.0, 0.025), ("HSD", 0.0, 0.025),
                                 ("power", 3.0, 0.025)):
            for tt in (0.25, 0.5, 0.75, 1.0):
                srow.append(dict(type=kind, param="" if par is None else par, alpha=alpha,
                                 t=tt, value=float(rg.spend(tt, alpha, kind, par))))
        write("gs_spend.csv", srow)

    if "power" in which:
        # ms_power reference for scenario B: OS looks at 50, 75, 100 percent of
        # 120 events, OF efficacy, non-binding beta spending HSD(-2), beta 0.1
        scn = scenario_b()
        n = 300
        targets = [60.0, 90.0, 120.0]
        times = [rl.event_time(scn, n, d, "os") for d in targets]
        mos = [rl.moments(scn, n, t, "os") for t in times]
        d = np.array([m["events"] for m in mos])
        info = np.array([m["info"] for m in mos])
        mu = np.array([m["z_direct"] for m in mos])
        corr0 = np.sqrt(np.minimum.outer(d, d) / np.maximum.outer(d, d))
        sig1 = np.sqrt(np.minimum.outer(info, info) / np.maximum.outer(info, info))
        tfrac = d / d[-1]
        lo, up = rg.bounds(tfrac, futility="beta", sf_fut="HSD", param_fut=-2, beta=0.1,
                           mean=mu, corr=corr0, sigma_alt=sig1)
        pu, pl = rg.stop_probs(lo, up, mu, sig1)
        # type I error under the null scenario at the same numbers of events
        s0 = rl.Scenario(scn.control, scn.control, 2.0, scn.accrual)
        times0 = [rl.event_time(s0, n, dd, "os") for dd in targets]
        v0 = np.array([rl.moments(s0, n, t, "os")["info"] for t in times0])
        c0 = np.sqrt(np.minimum.outer(v0, v0) / np.maximum.outer(v0, v0))
        pu0, _ = rg.stop_probs(np.full(3, -np.inf), up, np.zeros(3), c0)
        lo_b = lo.copy()
        pu0b, _ = rg.stop_probs(lo_b, up, np.zeros(3), c0)
        rows = []
        for k in range(3):
            rows.append(dict(look=k + 1, time=times[k], events=d[k], info=info[k],
                             z_mean=mu[k], upper=up[k], lower=lo[k], prob_upper=pu[k],
                             prob_lower=pl[k], time0=times0[k], prob_upper0_nb=pu0[k],
                             prob_upper0_b=pu0b[k]))
        write("power_B.csv", rows)
        print("power", pu.sum(), "type1", pu0.sum(), pu0b.sum(), flush=True)

    if "size" in which:
        rows = []
        za = stats.norm.isf(0.025)
        zb = stats.norm.isf(0.1)
        d_sch = (za + zb) ** 2 / (0.25 * np.log(0.7) ** 2)
        rows.append(dict(case="A_K1_schoenfeld", events=d_sch))
        # scenario A, two looks at 50 and 100 percent, OF, mean direct, n = 500:
        # solve the final calendar time for power 0.9
        scn = scenario_a()
        n = 500

        def power_at(tf):
            dfin = rl.events(scn, n, tf, "os")["events"]
            tm = [rl.event_time(scn, n, 0.5 * dfin, "os"), tf]
            mos = [rl.moments(scn, n, t, "os") for t in tm]
            dd = np.array([m["events"] for m in mos])
            info = np.array([m["info"] for m in mos])
            mu = np.array([m["z_direct"] for m in mos])
            c0 = np.sqrt(np.minimum.outer(dd, dd) / np.maximum.outer(dd, dd))
            s1 = np.sqrt(np.minimum.outer(info, info) / np.maximum.outer(info, info))
            lo, up = rg.bounds(dd / dd[-1], corr=c0)
            pu, _ = rg.stop_probs(lo, up, mu, s1)
            return pu.sum(), dfin

        tf = optimize.brentq(lambda x: power_at(x)[0] - 0.9, 20.0, 80.0, xtol=1e-9)
        rows.append(dict(case="A_K2_direct_n500", events=power_at(tf)[1]))
        rows.append(dict(case="A_K2_direct_n500_time", events=tf))
        write("size.csv", rows)

    print("total time", time.time() - t_start, flush=True)


if __name__ == "__main__":
    which = sys.argv[1:] or ["model", "events", "event_time", "moments", "alt",
                             "alt_check", "gs", "power", "size"]
    main(which)
