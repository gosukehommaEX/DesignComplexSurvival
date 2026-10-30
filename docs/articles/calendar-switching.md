# Switching from a calendar time

``` r

library(DesignComplexSurvival)
```

Switching to the experimental treatment is sometimes allowed only from a
calendar time, for example after an interim analysis has shown a
benefit. With `switch_start` in
[`ms_scenario()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_scenario.md),
nobody switches before that calendar time (time since the start of
accrual). A patient enrolled at time e can switch from the follow-up
time max(`switch_start` - e, 0) on: a progression before that time does
not lead to a switch at progression, and a patient who progressed before
it can switch with `switch_rate` from that time on.

## Scenario

The models are those of
[`vignette("switching")`](https://gosukehommaEX.github.io/DesignComplexSurvival/articles/switching.md):
in the control group 30 percent of the patients switch at progression
and the others may switch later with a hazard of 0.04 per month; the
death hazard after progression is 0.06 without switch and 0.035 after a
switch. The post-progression intensities are constant, so that the reset
and forward clocks give the same model; the forward clock is used here
because its computations with `switch_start` are faster.

``` r

ctl <- ms_model(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01,
                switch_prob = 0.3, death_prog = 0.06, switch_rate = 0.04,
                death_switch = 0.035, dropout = 0.003, clock = "forward")
trt <- ms_model(time = c(0, 4), prog = c(0.06, 0.05), death = 0.008,
                death_prog = 0.055, dropout = 0.003, clock = "forward")
acc <- ms_accrual(rate = 25)
looks <- ms_looks("os", fraction = c(0.6, 1))
```

## Power by the start of switching

The design enrolls 500 patients and analyses OS at 60 percent of the
final 300 deaths and at 300 deaths. The power and the type I error for
switching from randomization (`switch_start = 0`) and from months 6 to
24:

``` r

starts <- c(0, 6, 12, 18, 24)
res <- lapply(starts, function(s)
  ms_power(ms_scenario(ctl, trt, accrual = acc, switch_start = s), n = 500,
           looks = looks, test = ms_test(), events = 300))
tab <- data.frame(switch_start = starts,
                  power = sapply(res, function(x) x$power),
                  type1 = sapply(res, function(x) x$type1[["nonbinding"]]),
                  interim = sapply(res, function(x) x$looks$time[1]),
                  final = sapply(res, function(x) x$study_duration))
knitr::kable(tab, digits = c(0, 4, 5, 1, 1))
```

| switch_start |  power |   type1 | interim | final |
|-------------:|-------:|--------:|--------:|------:|
|            0 | 0.2206 | 0.02500 |    28.8 |  43.9 |
|            6 | 0.2507 | 0.02500 |    28.7 |  43.8 |
|           12 | 0.3446 | 0.02502 |    28.5 |  43.5 |
|           18 | 0.5033 | 0.02504 |    28.0 |  42.9 |
|           24 | 0.6736 | 0.02504 |    27.7 |  42.2 |

``` r

plot(tab$switch_start, tab$power, type = "b", pch = 19,
     xlab = "Calendar month from which switching is possible", ylab = "Power")
```

![Power of the OS test by the calendar time from which switching is
possible](calendar-switching_files/figure-html/power-plot-1.png)

When switching starts later, fewer control patients have switched by the
time of the analyses. The power is 0.221 with switching from
randomization and 0.674 with switching from month 24, and the interim
analysis takes place at month 28.8 and month 27.7, respectively.

## Type I error

The distribution of OS now depends on the entry time. Under the null
hypothesis
([`ms_null()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_null.md),
the control model with its switching in both groups) the switching still
changes the death hazard, so that the score no longer has independent
increments. The correlation between the looks is computed from the
per-patient influence functions, integrated over the entry time (see
[`?lr_cov`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/lr_cov.md)).
The efficacy bounds still come from the event fractions, so that the
type I error computed with the exact correlation can differ from 0.025:
it ranges from 0.025 to 0.02504 in the table above.

## Scope

The computations with `switch_start` cover the OS statistic alone with
`mean = "direct"`. PFS is not affected by the switching. Two endpoints
tested together, an OS test at looks driven by PFS events, switching
before progression and the Schoenfeld approximation of the mean give an
error:

``` r

scn12 <- ms_scenario(ctl, trt, accrual = acc, switch_start = 12)
err <- function(expr) tryCatch({
  expr
  "no error"
}, error = conditionMessage)
err(ms_power(scn12, n = 500, looks = looks, test = ms_test(mean = "schoenfeld"),
             events = 300))
#> [1] "Switching from a calendar time ('switch_start') needs mean = \"direct\" for OS."
err(ms_power(scn12, n = 500, looks = ms_looks("pfs", fraction = c(0.6, 1)),
             test = ms_test(), events = 300))
#> [1] "With switching from a calendar time ('switch_start'), the looks of an OS test cannot be driven by PFS events."
```

[`ms_size()`](https://gosukehommaEX.github.io/DesignComplexSurvival/reference/ms_size.md)
solves for the number of events, the sample size, the accrual duration
or the follow-up with `switch_start` in the same way as without it.
