test_that(".plot_survival draws OS of the groups when PFS is OS", {
  x <- plot_design_b()
  d <- with_pdf(.plot_survival(x, NULL, list()))
  expect_gt(attr(d, "pdf_size"), 0)
  expect_identical(unique(d$endpoint), "os")
  expect_identical(unique(d$group), c("control", "treatment"))
  tt <- seq(0, x$study_duration, length.out = 201)
  for (g in c("control", "treatment")) {
    s <- d[d$group == g, ]
    expect_identical(s$time, tt)
    expect_identical(s$surv, ms_surv(x$scenario[[g]], tt, "os")$surv)
  }
  expect_identical(unique(d$label), c("Control, OS", "Treatment, OS"))
  d <- with_pdf(.plot_survival(x, 48, list(main = "OS")))
  expect_equal(max(d$time), 48)
  expect_error(.plot_survival(x, -1, list()), "positive number")
  expect_error(.plot_survival(x, c(1, 2), list()), "positive number")
})

test_that(".plot_survival draws PFS and OS with switching from a calendar time", {
  x <- ms_power(scn_cal(), n = 500, looks = ms_looks("calendar", time = c(24, 36)),
                test = ms_test())
  d <- with_pdf(.plot_survival(x, 30, list()))
  expect_identical(unique(d$endpoint), c("os", "pfs"))
  tt <- seq(0, 30, length.out = 201)
  for (g in c("control", "treatment")) {
    for (e in c("os", "pfs")) {
      s <- d[d$group == g & d$endpoint == e, ]
      expect_identical(s$surv, ms_surv(x$scenario[[g]], tt, e)$surv)
    }
  }
  expect_identical(unique(d$label),
                   c("Control, OS (entry after the start of switching)",
                     "Treatment, OS (entry after the start of switching)",
                     "Control, PFS", "Treatment, PFS"))
})
