test_that(".check_calendar accepts OS alone and PFS alone", {
  s <- scn_cal()
  expect_true(.check_calendar(s, "os", "direct", "studentized", c("os", "calendar")))
  expect_true(.check_calendar(s, "pfs", "schoenfeld", "studentized", "pfs"))
  expect_true(.check_calendar(s, "pfs", "direct", "null", "os"))
  expect_true(.check_calendar(scn_c(), c("os", "pfs"), "schoenfeld", "studentized", "pfs"))
})

test_that(".check_calendar rejects what the calendar-time computations do not cover", {
  s <- scn_cal()
  expect_error(.check_calendar(s, c("pfs", "os")), "one endpoint only")
  expect_error(.check_calendar(s, "os", "schoenfeld"), "mean = \"direct\"")
  expect_error(.check_calendar(s, "os", "direct", "null", c("pfs", "os")),
               "cannot be driven by PFS events")
  expect_error(.check_calendar(s, "pfs", "direct", "studentized", "os"),
               "needs variance = \"null\"")
})
