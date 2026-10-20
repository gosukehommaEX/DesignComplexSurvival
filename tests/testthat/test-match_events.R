test_that(".match_events orders the final numbers of events by trigger", {
  expect_equal(.match_events(120, "os"), c(os = 120))
  expect_equal(.match_events(c(os = 220, pfs = 250), c("pfs", "os")), c(pfs = 250, os = 220))
  expect_error(.match_events(120, c("pfs", "os")), "named")
  expect_error(.match_events(c(os = 220), c("pfs", "os")), "PFS and OS")
  expect_error(.match_events(c(os = 220, pfs = 250), "os"), "OS")
})
