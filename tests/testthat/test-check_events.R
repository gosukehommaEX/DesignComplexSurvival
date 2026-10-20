test_that(".check_events accepts one number or numbers named by endpoint", {
  expect_true(.check_events(120))
  expect_true(.check_events(c(pfs = 250, os = 220)))
  expect_error(.check_events(-1), "positive number")
  expect_error(.check_events(c(120, 130)), "positive number")
  expect_error(.check_events(c(pfs = 250, dfs = 220)), "named")
  expect_error(.check_events(c(pfs = 250, pfs = 220)), "named")
  expect_error(.check_events(c(pfs = NA, os = 220)), "named")
})
