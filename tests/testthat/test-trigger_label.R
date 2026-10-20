test_that(".trigger_label describes the triggers", {
  expect_identical(.trigger_label("os"), "OS events")
  expect_identical(.trigger_label("calendar"), "calendar time")
  expect_identical(.trigger_label(c("pfs", "pfs", "os")), "PFS events and OS events")
  expect_identical(.trigger_label(c("calendar", "os")), "calendar time and OS events")
})
