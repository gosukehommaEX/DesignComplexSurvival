test_that(".check_trigger recycles and checks the trigger", {
  expect_null(.check_trigger(NULL, 3L))
  expect_identical(.check_trigger("os", 3L), rep("os", 3))
  expect_identical(.check_trigger(c("pfs", NA, "os"), 3L), c("pfs", NA, "os"))
  expect_identical(.check_trigger(NA, 2L), c(NA_character_, NA_character_))
  expect_error(.check_trigger("dfs", 2L), "trigger")
  expect_error(.check_trigger(c("os", "pfs"), 3L), "trigger")
  expect_error(.check_trigger(1, 1L), "trigger")
})
