test_that(".check_look_test accepts lists and character vectors", {
  expect_null(.check_look_test(NULL, 2L))
  expect_identical(.check_look_test(list("pfs", c("pfs", "os")), 2L),
                   list("pfs", c("pfs", "os")))
  expect_identical(.check_look_test(c("pfs", "os"), 2L), list("pfs", "os"))
  expect_error(.check_look_test(list("pfs"), 2L), "one element per look")
  expect_error(.check_look_test(list("pfs", character(0)), 2L), "one element per look")
  expect_error(.check_look_test(list("pfs", NA_character_), 2L), "one element per look")
})
