test_that("print.ms_power shows the overview and the table", {
  res <- ms_power(scn_a(), n = 400, looks = ms_looks("os", fraction = c(0.5, 1)),
                  test = ms_test(), events = 300)
  expect_output(print(res), "Power:")
  expect_output(print(res), "Type I error")
  expect_output(print(res), "z_upper")
  expect_output(expect_invisible(print(res)), "Power:")
})

test_that("print.ms_power shows a design with two endpoints", {
  res <- run_design("D2")
  expect_output(print(res), "two endpoints")
  expect_output(print(res), "Bonferroni")
  expect_output(print(res), "FWER under global null")
  expect_output(print(res), "null_scenarios\\[\\[1\\]\\] \\(PFS null\\)")
  expect_output(expect_invisible(print(res)), "at least one")
})
