test_that(".mtp_patterns gives the looks of the second endpoint for each procedure", {
  # first endpoint at analyses 1, 2, 4; second endpoint at analyses 2, 3, 4
  a1 <- c(1, 2, 4)
  a2 <- c(2, 3, 4)
  fut <- c(TRUE, FALSE, TRUE)
  elig <- function(p) lapply(p, `[[`, "eligible")
  type <- function(p) vapply(p, `[[`, character(1), "type")
  look <- function(p) vapply(p, `[[`, integer(1), "look")

  ov <- .mtp_patterns(a1, a2, fut, "fixed_sequence", "overall")
  expect_identical(type(ov), c("reject", "reject", "reject", "futility", "none"))
  expect_identical(look(ov), c(1L, 2L, 3L, 1L, 3L))
  expect_identical(elig(ov), list(1:3, 1:3, 3L, integer(0), integer(0)))

  sw <- .mtp_patterns(a1, a2, fut, "fixed_sequence", "stagewise")
  expect_identical(elig(sw), list(integer(0), 1L, 3L, integer(0), integer(0)))

  bf <- .mtp_patterns(a1, a2, fut, "bonferroni", "overall")
  expect_identical(elig(bf), list(1:3, 1:3, 1:3, integer(0), 1:3))

  bf2 <- .mtp_patterns(c(1, 3), c(1, 2, 3), c(TRUE, TRUE), "bonferroni", "overall")
  expect_identical(type(bf2), c("reject", "reject", "futility", "none"))
  expect_identical(elig(bf2)[[3L]], 1L)
})
