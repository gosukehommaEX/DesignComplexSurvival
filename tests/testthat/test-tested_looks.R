test_that(".tested_looks gives the looks of each endpoint", {
  l <- ms_looks(c("pfs", "pfs", "os"), fraction = c(0.6, 1, 1),
                test = list("pfs", c("pfs", "os"), "os"))
  expect_identical(.tested_looks(l, c("pfs", "os")), list(pfs = 1:2, os = 2:3))
  expect_identical(.tested_looks(l, c("os", "pfs")), list(os = 2:3, pfs = 1:2))
  l0 <- ms_looks("os", fraction = c(0.5, 1))
  expect_identical(.tested_looks(l0, "os"), list(os = 1:2))
  expect_identical(.tested_looks(l0, c("pfs", "os")), list(pfs = 1:2, os = 1:2))
  expect_error(.tested_looks(l, "os"), "PFS")
  expect_error(.tested_looks(ms_looks("os", fraction = c(0.5, 1), test = c("os", "os")),
                             c("pfs", "os")), "PFS is not tested")
  l11 <- ms_looks("os", fraction = (1:11) / 11)
  expect_error(.tested_looks(l11, c("pfs", "os")), "20 times")
})
