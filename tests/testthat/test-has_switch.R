test_that(".has_switch detects switching after progression", {
  m <- models_c()
  expect_true(.has_switch(m$control))
  expect_false(.has_switch(m$treatment))
  expect_false(.has_switch(ms_model(death = 0.1)))
  expect_true(.has_switch(ms_model(prog = 0.1, death = 0.01, death_prog = 0.05,
                                   switch_rate = 0.02)))
  expect_false(.has_switch(ms_model(death = 0.01, switch_pre = 0.02)))
})
