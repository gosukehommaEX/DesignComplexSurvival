# Expected values: inst/validation/python/results/model_states.csv
# (forward clock by products of matrix exponentials, reset clock by adaptive
# quadrature of the convolution over the progression time)

test_that("ms_prob matches the reference values for the reset clock", {
  m <- models_c("reset")
  p <- ms_prob(m$control, c(2, 5, 10, 20))
  ref <- rbind(
    c(0.8025187979624785, 0.11462994707077737, 0.057023865257011544, 0.0258273897097326),
    c(0.5886049696783552, 0.20366798276119408, 0.1263854427451381, 0.08134160481531263),
    c(0.3753110988513995, 0.22921654229747046, 0.1998563971004928, 0.19561596175063728),
    c(0.1525901057568839, 0.16808042268412016, 0.25874657626960607, 0.42058289528938986))
  got <- as.matrix(p[, c("progression_free", "progressed", "switched", "dead")])
  expect_lt(max(abs(got - ref)), 1e-10)
  p_t <- ms_prob(m$treatment, c(2, 20))
  expect_equal(p_t$progressed, c(0.10718752899784921, 0.3470709458139385), tolerance = 1e-10)
  expect_equal(p_t$switched, c(0, 0))
})

test_that("ms_prob matches the reference values for the forward clock", {
  m <- models_c("forward")
  p <- ms_prob(m$control, c(2, 5, 10, 20))
  ref <- rbind(
    c(0.8025187979624785, 0.11462994707077737, 0.05702386525701154, 0.025827389709732584),
    c(0.5886049696783551, 0.20118518315021522, 0.12556939829355684, 0.08464044887787266),
    c(0.3753110988513994, 0.2245281562411191, 0.19744321918484906, 0.2027175257226321),
    c(0.15259010575688384, 0.16391615594278708, 0.2545609888742223, 0.42893274942610654))
  got <- as.matrix(p[, c("progression_free", "progressed", "switched", "dead")])
  expect_lt(max(abs(got - ref)), 1e-10)
})

test_that("both clocks agree when the post-progression intensities are constant", {
  args <- list(time = c(0, 5), prog = c(0.1, 0.05), death = 0.01,
               switch_prob = c(0.2, 0.4), death_prog = 0.06,
               switch_rate = 0.03, death_switch = 0.02, dropout = 0)
  m_r <- do.call(ms_model, c(args, clock = "reset"))
  m_f <- do.call(ms_model, c(args, clock = "forward"))
  tt <- c(1, 4, 8, 15, 30)
  expect_lt(max(abs(as.matrix(ms_prob(m_r, tt)) - as.matrix(ms_prob(m_f, tt)))), 1e-12)
  expect_lt(max(abs(rowSums(ms_prob(m_r, tt)[, -1]) - 1)), 1e-12)
})

test_that("ms_prob without progression is exponential", {
  p <- ms_prob(ms_model(death = 0.1), c(0, 5))
  expect_equal(p$progression_free, exp(-0.1 * c(0, 5)), tolerance = 1e-14)
  expect_equal(p$dead, 1 - exp(-0.1 * c(0, 5)), tolerance = 1e-14)
  expect_error(ms_prob(list(), 1), "ms_model")
})
