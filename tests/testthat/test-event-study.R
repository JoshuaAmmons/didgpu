test_that("didgpu_event_study_data returns the right shape", {
  p <- didgpu_simulate_panel(n_units = 50L, n_periods = 14L, seed = 5L,
                              min_treat_period = 5L, max_treat_period = 9L,
                              tau_profile = c(0.5, 1.0, 1.2))
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 3L, placebo = 2L,
                 bootstrap_reps = 5L, seed = 1L,
                 backend = "r", verbose = FALSE)
  ev <- didgpu_event_study_data(fit)
  expect_s3_class(ev, "data.frame")
  expect_equal(names(ev),
               c("event_time", "estimate", "std.error",
                 "conf.low", "conf.high", "kind"))
  # Placebos at -2, -1; effects at 1, 2, 3.
  expect_equal(ev$event_time, c(-2, -1, 1, 2, 3))
  expect_equal(ev$kind,
               c("placebo", "placebo", "effect", "effect", "effect"))
})

test_that("event_study_data is robust to fit without placebos", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 5L,
                              min_treat_period = 3L, max_treat_period = 5L)
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 2L, placebo = 0L,
                 bootstrap_reps = 3L, seed = 1L,
                 backend = "r", verbose = FALSE)
  ev <- didgpu_event_study_data(fit)
  expect_equal(nrow(ev), 2L)
  expect_true(all(ev$kind == "effect"))
})
