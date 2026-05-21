test_that("only_never_switchers matches reference bit-for-bit", {
  skip_if_no_reference()
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 15L,
                              frac_treated = 0.6,
                              min_treat_period = 5L, max_treat_period = 10L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              sigma = 0.4, seed = 17L)

  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 1, graph_off = TRUE,
      only_never_switchers = TRUE
    )
  ))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 1L, only_never_switchers = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)

  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))), 1e-10)
  expect_lt(max(abs(as.numeric(us$results$Placebos[, "Estimate"]) -
                    as.numeric(ref$results$Placebos[, 1]))), 1e-10)
})

test_that("same_switchers matches reference bit-for-bit", {
  skip_if_no_reference()
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 15L,
                              frac_treated = 0.6,
                              min_treat_period = 5L, max_treat_period = 10L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              sigma = 0.4, seed = 17L)

  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE, same_switchers = TRUE
    )
  ))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, same_switchers = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))), 1e-10)
})

test_that("only_never_switchers shifts estimates relative to default", {
  # On panels where switchers' pre-switch rows materially differ from
  # never-treated controls, the two should give different answers.
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 15L,
                              frac_treated = 0.6,
                              min_treat_period = 5L, max_treat_period = 10L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              sigma = 0.4, seed = 17L)
  fit_default <- didgpu(p, "Y", "unit", "period", "D",
                         effects = 3L, bootstrap_reps = 0L,
                         backend = "r", verbose = FALSE)
  fit_only_never <- didgpu(p, "Y", "unit", "period", "D",
                            effects = 3L, only_never_switchers = TRUE,
                            bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  # At least one event-time estimate should differ.
  expect_true(any(abs(fit_default$results$Effects[, "Estimate"] -
                      fit_only_never$results$Effects[, "Estimate"]) > 1e-6))
})
