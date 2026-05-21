# didgpu_fect: scaffolded counterfactual-prediction estimators.
# Tests verify the public API surface exists and errors gracefully
# until the real implementation lands.

test_that("didgpu_fect dispatches all three methods successfully", {
  p <- didgpu_simulate_panel(n_units = 50L, n_periods = 10L, seed = 1L)
  for (method in c("fe", "ife", "mc")) {
    fit <- didgpu_fect(p, "Y", "unit", "period", "D",
                        method = method,
                        bootstrap_reps = 0L,
                        backend = "r", verbose = FALSE)
    expect_s3_class(fit, "didgpu_fect_result")
    expect_equal(fit$method, method)
  }
})

test_that("didgpu_fect validates required args before dispatching", {
  p <- didgpu_simulate_panel(n_units = 30L, n_periods = 8L, seed = 1L)
  expect_error(
    didgpu_fect(p, "Ynope", "unit", "period", "D"),
    "column not in df: Ynope"
  )
  expect_error(
    didgpu_fect(p, "Y", "unit", "period", "D", method = "wrong"),
    "should be one of"
  )
})

test_that("didgpu_fect default method is 'fe' (which is implemented, not stubbed)", {
  # The default method 'fe' is now actually implemented, so the call
  # succeeds rather than erroring. ife and mc still stub out.
  p <- didgpu_simulate_panel(n_units = 30L, n_periods = 8L, seed = 1L)
  fit <- didgpu_fect(p, "Y", "unit", "period", "D",
                      bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_s3_class(fit, "didgpu_fect_result")
  expect_equal(fit$method, "fe")
})
