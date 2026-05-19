# Placebo + equivalence tests for fect estimators.

test_that("didgpu_fect_placebo returns a data.frame with expected columns", {
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  pl <- didgpu_fect_placebo(p, "Y", "unit", "period", "D",
                              method = "fe", n_placebos = 2L)
  expect_s3_class(pl, "didgpu_fect_placebo")
  expect_s3_class(pl, "data.frame")
  expect_true(all(c("horizon", "estimate", "se",
                     "p_value", "n_cells") %in% names(pl)))
  expect_equal(pl$horizon, c(-1L, -2L))
})

test_that("didgpu_fect_placebo works for all three methods", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  for (method in c("fe", "ife", "mc")) {
    pl <- didgpu_fect_placebo(p, "Y", "unit", "period", "D",
                                method = method, n_placebos = 1L)
    expect_s3_class(pl, "didgpu_fect_placebo")
    expect_true(is.finite(pl$estimate[1]))
  }
})

test_that("didgpu_fect_placebo recovers ~0 effect on a clean panel (no anticipation)", {
  # With a panel that satisfies parallel trends + no anticipation,
  # the placebo effect should be small (much less than the post-treatment
  # effect of 1.0).
  p <- didgpu_simulate_panel(n_units = 200L, n_periods = 15L,
                              tau_profile = c(1.0, 1.0, 1.0),
                              sigma = 0.2, seed = 17L)
  pl <- didgpu_fect_placebo(p, "Y", "unit", "period", "D",
                              method = "fe", n_placebos = 2L)
  # Placebo estimates should be smaller in magnitude than the true ATT.
  expect_true(all(abs(pl$estimate) < 0.5))
})

test_that("didgpu_fect_placebo with bootstrap returns finite SEs", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  pl <- didgpu_fect_placebo(p, "Y", "unit", "period", "D",
                              method = "fe", n_placebos = 2L,
                              bootstrap_reps = 10L, seed = 1L)
  expect_true(all(is.finite(pl$se)))
  expect_true(all(pl$se > 0))
})

test_that("didgpu_fect_equivalence augments placebo result with eq columns", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  pl <- didgpu_fect_placebo(p, "Y", "unit", "period", "D",
                              method = "fe", n_placebos = 2L,
                              bootstrap_reps = 10L, seed = 1L)
  eq <- didgpu_fect_equivalence(pl, delta = 1.0)
  expect_true("equivalence_p" %in% names(eq))
  expect_true("passes_at_delta" %in% names(eq))
  # Large delta = generous tolerance = should pass.
  expect_true(any(eq$passes_at_delta))
})

test_that("didgpu_fect_equivalence with tiny delta fails (correctly strict)", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  pl <- didgpu_fect_placebo(p, "Y", "unit", "period", "D",
                              method = "fe", n_placebos = 2L,
                              bootstrap_reps = 10L, seed = 1L)
  eq_tight <- didgpu_fect_equivalence(pl, delta = 1e-6)
  # With near-zero tolerance, no horizon should pass.
  expect_false(any(eq_tight$passes_at_delta))
})

test_that("didgpu_fect_equivalence errors on non-placebo input", {
  expect_error(
    didgpu_fect_equivalence(data.frame(x = 1), delta = 1),
    "didgpu_fect_placebo"
  )
})

# -------- CV-based lambda for fect_mc --------

test_that("fect_mc CV picks a positive lambda on a clean panel", {
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0),
                              seed = 17L)
  mats <- didgpu:::.fect_build_matrices(p, "Y", "unit", "period", "D")
  l <- didgpu:::.fect_mc_cv_lambda(mats$Y, mats$M, K = 3L, n_grid = 4L)
  expect_true(is.finite(l))
  expect_gt(l, 0)
})

test_that("fect_mc CV path runs end-to-end via didgpu_fect", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                              tau_profile = c(0.5, 1.0),
                              seed = 17L)
  # No lambda passed = CV is invoked on iter 0.
  fit <- didgpu_fect(p, "Y", "unit", "period", "D",
                      method = "mc", effects = 1L,
                      bootstrap_reps = 0L,
                      backend = "r", verbose = FALSE)
  expect_s3_class(fit, "didgpu_fect_result")
})
