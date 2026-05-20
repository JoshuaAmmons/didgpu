# Tests for didgpu_equivalence() — pre-trends TOST equivalence (#99).
# Uses backend = "r" so it runs without a GPU or the reference package.

make_fit <- function(seed = 7L, placebo = 2L, reps = 60L, effects = 2L) {
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 12L,
                             tau_profile = c(0.5, 1.0), seed = seed)
  p$D <- as.integer(p$D >= 0.5)
  didgpu(p, "Y", "unit", "period", "D", effects = effects, placebo = placebo,
         bootstrap_reps = reps, backend = "r", seed = 1L, verbose = FALSE)
}

test_that("didgpu_equivalence returns the documented structure", {
  fit <- make_fit()
  eq <- didgpu_equivalence(fit, delta = 0.5)
  expect_s3_class(eq, "didgpu_equivalence")
  expect_true(all(c("event_time", "estimate", "std.error",
                    "equivalence_p", "passes_at_delta") %in% names(eq)))
  expect_equal(nrow(eq), nrow(fit$results$Placebos))
  expect_true(all(eq$event_time < 0))                 # placebos are pre-period
  expect_false(is.null(attr(eq, "delta")))
  expect_false(is.null(attr(eq, "breakdown_delta")))
  expect_true(is.finite(attr(eq, "breakdown_delta")))
})

test_that("TOST p-value matches the manual two-one-sided formula", {
  fit <- make_fit()
  delta <- 0.7
  eq <- didgpu_equivalence(fit, delta = delta)
  est <- as.numeric(fit$results$Placebos[, "Estimate"])
  se  <- as.numeric(fit$results$Placebos[, "SE"])
  man <- pmax(stats::pnorm((delta - est) / se, lower.tail = FALSE),
              stats::pnorm((delta + est) / se, lower.tail = FALSE))
  expect_equal(eq$equivalence_p, man, tolerance = 1e-12)
})

test_that("large delta -> joint PASS; tiny delta -> joint FAIL", {
  fit <- make_fit()
  big   <- didgpu_equivalence(fit, delta = 50)     # everything well within
  small <- didgpu_equivalence(fit, delta = 1e-4)   # nothing within
  expect_true(isTRUE(attr(big, "joint_pass")))
  expect_true(all(big$passes_at_delta))
  expect_false(isTRUE(attr(small, "joint_pass")))
})

test_that("breakdown_delta is the smallest margin that makes the joint test pass", {
  fit <- make_fit()
  bd <- attr(didgpu_equivalence(fit, delta = 1), "breakdown_delta")
  # Just above the breakdown -> joint pass; just below -> joint fail.
  expect_true(isTRUE(attr(didgpu_equivalence(fit, delta = bd * 1.001), "joint_pass")))
  expect_false(isTRUE(attr(didgpu_equivalence(fit, delta = bd * 0.999), "joint_pass")))
})

test_that("informative errors on misuse", {
  fit <- make_fit()
  expect_error(didgpu_equivalence(fit, delta = -1), "positive")
  expect_error(didgpu_equivalence(fit, delta = 0.5, alpha = 1.5), "alpha")
  expect_error(didgpu_equivalence(list(), delta = 0.5), "didgpu_result")

  # placebo > 0 but no bootstrap -> SEs are NA -> clear error.
  fit_nose <- make_fit(reps = 0L)
  expect_error(didgpu_equivalence(fit_nose, delta = 0.5), "bootstrap_reps")

  # no placebos at all -> clear error.
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                             tau_profile = c(0.5, 1.0), seed = 3L)
  p$D <- as.integer(p$D >= 0.5)
  fit_nopl <- didgpu(p, "Y", "unit", "period", "D", effects = 2L, placebo = 0L,
                     bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_error(didgpu_equivalence(fit_nopl, delta = 0.5), "placebo")
})

test_that("print method runs and returns invisibly", {
  fit <- make_fit()
  eq <- didgpu_equivalence(fit, delta = 0.5)
  expect_output(print(eq), "equivalence test")
  expect_identical(withVisible(print(eq))$visible, FALSE)
})
