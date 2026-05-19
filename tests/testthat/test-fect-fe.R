# fect_fe (two-way fixed effects counterfactual). Reference impl in R;
# CUDA kernel in src/cuda_fect_fe.cu compiles when nvcc is available.

# Build a simple panel with a known constant treatment effect.
build_fect_panel <- function(n_units = 80L, n_periods = 12L,
                              frac_treated = 0.5,
                              true_att = 1.0, seed = 7L) {
  set.seed(seed)
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units),
                         as.integer(n_units * frac_treated)))
  F_g[treated] <- sample(seq(3L, n_periods - 1L), length(treated),
                          replace = TRUE)
  panel <- expand.grid(unit = seq_len(n_units), period = seq_len(n_periods))
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g)
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] +
             true_att * panel$D +
             rnorm(nrow(panel), 0, 0.2)
  panel <- panel[order(panel$unit, panel$period),
                 c("unit", "period", "D", "Y")]
  attr(panel, "true_att") <- true_att
  panel
}

test_that("fect_fe recovers the true ATT within sampling error", {
  p <- build_fect_panel(n_units = 100L, n_periods = 15L,
                         true_att = 1.0, seed = 17L)
  fit <- didgpu_fect(p, "Y", "unit", "period", "D",
                      method = "fe", effects = 3L,
                      bootstrap_reps = 0L,
                      backend = "r", verbose = FALSE)
  ate_est <- as.numeric(fit$results$ATE[1, "Estimate"])
  expect_true(is.finite(ate_est))
  # With this panel size and noise, the ATE should be within ~0.15 of truth.
  expect_lt(abs(ate_est - 1.0), 0.15)
})

test_that("fect_fe produces a per-event-time effect breakdown", {
  p <- build_fect_panel(n_units = 100L, n_periods = 15L,
                         true_att = 1.0, seed = 17L)
  fit <- didgpu_fect(p, "Y", "unit", "period", "D",
                      method = "fe", effects = 4L,
                      bootstrap_reps = 0L,
                      backend = "r", verbose = FALSE)
  effects <- as.numeric(fit$results$Effects[, "Estimate"])
  expect_length(effects, 4L)
  # Each event-time should also be close to 1.0 (constant effect DGP).
  expect_true(all(abs(effects - 1.0) < 0.3))
})

test_that("fect_fe bootstrap produces finite SEs", {
  p <- build_fect_panel(n_units = 80L, n_periods = 12L,
                         true_att = 1.0, seed = 17L)
  fit <- suppressMessages(didgpu_fect(
    p, "Y", "unit", "period", "D",
    method = "fe", effects = 2L,
    bootstrap_reps = 10L, seed = 1L,
    backend = "r", verbose = FALSE))
  se <- as.numeric(fit$results$Effects[, "SE"])
  expect_true(all(is.finite(se)))
  expect_true(all(se > 0))
})

test_that("fect_fe checkpoint + resume reproduces the same result", {
  p <- build_fect_panel(n_units = 60L, n_periods = 10L,
                         true_att = 1.0, seed = 17L)
  cdir <- tempfile("fect_resume_")
  on.exit(unlink(cdir, recursive = TRUE), add = TRUE)
  fit1 <- suppressMessages(didgpu_fect(
    p, "Y", "unit", "period", "D",
    method = "fe", effects = 2L,
    bootstrap_reps = 3L, seed = 1L,
    checkpoint_dir = cdir,
    backend = "r", verbose = FALSE))
  fit2 <- suppressMessages(didgpu_fect(
    p, "Y", "unit", "period", "D",
    method = "fe", effects = 2L,
    bootstrap_reps = 3L, seed = 1L,
    checkpoint_dir = cdir,
    backend = "r", verbose = FALSE))
  expect_equal(as.numeric(fit1$results$Effects[, "Estimate"]),
               as.numeric(fit2$results$Effects[, "Estimate"]))
})

test_that("fect_fe returns the expected didgpu_fect_result class", {
  p <- build_fect_panel(n_units = 50L, n_periods = 8L, seed = 1L)
  fit <- didgpu_fect(p, "Y", "unit", "period", "D",
                      method = "fe", effects = 1L,
                      bootstrap_reps = 0L,
                      backend = "r", verbose = FALSE)
  expect_s3_class(fit, "didgpu_fect_result")
  expect_s3_class(fit, "didgpu_result")
  expect_equal(fit$method, "fe")
})

test_that("ife and mc now run successfully (no longer stubbed)", {
  p <- build_fect_panel(n_units = 50L, n_periods = 10L, seed = 1L)
  for (method in c("ife", "mc")) {
    fit <- didgpu_fect(p, "Y", "unit", "period", "D",
                        method = method, bootstrap_reps = 0L,
                        backend = "r", verbose = FALSE)
    expect_s3_class(fit, "didgpu_fect_result")
    expect_equal(fit$method, method)
  }
})
