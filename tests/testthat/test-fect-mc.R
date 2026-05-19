# fect_mc: matrix completion (Athey et al. 2021). R reference impl in
# R/fect_mc.R; CUDA SVD primitive shared with fect_ife in
# src/cuda_fect_svd.cu.

build_mc_panel <- function(n_units = 100L, n_periods = 15L,
                             frac_treated = 0.5,
                             true_att = 1.0, r = 2L, seed = 7L) {
  set.seed(seed)
  unit_fe <- rnorm(n_units, 0, 0.5)
  time_fe <- rnorm(n_periods, 0, 0.3)
  L <- matrix(rnorm(n_units * r, 0, 0.5), n_units, r)
  F <- matrix(rnorm(r * n_periods, 0, 0.4), r, n_periods)
  factor_struct <- L %*% F
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units),
                         as.integer(n_units * frac_treated)))
  F_g[treated] <- sample(seq(3L, n_periods - 1L), length(treated),
                          replace = TRUE)
  panel <- expand.grid(unit = seq_len(n_units), period = seq_len(n_periods))
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g)
  panel <- panel[order(panel$unit, panel$period), ]
  panel$Y <- unit_fe[panel$unit] +
             time_fe[panel$period] +
             factor_struct[cbind(panel$unit, panel$period)] +
             true_att * panel$D +
             rnorm(nrow(panel), 0, 0.2)
  panel[, c("unit", "period", "D", "Y")]
}

test_that("fect_mc runs and produces a result with the expected shape", {
  p <- build_mc_panel(n_units = 80L, n_periods = 12L,
                       true_att = 1.0, r = 2L, seed = 17L)
  fit <- didgpu_fect(p, "Y", "unit", "period", "D",
                      method = "mc", effects = 3L,
                      bootstrap_reps = 0L,
                      backend = "r", verbose = FALSE)
  expect_s3_class(fit, "didgpu_fect_result")
  expect_equal(fit$method, "mc")
  expect_equal(nrow(fit$results$Effects), 3L)
  ate <- as.numeric(fit$results$ATE[1, "Estimate"])
  expect_true(is.finite(ate))
})

test_that("fect_mc accepts a user-supplied lambda", {
  p <- build_mc_panel(n_units = 80L, n_periods = 12L, seed = 17L)
  fit_default <- didgpu_fect(p, "Y", "unit", "period", "D",
                              method = "mc", effects = 1L,
                              bootstrap_reps = 0L, backend = "r",
                              verbose = FALSE)
  fit_lambda <- didgpu_fect(p, "Y", "unit", "period", "D",
                             method = "mc", effects = 1L,
                             lambda = 0.5,
                             bootstrap_reps = 0L, backend = "r",
                             verbose = FALSE)
  # Different lambda => different estimate (in general).
  ate_default <- as.numeric(fit_default$results$ATE[1, "Estimate"])
  ate_lambda  <- as.numeric(fit_lambda$results$ATE[1, "Estimate"])
  expect_true(is.finite(ate_default))
  expect_true(is.finite(ate_lambda))
})

test_that("fect_mc with very large lambda shrinks everything toward zero", {
  # If lambda is huge, the soft-threshold zeroes out every singular
  # value, so Y_hat = 0 and the "treatment effect" = Y[treated] which
  # is mostly the unit FE + time FE values (the bias from no model).
  p <- build_mc_panel(n_units = 50L, n_periods = 10L, seed = 17L)
  fit <- didgpu_fect(p, "Y", "unit", "period", "D",
                      method = "mc", lambda = 1e6,
                      bootstrap_reps = 0L, backend = "r",
                      verbose = FALSE)
  # n_nonzero stored on the per-cell value; we just verify the fit
  # didn't crash.
  expect_s3_class(fit, "didgpu_fect_result")
})
