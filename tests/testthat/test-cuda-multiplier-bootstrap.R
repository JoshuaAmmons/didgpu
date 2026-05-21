# Phase-2 #83 tests: GPU multiplier (wild) bootstrap.

test_that("didgpu_cuda_multiplier_bootstrap_r returns the right shape", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  set.seed(13L)
  n_units <- 40L; n_dims <- 6L
  IF <- matrix(rnorm(n_units * n_dims), nrow = n_units, ncol = n_dims)
  out <- didgpu:::didgpu_cuda_multiplier_bootstrap_r(
    IF = IF, B = 32L, mult_kind = 0L, seed = 17L)
  expect_equal(dim(out), c(32L, n_dims))
  expect_true(all(is.finite(out)))
})

test_that("Rademacher weights produce zero-mean bootstrap deviations on average", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  # E[xi] = 0 under Rademacher, so the cross-replicate MEAN of
  # bootstrap deviations should be ~0 for any fixed IF column.
  set.seed(21L)
  n_units <- 50L; n_dims <- 3L
  IF <- matrix(rnorm(n_units * n_dims), nrow = n_units, ncol = n_dims)
  B <- 2000L
  out <- didgpu:::didgpu_cuda_multiplier_bootstrap_r(
    IF = IF, B = B, mult_kind = 0L, seed = 17L)
  # Population SD of column j is sd(IF[, j]) * sqrt(n_units). MC SE
  # of column-mean estimator is sd / sqrt(B) ~ 0.1 * sd / 45 ~ small.
  col_means <- colMeans(out)
  col_sds   <- apply(out, 2L, stats::sd)
  expect_true(all(abs(col_means) < 0.3 * col_sds),
              info = paste("col_means:", paste(round(col_means, 3), collapse = ", ")))
})

test_that("N(0,1) variant returns finite values", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  set.seed(27L)
  IF <- matrix(rnorm(60), nrow = 20, ncol = 3)
  out <- didgpu:::didgpu_cuda_multiplier_bootstrap_r(
    IF = IF, B = 64L, mult_kind = 1L, seed = 17L)
  expect_equal(dim(out), c(64L, 3L))
  expect_true(all(is.finite(out)))
})

test_that("GPU multiplier-bootstrap SDs match a CPU IF-shortcut reference within MC error", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  set.seed(33L)
  n_units <- 60L; n_dims <- 4L
  IF <- matrix(rnorm(n_units * n_dims), nrow = n_units, ncol = n_dims)
  B <- 5000L

  gpu_mat <- didgpu:::didgpu_cuda_multiplier_bootstrap_r(
    IF = IF, B = B, mult_kind = 0L, seed = 17L)
  gpu_sds <- apply(gpu_mat, 2L, stats::sd)

  set.seed(17L, kind = "Mersenne-Twister")
  cpu_mat <- matrix(NA_real_, nrow = B, ncol = n_dims)
  for (b in seq_len(B)) {
    xi <- sample(c(-1, 1), n_units, replace = TRUE)
    cpu_mat[b, ] <- as.numeric(crossprod(xi, IF))
  }
  cpu_sds <- apply(cpu_mat, 2L, stats::sd)

  rel_diff <- abs(gpu_sds - cpu_sds) / cpu_sds
  expect_true(max(rel_diff) < 0.05,
              info = paste("max rel diff =", round(max(rel_diff), 4)))
})

test_that("didgpu_cs(backend='cuda', bootstrap_kind='multiplier') populates SE", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  p <- didgpu_simulate_panel(n_units = 30L, n_periods = 6L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  p$D <- as.integer(p$D >= 0.5)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR", aggregation = "event",
                    bootstrap_reps = 200L,
                    bootstrap_kind = "multiplier",
                    backend = "cuda", verbose = FALSE)
  expect_true("se" %in% names(fit$att_gt))
  expect_true(any(!is.na(fit$att_gt$se)))
})
