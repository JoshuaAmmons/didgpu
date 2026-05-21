# Phase-1 wiring test for task #81.
#
# The CUDA testmechs bootstrap (cuda_testmechs_bootstrap.cu) uses
# cuRAND for the multinomial-resample step, whereas the R fallback
# uses base R's Mersenne-Twister. They produce DIFFERENT random
# sequences at the same seed, so per-replicate output differs — but
# both are valid bootstrap procedures and converge to the same
# population moments.
#
# These tests verify the wiring contract (shape, normalisation) and
# a Monte-Carlo equivalence on the bootstrap means.

test_that(".testmechs_bootstrap_cuda returns the right shape", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  set.seed(7L)
  n <- 200L; K <- 2L; dy <- 3L
  d <- sample.int(2L, n, replace = TRUE) - 1L      # {0, 1}
  m <- sample.int(K,  n, replace = TRUE)            # {1, .., K}
  y <- sample.int(dy, n, replace = TRUE)            # {1, .., dy}
  B <- 64L

  M <- didgpu:::.testmechs_bootstrap_cuda(d, m, y, B,
                                           method = "nonparametric",
                                           seed   = 17L)
  expect_equal(dim(M), c(B, 2L * K * dy))
  # Per-row, each per-D block (first K*dy and last K*dy entries) sums
  # to either 1 (saw at least one d=that) or 0 (saw none — possible
  # for small n with one D heavily represented). Tolerate 0 + 1.
  row_sums_d0 <- rowSums(M[, seq_len(K * dy), drop = FALSE])
  row_sums_d1 <- rowSums(M[, (K * dy + 1L):(2L * K * dy), drop = FALSE])
  expect_true(all(abs(row_sums_d0 - 1) < 1e-10 | abs(row_sums_d0) < 1e-10))
  expect_true(all(abs(row_sums_d1 - 1) < 1e-10 | abs(row_sums_d1) < 1e-10))
})

test_that(".testmechs_bootstrap_cuda falls back for method = 'bayes'", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  set.seed(8L)
  n <- 100L; K <- 2L; dy <- 2L
  d <- sample.int(2L, n, replace = TRUE) - 1L
  m <- sample.int(K,  n, replace = TRUE)
  y <- sample.int(dy, n, replace = TRUE)

  # Bayes path goes to R; verify the result matches the R-direct call.
  set.seed(99L); M_cuda_for_bayes <- didgpu:::.testmechs_bootstrap_cuda(
    d, m, y, 32L, method = "bayes", seed = 17L)
  set.seed(99L); M_r_direct <- didgpu:::.testmechs_bootstrap_r(
    d, m, y, 32L, method = "bayes", seed = 17L)
  expect_equal(M_cuda_for_bayes, M_r_direct)
})

test_that("CUDA bootstrap means agree with R bootstrap means within Monte-Carlo error", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  # Large-B Monte-Carlo equivalence. cuRAND and MT19937 give different
  # per-replicate values, but the empirical mean over B replicates is
  # an unbiased estimator of the same population proportion.
  set.seed(11L)
  n <- 500L; K <- 2L; dy <- 3L
  d <- sample.int(2L, n, replace = TRUE) - 1L
  m <- sample.int(K,  n, replace = TRUE)
  y <- sample.int(dy, n, replace = TRUE)
  B <- 1000L

  M_cuda <- didgpu:::.testmechs_bootstrap_cuda(d, m, y, B,
                                                method = "nonparametric",
                                                seed   = 17L)
  M_r <- didgpu:::.testmechs_bootstrap_r(d, m, y, B,
                                          method = "nonparametric",
                                          seed   = 17L)
  # Population values: count(d, m, y) / count(d) over the original sample.
  # Both bootstraps should average to those same proportions ± Monte-Carlo
  # SD of about sqrt(p(1-p)/B). For B=1000 and p ~ 0.1, SD ~ 0.01.
  # Allow 0.05 to be generous (3-5 SDs).
  diff_means <- colMeans(M_cuda) - colMeans(M_r)
  expect_true(max(abs(diff_means)) < 0.05,
              info = paste("max abs diff =", max(abs(diff_means))))
})
