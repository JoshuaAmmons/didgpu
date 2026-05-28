# Edge-case robustness tests for the GPU bootstrap kernels.
#
# Degenerate shapes (single cluster, single dim, B=1, all-one-cluster)
# are exactly where off-by-one indexing or empty-reduction bugs hide.
# These call the kernels directly so a regression surfaces immediately.

test_that("cluster bootstrap: single cluster (n_clusters = 1) is well-defined", {
  skip_if_no_cuda()
  # With one cluster, every replicate picks it n_clusters = 1 times, so
  # weight is always 1 and every replicate equals colSums(IF).
  IF <- matrix(c(1.0, 2.0, 3.0, -1.0, 0.5, 0.25), nrow = 3, byrow = TRUE)
  out <- didgpu:::didgpu_cuda_cluster_bootstrap_r(
    IF = IF, cluster_id = c(0L, 0L, 0L), n_clusters = 1L,
    B = 8L, seed = 17L)
  expect_equal(dim(out), c(8L, 2L))
  expected_col <- colSums(IF)
  for (b in seq_len(8L)) {
    expect_equal(as.numeric(out[b, ]), expected_col, tolerance = 1e-12)
  }
})

test_that("cluster bootstrap: single dim (n_dims = 1)", {
  skip_if_no_cuda()
  IF <- matrix(rnorm(20), ncol = 1)
  out <- didgpu:::didgpu_cuda_cluster_bootstrap_r(
    IF = IF, cluster_id = as.integer(0:19), n_clusters = 20L,
    B = 16L, seed = 3L)
  expect_equal(dim(out), c(16L, 1L))
  expect_true(all(is.finite(out)))
})

test_that("cluster bootstrap: B = 1 single replicate", {
  skip_if_no_cuda()
  IF <- matrix(rnorm(12), nrow = 4)
  out <- didgpu:::didgpu_cuda_cluster_bootstrap_r(
    IF = IF, cluster_id = as.integer(0:3), n_clusters = 4L,
    B = 1L, seed = 9L)
  expect_equal(dim(out), c(1L, 3L))
  expect_true(all(is.finite(out)))
})

test_that("multiplier bootstrap: single dim and B = 1", {
  skip_if_no_cuda()
  IF <- matrix(rnorm(30), ncol = 1)
  out1 <- didgpu:::didgpu_cuda_multiplier_bootstrap_r(
    IF = IF, B = 1L, mult_kind = 0L, seed = 5L)
  expect_equal(dim(out1), c(1L, 1L))
  # The single-rep Rademacher sum is sum(+/- IF), so |out| <= sum|IF|.
  expect_lte(abs(out1[1, 1]), sum(abs(IF)) + 1e-9)
})

test_that("multiplier bootstrap: normal weights have ~unit variance scaling", {
  skip_if_no_cuda()
  # For N(0,1) multipliers, Var(sum_i xi_i * IF_i) = sum_i IF_i^2.
  # Check the empirical column variance matches sum(IF^2) within MC error.
  set.seed(1L)
  IF <- matrix(rnorm(50), ncol = 1)
  B <- 20000L
  out <- didgpu:::didgpu_cuda_multiplier_bootstrap_r(
    IF = IF, B = B, mult_kind = 1L, seed = 7L)
  emp_var <- stats::var(out[, 1])
  theo_var <- sum(IF^2)
  expect_equal(emp_var, theo_var, tolerance = 0.05 * theo_var)
})

test_that("CS OR: single treated cohort still produces ATTs on GPU", {
  skip_if_no_cuda()
  # Construct a panel with exactly one treated cohort.
  set.seed(2L)
  n_units <- 40L; n_periods <- 8L
  F_g <- rep(Inf, n_units)
  treated <- sample(seq_len(n_units), 20L)
  F_g[treated] <- 5L  # all treated at the same period -> one cohort
  panel <- expand.grid(unit = seq_len(n_units), period = seq_len(n_periods))
  panel$D <- as.integer(panel$period >= F_g[panel$unit])
  panel$Y <- rnorm(nrow(panel)) + 1.0 * panel$D
  p <- panel[order(panel$unit, panel$period), c("unit", "period", "D", "Y")]
  fit_r <- didgpu_cs(p, "Y", "unit", "period", "D", est_method = "OR",
                      bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  fit_c <- didgpu_cs(p, "Y", "unit", "period", "D", est_method = "OR",
                      bootstrap_reps = 0L, backend = "cuda", verbose = FALSE)
  expect_equal(fit_r$att_gt$att, fit_c$att_gt$att, tolerance = 1e-12)
})
