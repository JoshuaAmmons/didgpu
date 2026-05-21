# Phase-2 #82 tests: GPU cluster bootstrap.
#
# The kernel implements the IF-shortcut variant of the cluster
# bootstrap. It's asymptotically equivalent to the per-rep R refit
# (.cs_bootstrap_se) but uses different RNG and a different finite-
# sample estimator, so per-rep matrices differ. The columnwise SDs
# should still converge as B grows.

test_that("didgpu_cuda_cluster_bootstrap_r returns the right shape", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  set.seed(11L)
  n_units <- 30L; n_dims <- 5L
  IF <- matrix(rnorm(n_units * n_dims), nrow = n_units, ncol = n_dims)
  cluster_id <- as.integer((seq_len(n_units) - 1L) %/% 3L)  # 10 clusters of 3
  out <- didgpu:::didgpu_cuda_cluster_bootstrap_r(
    IF         = IF,
    cluster_id = cluster_id,
    n_clusters = 10L,
    B          = 32L,
    seed       = 17L)
  expect_equal(dim(out), c(32L, n_dims))
  expect_true(all(is.finite(out)))
})

test_that("CUDA cluster bootstrap SEs agree with a CPU-side IF reference within MC error", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  # Build a reference IF-shortcut cluster bootstrap on the CPU using R's
  # MT19937, then compare to the CUDA kernel's columnwise SDs. Different
  # RNG, same algorithm — SDs should agree within ~2 / sqrt(B).
  set.seed(31L)
  n_units <- 50L; n_dims <- 4L; n_clusters <- 25L
  IF <- matrix(rnorm(n_units * n_dims), nrow = n_units, ncol = n_dims)
  cluster_id <- as.integer((seq_len(n_units) - 1L) %% n_clusters)
  B <- 2000L

  # CUDA bootstrap.
  cuda_mat <- didgpu:::didgpu_cuda_cluster_bootstrap_r(
    IF = IF, cluster_id = cluster_id, n_clusters = n_clusters,
    B = B, seed = 17L)
  cuda_sds <- apply(cuda_mat, 2L, stats::sd)

  # CPU reference: same IF-shortcut algorithm with R's RNG.
  set.seed(17L, kind = "Mersenne-Twister")
  IF_cluster <- matrix(0.0, nrow = n_clusters, ncol = n_dims)
  for (u in seq_len(n_units))
    IF_cluster[cluster_id[u] + 1L, ] <- IF_cluster[cluster_id[u] + 1L, ] + IF[u, ]
  cpu_mat <- matrix(NA_real_, nrow = B, ncol = n_dims)
  for (b in seq_len(B)) {
    picks <- sample.int(n_clusters, n_clusters, replace = TRUE)
    wt <- tabulate(picks, nbins = n_clusters)
    cpu_mat[b, ] <- as.numeric(wt %*% IF_cluster)
  }
  cpu_sds <- apply(cpu_mat, 2L, stats::sd)

  # Asymptotic SD of the SD estimator is ~ SD / sqrt(2 * (B - 1)),
  # which is ~0.016 * SD here. Allow 0.15 (about 9x SE) — the GPU and
  # CPU produce different per-replicate draws but the population SD
  # is the same.
  rel_diff <- abs(cuda_sds - cpu_sds) / cpu_sds
  expect_true(max(rel_diff) < 0.15,
              info = paste("max rel diff =", max(rel_diff)))
})

test_that("didgpu_cs(backend='cuda', bootstrap_reps>0) populates SE columns", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  p <- didgpu_simulate_panel(n_units = 30L, n_periods = 6L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  p$D <- as.integer(p$D >= 0.5)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR", aggregation = "event",
                    bootstrap_reps = 200L,
                    bootstrap_kind = "cluster",
                    backend = "cuda", verbose = FALSE)
  expect_true("se" %in% names(fit$att_gt))
  expect_true(any(!is.na(fit$att_gt$se)),
              info = "at least some SEs should be populated")
})
