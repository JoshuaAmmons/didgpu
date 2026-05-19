# ============================================================================
# Phase 4 #90: CUDA-vs-R equivalence grid.
#
# For every estimator path that actually computes on the GPU (not just
# falls back), assert backend = "cuda" agrees with backend = "r"
# within a documented tolerance, across several panel sizes.
#
# Tolerance rationale:
#   - CS OR, no covariates (p = 1): the kernel computes mean(Y_t) -
#     mean(Y_c) with the same arithmetic as R -> bit-exact (1e-12).
#   - CS OR, with covariates (p > 1): the kernel solves the normal
#     equations by Cholesky; R uses QR (qr.solve). Both are backward-
#     stable but round differently -> 1e-6.
#   - Bootstrap SEs use a different RNG stream on GPU vs CPU, so they
#     are NOT bit-comparable; their equivalence is covered
#     statistically in the dedicated bootstrap test files. This grid
#     covers the DETERMINISTIC point-estimate paths only.
# ============================================================================

cuda_available <- function() {
  isTRUE(tryCatch(didgpu_has_cuda_support(), error = function(e) FALSE))
}

make_cs_panel <- function(n_units, n_periods, seed = 17L, with_cov = FALSE) {
  p <- didgpu_simulate_panel(n_units = n_units, n_periods = n_periods,
                              tau_profile = c(0.5, 1.0), seed = seed)
  p$D <- as.integer(p$D >= 0.5)
  if (with_cov) {
    # Time-invariant unit-level covariate (CS uses the value at T_min).
    set.seed(seed + 1L)
    uvals <- stats::rnorm(length(unique(p$unit)))
    names(uvals) <- as.character(sort(unique(p$unit)))
    p$x1 <- uvals[as.character(p$unit)]
  }
  p
}

test_that("CS OR point estimate (no covariates) is bit-exact CUDA vs R", {
  skip_if_not(cuda_available(), "CUDA not compiled in")
  for (nu in c(40L, 100L, 200L)) {
    for (np in c(8L, 12L)) {
      p <- make_cs_panel(nu, np)
      fit_r <- didgpu_cs(p, "Y", "unit", "period", "D",
                          est_method = "OR", aggregation = "event",
                          bootstrap_reps = 0L, backend = "r", verbose = FALSE)
      fit_c <- didgpu_cs(p, "Y", "unit", "period", "D",
                          est_method = "OR", aggregation = "event",
                          bootstrap_reps = 0L, backend = "cuda", verbose = FALSE)
      info <- sprintf("nu=%d np=%d", nu, np)
      expect_equal(nrow(fit_r$att_gt), nrow(fit_c$att_gt), info = info)
      expect_equal(fit_r$att_gt$att, fit_c$att_gt$att,
                   tolerance = 1e-12, info = info)
      # Aggregated event-study estimates must also match.
      expect_equal(fit_r$aggregation$estimate, fit_c$aggregation$estimate,
                   tolerance = 1e-12, info = info)
    }
  }
})

test_that("CS OR point estimate (with covariates) matches CUDA vs R within 1e-6", {
  skip_if_not(cuda_available(), "CUDA not compiled in")
  for (nu in c(60L, 150L)) {
    for (np in c(8L, 12L)) {
      p <- make_cs_panel(nu, np, with_cov = TRUE)
      fit_r <- didgpu_cs(p, "Y", "unit", "period", "D",
                          covariates = "x1",
                          est_method = "OR", aggregation = "event",
                          bootstrap_reps = 0L, backend = "r", verbose = FALSE)
      fit_c <- didgpu_cs(p, "Y", "unit", "period", "D",
                          covariates = "x1",
                          est_method = "OR", aggregation = "event",
                          bootstrap_reps = 0L, backend = "cuda", verbose = FALSE)
      info <- sprintf("nu=%d np=%d", nu, np)
      expect_equal(nrow(fit_r$att_gt), nrow(fit_c$att_gt), info = info)
      # Cholesky (CUDA) vs QR (R): backward-stable, ~1e-6 agreement.
      expect_equal(fit_r$att_gt$att, fit_c$att_gt$att,
                   tolerance = 1e-6, info = info)
    }
  }
})

test_that("CS OR influence functions match CUDA vs R (no covariates)", {
  skip_if_not(cuda_available(), "CUDA not compiled in")
  # The per-cell IF drives the bootstrap SEs, so its equivalence is
  # what makes the CUDA cluster/multiplier bootstrap trustworthy.
  p <- make_cs_panel(100L, 10L)
  fit_r <- didgpu_cs(p, "Y", "unit", "period", "D",
                      est_method = "OR", aggregation = "event",
                      bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  fit_c <- didgpu_cs(p, "Y", "unit", "period", "D",
                      est_method = "OR", aggregation = "event",
                      bootstrap_reps = 0L, backend = "cuda", verbose = FALSE)
  IF_r <- attr(fit_r$att_gt, "IF_per_cell")
  IF_c <- attr(fit_c$att_gt, "IF_per_cell")
  expect_equal(length(IF_r), length(IF_c))
  for (i in seq_along(IF_r)) {
    # Cells appear in the same order; compare unit-aligned IF vectors.
    expect_equal(IF_r[[i]]$units, IF_c[[i]]$units,
                 info = sprintf("cell %d units", i))
    expect_equal(IF_r[[i]]$IF, IF_c[[i]]$IF, tolerance = 1e-12,
                 info = sprintf("cell %d IF", i))
  }
})

test_that("SAXPY smoke kernel is exact CUDA vs R reference", {
  skip_if_not(cuda_available(), "CUDA not compiled in")
  for (n in c(1L, 5L, 100L, 1000L)) {
    a <- 2.5
    x <- as.numeric(seq_len(n))
    y <- as.numeric(rev(seq_len(n)))
    got <- didgpu_run_saxpy(a, x, y)
    want <- a * x + y
    expect_equal(got, want, tolerance = 1e-5,
                 info = sprintf("n=%d", n))  # float32 kernel -> 1e-5
  }
})

test_that("fect size gate keeps CUDA == R for small panels (transparent fallback)", {
  skip_if_not(cuda_available(), "CUDA not compiled in")
  # Below the .fect_cuda_svd_worthwhile threshold, backend = "cuda"
  # uses R svd(), so results must be identical to backend = "r".
  for (method in c("fe", "ife", "mc")) {
    p <- make_cs_panel(80L, 12L)
    fit_r <- didgpu_fect(p, "Y", "unit", "period", "D", method = method,
                          backend = "r", bootstrap_reps = 0L, verbose = FALSE)
    fit_c <- didgpu_fect(p, "Y", "unit", "period", "D", method = method,
                          backend = "cuda", bootstrap_reps = 0L, verbose = FALSE)
    info <- sprintf("method=%s", method)
    expect_equal(fit_r$results$ATE[1, "Estimate"],
                 fit_c$results$ATE[1, "Estimate"],
                 tolerance = 1e-10, info = info)
  }
})
