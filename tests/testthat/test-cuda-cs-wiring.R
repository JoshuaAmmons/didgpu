# Phase-1 wiring test for task #79.
#
# The CUDA cs_inner_batched kernel is currently a scaffold that
# returns -1 ("not implemented"). The R-side dispatch in
# `.cs_compute_att_gt` is supposed to detect that and fall back to the
# per-cell R loop, producing bit-identical results to a backend = "r"
# run. These tests pin that contract so that:
#   - Phase 2 work on the kernel can't break the fallback by accident.
#   - Future signature changes to .cs_inner_batched_cuda fail loudly.

test_that("didgpu_cs(backend = 'cuda') falls back cleanly when kernel returns -1", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")

  p <- didgpu_simulate_panel(n_units = 20L, n_periods = 6L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  p$D <- as.integer(p$D >= 0.5)

  fit_r <- didgpu_cs(p, "Y", "unit", "period", "D",
                      est_method = "OR", aggregation = "event",
                      bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  fit_c <- didgpu_cs(p, "Y", "unit", "period", "D",
                      est_method = "OR", aggregation = "event",
                      bootstrap_reps = 0L, backend = "cuda", verbose = FALSE)

  expect_equal(nrow(fit_r$att_gt), nrow(fit_c$att_gt))
  expect_equal(fit_r$att_gt$att, fit_c$att_gt$att, tolerance = 1e-12)
  expect_equal(fit_r$att_gt$n_treated, fit_c$att_gt$n_treated)
  expect_equal(fit_r$att_gt$n_control, fit_c$att_gt$n_control)
})

test_that(".cs_inner_batched_cuda returns NULL while kernel is scaffolded", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  cells <- list(
    list(delta = c(1.0, 2.0, -0.5), D_mask = c(TRUE, FALSE, FALSE),
         X = NULL, n_total = 3L, units = c(1L, 2L, 3L)),
    list(delta = c(0.7, 1.3),       D_mask = c(TRUE, FALSE),
         X = NULL, n_total = 2L, units = c(1L, 2L))
  )
  result <- didgpu:::.cs_inner_batched_cuda(
    cells = cells, method = "OR", all_units = c(1L, 2L, 3L))
  # Kernel currently returns -1 -> Rcpp wrapper returns NULL ->
  # R-side helper returns NULL. Phase 2 will flip this to a list.
  expect_null(result)
})

test_that("Rcpp wrapper didgpu_cuda_cs_inner_batched_r exists and accepts inputs", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  # Smoke: a 2-cell, 1-covariate (just intercept) layout. Kernel
  # currently returns -1 so the Rcpp wrapper returns NULL.
  # The wrapper is an internal symbol — access via didgpu:::.
  result <- didgpu:::didgpu_cuda_cs_inner_batched_r(
    X_concat       = rep(1.0, 5),     # n_total = 5, p = 1 -> length 5
    X_offsets      = c(0L, 3L, 5L),   # cells: [0,3), [3,5)
    Y_concat       = c(1, 2, -0.5, 0.7, 1.3),
    W_concat       = c(1, 0, 0, 1, 0),
    p              = 1L,
    n_units        = 3L,
    est_method     = 0L,
    want_influence = TRUE)
  expect_null(result)
})
