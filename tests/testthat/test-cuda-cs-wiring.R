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
    # Cell 1: 2 treated, 2 control. Y_t = (1.0, 3.0), Y_c = (2.0, 0.0).
    list(delta = c(1.0, 3.0, 2.0, 0.0), D_mask = c(TRUE, TRUE, FALSE, FALSE),
         X = NULL, n_total = 4L, units = c(1L, 2L, 3L, 4L)),
    # Cell 2: 1 treated, 1 control.
    list(delta = c(0.7, 1.3),       D_mask = c(TRUE, FALSE),
         X = NULL, n_total = 2L, units = c(1L, 2L))
  )
  result <- didgpu:::.cs_inner_batched_cuda(
    cells = cells, method = "OR", all_units = c(1L, 2L, 3L, 4L))
  # Phase 2 #84: OR kernel now returns real ATTs (was NULL before).
  # Cell 1: mean_Yt = 2.0, mean_Yc = 1.0, ATT = 1.0.
  #   IF[unit 1] (treated row 0, Y=1.0) = Y - mean_Yt = 1.0 - 2.0 = -1.0
  #   IF[unit 2] (treated row 1, Y=3.0) = 3.0 - 2.0 =  1.0
  #   IF[units 3, 4] (control)          = 0.0
  # Cell 2: ATT = 0.7 - 1.3 = -0.6; one treated with Y = mean_Yt = 0.7
  #   so IF = 0 for everyone.
  expect_false(is.null(result))
  expect_equal(result$att, c(1.0, -0.6), tolerance = 1e-12)
  expect_equal(dim(result$influence), c(4L, 2L))
  expect_equal(result$influence[1, 1], -1.0, tolerance = 1e-12)
  expect_equal(result$influence[2, 1],  1.0, tolerance = 1e-12)
  expect_equal(result$influence[3, 1],  0.0)
  expect_equal(result$influence[4, 1],  0.0)
})

test_that("Rcpp wrapper didgpu_cuda_cs_inner_batched_r computes ATT and IF for OR", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  # 2 cells, intercept-only design (p = 1). Cell 0: rows 0..2 with
  # treated row 0; Cell 1: rows 3..4 with treated row 3.
  result <- didgpu:::didgpu_cuda_cs_inner_batched_r(
    X_concat        = rep(1.0, 5),     # p = 1 intercept column
    X_offsets       = c(0L, 3L, 5L),   # cells: [0,3), [3,5)
    Y_concat        = c(2.0, 1.0, -1.0, 5.0, 2.0),
    W_concat        = c(1, 0, 0, 1, 0),
    unit_id_per_row = c(0L, 1L, 2L, 0L, 2L),
    p               = 1L,
    n_units         = 3L,
    est_method      = 0L,
    want_influence  = TRUE)
  expect_false(is.null(result))
  expect_named(result, c("att", "status", "influence"), ignore.order = TRUE)
  # Cell 0: ATT = mean(Y_t) - mean(Y_c) = 2.0 - mean(1, -1) = 2.0 - 0 = 2.0.
  # Cell 1: ATT = 5.0 - 2.0 = 3.0.
  expect_equal(result$att, c(2.0, 3.0), tolerance = 1e-10)
  expect_equal(dim(result$influence), c(3L, 2L))
})

test_that("Rcpp wrapper computes IPW / DR (no-cov closed form)", {
  skip_if_not(isTRUE(tryCatch(didgpu_has_cuda_support(),
                              error = function(e) FALSE)),
              "CUDA support not compiled into this build")
  # p = 1 (intercept only) -> IPW and DR both reduce to the simple
  # mean-difference with the SPECIAL no-cov influence function:
  #   att = mean(Y_t) - mean(Y_c)
  #   IF[treated]  = Y - mean_t - att/2
  #   IF[control]  = -(Y - mean_c) - att/2
  # Cell: treated rows 0 & 3, control rows 1, 2, 4.
  # Y_t = (1.0, 0.7) -> mean_t = 0.85; Y_c = (2, -0.5, 1.3) -> mean_c = 0.9333..
  # But these are TWO cells: [0,3) and [3,5). Compute per cell.
  for (m in c(1L, 2L)) {
    res <- didgpu:::didgpu_cuda_cs_inner_batched_r(
      X_concat        = rep(1.0, 5), X_offsets = c(0L, 3L, 5L),
      Y_concat        = c(1, 2, -0.5, 0.7, 1.3),
      W_concat        = c(1, 0, 0, 1, 0),
      unit_id_per_row = c(0L, 1L, 2L, 0L, 2L),
      p = 1L, n_units = 3L, est_method = m, want_influence = TRUE)
    expect_false(is.null(res))
    # Cell 1: treated {1.0}, control {2, -0.5} -> att = 1.0 - 0.75 = 0.25.
    # Cell 2: treated {0.7}, control {1.3}     -> att = 0.7 - 1.3  = -0.6.
    expect_equal(res$att, c(0.25, -0.6), tolerance = 1e-12,
                 info = sprintf("est_method=%d", m))
  }
})
