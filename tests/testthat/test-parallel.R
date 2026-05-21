# Parallel bootstrap should produce bit-identical results to sequential.

# CRAN's R CMD check sets _R_CHECK_LIMIT_CORES_=TRUE and forbids spawning more
# than 2 worker processes; use 2 there, 4 in normal/dev runs. Either way this
# exercises the n_workers > 1 path against the sequential baseline.
nw_multi <- if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) 2L else 4L

test_that("n_workers > 1 produces identical Effects/Placebos/SEs to n_workers = 1", {
  p <- didgpu_simulate_panel(n_units = 50L, n_periods = 12L,
                              frac_treated = 0.6,
                              min_treat_period = 4L, max_treat_period = 8L,
                              tau_profile = c(0.5, 1.0),
                              sigma = 0.4, seed = 17L)

  fit_seq <- didgpu(p, "Y", "unit", "period", "D",
                     effects = 2L, placebo = 1L,
                     bootstrap_reps = 12L, seed = 1L,
                     backend = "r", verbose = FALSE, n_workers = 1L)
  fit_par <- didgpu(p, "Y", "unit", "period", "D",
                     effects = 2L, placebo = 1L,
                     bootstrap_reps = 12L, seed = 1L,
                     backend = "r", verbose = FALSE, n_workers = nw_multi)

  expect_equal(fit_seq$results$Effects[, "Estimate"],
               fit_par$results$Effects[, "Estimate"])
  expect_equal(fit_seq$results$Effects[, "SE"],
               fit_par$results$Effects[, "SE"])
  expect_equal(fit_seq$results$Placebos[, "Estimate"],
               fit_par$results$Placebos[, "Estimate"])
  expect_equal(fit_seq$results$ATE[1, "Estimate"],
               fit_par$results$ATE[1, "Estimate"])
})

test_that("parallel run with checkpointing supports resume", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 5L,
                              min_treat_period = 3L, max_treat_period = 6L)
  cdir <- tempfile("didgpu_par_resume_")
  on.exit(unlink(cdir, recursive = TRUE), add = TRUE)

  fit1 <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 2L, bootstrap_reps = 8L, seed = 1L,
                  checkpoint_dir = cdir,
                  backend = "r", verbose = FALSE, n_workers = 2L)
  # Resume with the SAME config — should add no new cells, return identical.
  fit2 <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 2L, bootstrap_reps = 8L, seed = 1L,
                  checkpoint_dir = cdir,
                  backend = "r", verbose = FALSE, n_workers = nw_multi)
  expect_equal(fit1$results$Effects[, "Estimate"],
               fit2$results$Effects[, "Estimate"])
  expect_equal(fit1$results$Effects[, "SE"],
               fit2$results$Effects[, "SE"])
})
