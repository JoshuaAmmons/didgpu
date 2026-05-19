# Leave-one-out robustness analysis across all five estimator families.

test_that("didgpu_loo on didgpu_cs drops one cohort at a time", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR", aggregation = "overall",
                    bootstrap_reps = 0L,
                    backend = "r", verbose = FALSE)
  loo <- didgpu_loo(fit, by = "cohort", df = p, verbose = FALSE)
  expect_s3_class(loo, "didgpu_loo_result")
  expect_true(all(c("leave_out", "estimate", "delta", "delta_pct") %in%
                   names(loo)))
  # One row per cohort.
  n_cohorts <- length(unique(stats::na.omit(p$D * p$period)[p$D == 1L]))
  expect_gt(nrow(loo), 1L)
  expect_equal(attr(loo, "by"), "cohort (F_g)")
})

test_that("didgpu_loo on didgpu drops cohorts and reports finite delta", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 2L, bootstrap_reps = 0L,
                 backend = "r", verbose = FALSE)
  loo <- didgpu_loo(fit, by = "cohort", df = p, verbose = FALSE)
  expect_true(any(is.finite(loo$delta)))
})

test_that("didgpu_loo on didgpu_fect works", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  fit <- didgpu_fect(p, "Y", "unit", "period", "D",
                      method = "fe", effects = 1L,
                      bootstrap_reps = 0L,
                      backend = "r", verbose = FALSE)
  loo <- didgpu_loo(fit, by = "cohort", df = p, verbose = FALSE)
  expect_s3_class(loo, "didgpu_loo_result")
  expect_gt(nrow(loo), 1L)
})

test_that("didgpu_loo accepts an arbitrary column as the drop key", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  p$region <- ifelse(p$unit %% 2L == 0L, "north", "south")
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR", aggregation = "overall",
                    bootstrap_reps = 0L,
                    backend = "r", verbose = FALSE)
  loo <- didgpu_loo(fit, by = "region", df = p, verbose = FALSE)
  # Two levels of region; two rows in the LOO.
  expect_equal(nrow(loo), 2L)
  expect_setequal(loo$leave_out, c("north", "south"))
})

test_that("didgpu_loo errors clearly without df", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 8L, seed = 17L)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR",
                    bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_error(didgpu_loo(fit, by = "cohort"), "original panel")
})

test_that("didgpu_loo errors on unrecognised fit class", {
  expect_error(didgpu_loo(list(x = 1), by = "cohort", df = data.frame()),
                "didgpu_result")
})

test_that("didgpu_loo sorts results by abs(delta) descending", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR", aggregation = "overall",
                    bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  loo <- didgpu_loo(fit, by = "cohort", df = p, verbose = FALSE)
  abs_d <- abs(loo$delta)
  ok <- !is.na(abs_d)
  if (sum(ok) >= 2L) {
    diffs <- diff(abs_d[ok])
    expect_true(all(diffs <= 1e-10))   # monotonically non-increasing
  }
})

test_that("print.didgpu_loo_result runs without error", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 8L, seed = 17L)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR", aggregation = "overall",
                    bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  loo <- didgpu_loo(fit, by = "cohort", df = p, verbose = FALSE)
  expect_output(print(loo), "Leave-one-out")
  expect_output(print(loo), "full-sample estimate")
})
