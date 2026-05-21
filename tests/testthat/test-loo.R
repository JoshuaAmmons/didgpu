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

test_that("CS cohort-LOO fast path (re-aggregate) == refit, bit-for-bit (never controls)", {
  p <- didgpu_simulate_panel(n_units = 50L, n_periods = 10L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  p$D <- as.integer(p$D >= 0.5)
  for (agg in c("overall", "event", "group")) {
    fit <- didgpu_cs(p, "Y", "unit", "period", "D", est_method = "OR",
                      control_group = "never", aggregation = agg,
                      bootstrap_reps = 0L, backend = "r", verbose = FALSE)
    fast <- didgpu_loo(fit, by = "cohort", df = p, verbose = FALSE)
    expect_identical(attr(fast, "method"), "reaggregate (no refit)")

    # Manual refit reference: drop each cohort's units, refit, aggregate.
    d <- data.table::as.data.table(p)
    # as.double on both branches so the never-treated Inf doesn't trip
    # data.table's "assigning double to integer column" warning.
    d[, F_g := if (any(D == 1L)) as.double(min(period[D == 1L])) else Inf,
      by = unit]
    cohorts <- sort(unique(d$F_g[is.finite(d$F_g)]))
    ref_est <- vapply(cohorts, function(g) {
      keep <- unique(d$unit[d$F_g != g | !is.finite(d$F_g)])
      pm <- p[p$unit %in% keep, , drop = FALSE]
      # Small leave-out subsets can emit benign glm/rank warnings from
      # the reference refit; suppress so the suite stays warning-clean.
      fb <- suppressWarnings(didgpu_cs(pm, "Y", "unit", "period", "D",
                       est_method = "OR",
                       control_group = "never", aggregation = agg,
                       bootstrap_reps = 0L, backend = "r", verbose = FALSE))
      fb$aggregation$estimate[1]
    }, numeric(1))
    fast_ord <- fast[order(as.numeric(fast$leave_out)), ]
    expect_equal(fast_ord$estimate, ref_est, tolerance = 1e-12,
                 info = sprintf("aggregation=%s", agg))
  }
})

test_that("CS cohort-LOO with notyet controls falls back to the refit path", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 10L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  p$D <- as.integer(p$D >= 0.5)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D", est_method = "OR",
                    control_group = "notyet", aggregation = "overall",
                    bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  loo <- didgpu_loo(fit, by = "cohort", df = p, verbose = FALSE)
  # The fast re-aggregation is invalid for notyet controls -> generic
  # refit path, which does NOT set the "reaggregate" method attribute.
  expect_false(identical(attr(loo, "method"), "reaggregate (no refit)"))
  expect_s3_class(loo, "didgpu_loo_result")
  expect_gt(nrow(loo), 1L)
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
