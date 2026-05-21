test_that("didgpu_tidy returns expected columns and rows", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 5L,
                              min_treat_period = 3L, max_treat_period = 7L,
                              tau_profile = c(0.5, 1.0))
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 2L, placebo = 1L,
                 bootstrap_reps = 5L, seed = 1L,
                 backend = "r", verbose = FALSE)

  td <- didgpu_tidy(fit)
  expect_s3_class(td, "data.frame")
  expect_equal(names(td),
               c("term", "estimate", "std.error", "statistic",
                 "p.value", "conf.low", "conf.high", "kind"))
  # 2 effects + 1 ATE (since effects > 1, ATE is included) + 1 placebo = 4 rows.
  expect_equal(nrow(td), 4L)
  expect_setequal(td$kind, c("effect", "ate", "placebo"))
  # Estimates from tidy must equal the matrix Estimates.
  expect_equal(td$estimate[1:2],
               as.numeric(fit$results$Effects[, "Estimate"]))
  expect_equal(td$estimate[td$kind == "placebo"],
               as.numeric(fit$results$Placebos[, "Estimate"]))
})

test_that("didgpu_tidy can suppress CI columns", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 5L,
                              min_treat_period = 3L, max_treat_period = 7L)
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 1L, bootstrap_reps = 2L, seed = 1L,
                 backend = "r", verbose = FALSE)
  td <- didgpu_tidy(fit, conf.int = FALSE)
  expect_false("conf.low" %in% names(td))
  expect_false("conf.high" %in% names(td))
})

test_that("didgpu_glance returns a one-row summary", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 5L,
                              min_treat_period = 3L, max_treat_period = 7L)
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 2L, placebo = 1L,
                 bootstrap_reps = 5L, seed = 1L,
                 backend = "r", verbose = FALSE)
  gl <- didgpu_glance(fit)
  expect_equal(nrow(gl), 1L)
  expect_true(all(c("n_effects", "n_placebos", "n_switchers_e1",
                    "p_jointeffects", "n_boot", "backend") %in% names(gl)))
  expect_equal(gl$n_effects, 2L)
  expect_equal(gl$n_placebos, 1L)
  expect_equal(gl$n_boot, 5L)
  expect_equal(gl$backend, "r")
})

test_that("tidy.didgpu_result S3 method dispatches correctly", {
  p <- didgpu_simulate_panel(n_units = 30L, n_periods = 8L, seed = 5L,
                              min_treat_period = 3L, max_treat_period = 5L)
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 1L, bootstrap_reps = 3L, seed = 1L,
                 backend = "r", verbose = FALSE)
  via_method <- tidy.didgpu_result(fit)
  via_direct <- didgpu_tidy(fit)
  expect_equal(via_method, via_direct)
})
