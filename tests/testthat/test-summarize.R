test_that("didgpu_summarize_panel returns the right shape on a known panel", {
  p <- didgpu_simulate_panel_bidir(n_units = 80L, n_periods = 15L,
                                    frac_treated = 0.6, frac_in = 0.5,
                                    min_treat_period = 5L, max_treat_period = 10L,
                                    seed = 11L)
  s <- didgpu_summarize_panel(p, "Y", "unit", "period", "D", verbose = FALSE)
  expect_equal(s$n_units, 80L)
  expect_equal(s$n_periods, 15L)
  expect_equal(s$n_rows, 1200L)
  expect_true(s$is_balanced)
  # 60% of 80 = 48 switchers, half in / half out
  expect_equal(s$n_switchers_in, 24L)
  expect_equal(s$n_switchers_out, 24L)
  expect_equal(s$n_never_change, 32L)
  # Bidir simulator's d_sq: in-switchers + never-treated have d_sq=0
  # (24+32=56); out-switchers have d_sq=1 (24).
  expect_equal(as.integer(s$d_sq_dist), c(56L, 24L))
  # Sanity bounds: each individually <= n_periods. (The two can come
  # from different groups so their sum can exceed n_periods.)
  expect_lte(s$max_effects, 15L)
  expect_lte(s$max_placebo, 15L)
})

test_that("didgpu_summarize_panel handles unbalanced panel", {
  p <- didgpu_simulate_panel(n_units = 30L, n_periods = 10L, seed = 5L,
                              min_treat_period = 3L, max_treat_period = 6L)
  # Drop 5 random rows.
  p <- p[-sample(seq_len(nrow(p)), 5L), ]
  s <- didgpu_summarize_panel(p, "Y", "unit", "period", "D", verbose = FALSE)
  expect_false(s$is_balanced)
  expect_lt(s$n_rows, s$n_units * s$n_periods)
})

test_that("didgpu_summarize_panel errors on missing columns", {
  p <- didgpu_simulate_panel(n_units = 20L, n_periods = 8L)
  expect_error(
    didgpu_summarize_panel(p, "nonexistent", "unit", "period", "D",
                            verbose = FALSE),
    "not in df"
  )
})

test_that("didgpu_estimate_runtime returns scaled estimate by n_workers", {
  p <- didgpu_simulate_panel(n_units = 30L, n_periods = 8L, seed = 5L,
                              min_treat_period = 3L, max_treat_period = 6L)
  e1 <- didgpu_estimate_runtime(p, "Y", "unit", "period", "D",
                                  effects = 1L, bootstrap_reps = 10L,
                                  n_workers = 1L, probes = 1L)
  e4 <- didgpu_estimate_runtime(p, "Y", "unit", "period", "D",
                                  effects = 1L, bootstrap_reps = 10L,
                                  n_workers = 4L, probes = 1L)
  expect_equal(e1$n_remaining, 11L)
  expect_equal(e4$n_remaining, 11L)
  # n_workers=4 should give roughly 4x less total_seconds (plus the
  # 1.5s cluster overhead). At minimum, less than the sequential time.
  expect_lt(e4$total_seconds, e1$total_seconds + 2)
})
