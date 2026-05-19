test_that("simulator produces a balanced panel with the requested treatment pattern", {
  p <- didgpu_simulate_panel(n_units = 50L, n_periods = 10L,
                              frac_treated = 0.5,
                              min_treat_period = 4L, max_treat_period = 7L,
                              tau_profile = c(1.0, 1.0, 1.0),
                              sigma = 0.0,            # noiseless
                              unit_fe_sd = 0, time_fe_sd = 0,
                              seed = 1L)
  expect_equal(nrow(p), 50L * 10L)
  expect_true(all(c("unit", "period", "D", "Y") %in% names(p)))
  expect_true(all(p$D %in% c(0L, 1L)))
  truth <- attr(p, "truth")
  expect_true(!is.null(truth))
  expect_equal(length(truth$F_g), 50L)
  expect_true(any(is.infinite(truth$F_g)))   # at least one never-treated
})

test_that("simulator respects auto-scaling defaults on tiny panels", {
  # n_periods = 8 used to fail because defaults said max_treat_period = 15.
  expect_silent(p <- didgpu_simulate_panel(n_units = 30L, n_periods = 8L,
                                            seed = 2L))
  expect_equal(nrow(p), 30L * 8L)
})

test_that("simulator recovers tau when noise is zero", {
  p <- didgpu_simulate_panel(n_units = 200L, n_periods = 12L,
                              frac_treated = 0.5,
                              min_treat_period = 5L, max_treat_period = 6L,
                              tau_profile = c(0.5, 1.0, 1.5),
                              sigma = 0, unit_fe_sd = 0, time_fe_sd = 0,
                              seed = 7L)
  # In the noiseless case, Y at (unit, period) == tau(period - F_g) for
  # treated units with t >= F_g. Spot check a few rows.
  truth <- attr(p, "truth")
  for (u in head(which(is.finite(truth$F_g)), 3L)) {
    f <- truth$F_g[as.character(u)]
    rows <- p[p$unit == u, , drop = FALSE]
    for (t in seq_len(nrow(rows))) {
      per <- rows$period[t]
      if (per < f) {
        expect_equal(rows$Y[t], 0)
      } else {
        k <- per - f
        expected_tau <- truth$tau_profile[min(k + 1L, length(truth$tau_profile))]
        expect_equal(rows$Y[t], expected_tau, tolerance = 1e-10)
      }
    }
  }
})
