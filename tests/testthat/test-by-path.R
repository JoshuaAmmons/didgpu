# didgpu_compute_paths + didgpu_by_path: treatment-trajectory subgroup
# analysis. Augments the panel with a per-group `path` column built
# from the treatment values at (F_g - 1, F_g, ..., F_g + effects - 1),
# then runs didgpu_by on the path column.

test_that("didgpu_compute_paths adds a non-empty path column", {
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0),
                              seed = 17L)
  aug <- didgpu_compute_paths(p, "Y", "unit", "period", "D",
                                effects = 2L)
  expect_true("path" %in% names(aug))
  expect_equal(nrow(aug), nrow(p))
  # Every row of a given unit has the same path.
  paths_per_unit <- aggregate(path ~ unit, data = aug,
                               function(x) length(unique(x)))
  expect_true(all(paths_per_unit$path == 1L))
  # At least two distinct paths (treated vs. never-switched).
  expect_gt(length(unique(aug$path)), 1L)
})

test_that("didgpu_compute_paths respects top_n", {
  # Build a panel with multivalued treatment so we get >2 distinct paths
  # (binary panels only ever produce "0,1,..." and "no_switch").
  set.seed(23L)
  n_units <- 80L; n_periods <- 12L
  fg <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), 60L))
  fg[treated] <- sample(3L:8L, 60L, replace = TRUE)
  doses <- sample(c(1L, 2L, 3L), n_units, replace = TRUE)
  panel <- data.table::data.table(
    unit = rep(seq_len(n_units), each = n_periods),
    period = rep(seq_len(n_periods), n_units))
  panel[, F_g := fg[unit]]
  panel[, dose := doses[unit]]
  panel[, D := as.integer(ifelse(period >= F_g, dose, 0L))]
  panel[, Y := rnorm(.N)]
  p <- as.data.frame(panel[order(unit, period), .(unit, period, D, Y)])
  aug <- didgpu_compute_paths(p, "Y", "unit", "period", "D",
                                effects = 2L, top_n = 2L)
  # At most 2 non-NA paths.
  expect_lte(length(unique(stats::na.omit(aug$path))), 2L)
  # The multivalued panel has 3 doses × possibly several F_g periods, so
  # there should be more than 2 paths total — some must have been NA'd.
  expect_true(any(is.na(aug$path)))
})

test_that("never-switchers get a 'no_switch' path", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L,
                              frac_treated = 0.5,  # half never switch
                              tau_profile = c(0.5),
                              seed = 17L)
  aug <- didgpu_compute_paths(p, "Y", "unit", "period", "D",
                                effects = 1L)
  expect_true("no_switch" %in% aug$path)
})

test_that("didgpu_by_path runs per-trajectory estimation", {
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0),
                              seed = 17L)
  fit <- suppressWarnings(didgpu_by_path(
    p, outcome = "Y", group = "unit", time = "period", treatment = "D",
    effects = 2L, top_n = 3L,
    bootstrap_reps = 0L, backend = "r", verbose = FALSE
  ))
  expect_s3_class(fit, "didgpu_by_result")
  expect_gt(length(fit), 0L)
  # Each subgroup is either a didgpu_result or a didgpu_by_failed
  # sentinel (the latter for degenerate paths like "no_switch" which
  # have no switchers and can't produce effect estimates).
  for (nm in names(fit)) {
    expect_true(inherits(fit[[nm]], "didgpu_result") ||
                  inherits(fit[[nm]], "didgpu_by_failed"))
  }
  # At least one subgroup should succeed (the actually-switching path).
  succeeded <- vapply(fit, function(r) inherits(r, "didgpu_result"),
                       logical(1))
  expect_true(any(succeeded))
})

test_that("didgpu_by_path errors if all paths get filtered out", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L,
                              tau_profile = c(0.5),
                              seed = 17L)
  expect_error(
    didgpu_by_path(p, outcome = "Y", group = "unit", time = "period",
                    treatment = "D", effects = 1L, top_n = 0L,
                    bootstrap_reps = 0L, backend = "r", verbose = FALSE),
    "no rows remain"
  )
})
