# Reference parity on UNBALANCED panels.
#
# Regression test for the baseline-treatment bug: d_sq_XX was keyed off the
# GLOBAL first period (min(time_XX)), so any group entering the panel after
# that period had treatment_XX = NA there, got d_sq_XX = NA, fell through
# F_g_XX to T_max + 1, and was silently reclassified as a NEVER-SWITCHER.
# Every late-entrant switcher was dropped from the estimand.
#
# Balanced-panel parity was exact throughout (on a balanced panel a group's
# own first period IS the global first period), which is why the randomized
# differential suite never caught this: didgpu_simulate_panel() could only
# produce balanced panels until late_entry_frac was added alongside this fix.

# ---- simulator support -----------------------------------------------------

test_that("simulator can produce a genuinely unbalanced panel", {
  p <- didgpu_simulate_panel(
    n_units = 60L, n_periods = 12L, frac_treated = 0.6,
    min_treat_period = 4L, max_treat_period = 9L,
    seed = 17L, late_entry_frac = 0.4
  )
  expect_lt(nrow(p), 60L * 12L)
  expect_gt(length(unique(as.integer(table(p$unit)))), 1L)
  expect_false(is.null(attr(p, "truth")$entry_period))
})

test_that("late_entry_frac = 0 leaves seeded fixtures bit-identical", {
  # The unbalancing block must not touch the RNG stream when disabled,
  # or every existing seeded test fixture silently changes.
  a <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 3L)
  b <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 3L,
                             late_entry_frac = 0)
  expect_identical(a, b)
})

# ---- the core regression, no reference package required --------------------

test_that("a late-entrant switcher is estimated, not silently discarded", {
  # 4 late entrants (absent at t = 1) that switch at t = 4, and 8
  # never-treated controls present throughout. Under the bug the 4
  # switchers were reclassified as never-switchers, leaving NO switchers
  # at all.
  set.seed(11L)
  late <- do.call(rbind, lapply(1:4, function(u) {
    data.frame(unit = u, period = 2:6,
               D = as.integer(2:6 >= 4L), Y = rnorm(5L))
  }))
  ctrl <- do.call(rbind, lapply(5:12, function(u) {
    data.frame(unit = u, period = 1:6, D = 0L, Y = rnorm(6L))
  }))
  p <- rbind(late, ctrl)

  for (bk in c("reference", "r")) {
    fit <- didgpu(df = p, outcome = "Y", group = "unit", time = "period",
                  treatment = "D", effects = 1L, placebo = 0L,
                  bootstrap_reps = 0L, backend = bk, verbose = FALSE)
    expect_true(is.finite(fit$results$Effects[1L, "Estimate"]), info = bk)
  }

  # And they must be COUNTED as switchers. Under the bug d_sq_XX was NA
  # for all four, so they were reclassified as never-switchers.
  expect_equal(unname(fit$results$Effects[1L, "Switchers"]), 4)
})

test_that("baseline treatment uses each group's own first observed period", {
  # A late entrant whose treatment is already 1 at its first observed
  # period is an always-treated unit from its own baseline, NOT a
  # switcher-in. Getting this wrong flips its classification.
  set.seed(12L)
  late <- do.call(rbind, lapply(1:4, function(u) {
    data.frame(unit = u, period = 3:8, D = 1L, Y = rnorm(6L))
  }))
  sw <- do.call(rbind, lapply(5:10, function(u) {
    data.frame(unit = u, period = 1:8,
               D = as.integer(1:8 >= 5L), Y = rnorm(8L))
  }))
  ctrl <- do.call(rbind, lapply(11:18, function(u) {
    data.frame(unit = u, period = 1:8, D = 0L, Y = rnorm(8L))
  }))
  p <- rbind(late, sw, ctrl)

  fit <- didgpu(df = p, outcome = "Y", group = "unit", time = "period",
                treatment = "D", effects = 2L, placebo = 0L,
                bootstrap_reps = 0L, backend = "reference", verbose = FALSE)
  expect_true(all(is.finite(fit$results$Effects[, "Estimate"])))
})

# ---- parity against DIDmultiplegtDYN ---------------------------------------

test_that("all CPU backends match DIDmultiplegtDYN on an unbalanced panel", {
  skip_if_no_reference()

  # Reproduces the reported failure. NOTE the backend loop: backend
  # "reference" agreed with the oracle even on the BROKEN build, because
  # the NA d_sq only corrupts results via the fast backends' (time, d_sq)
  # cohort-key encoding (backend.R). A parity test that pinned only
  # backend = "reference" passed throughout and caught nothing.
  sim <- as.data.frame(didgpu_simulate_panel(n_units = 80L, n_periods = 12L,
                                             seed = 11L))
  set.seed(99)
  drop_idx <- unlist(lapply(split(seq_len(nrow(sim)), sim$unit), function(ix) {
    if (runif(1) < 0.4) utils::head(ix, sample(1:3, 1)) else integer(0)
  }))
  p <- sim[-drop_idx, ]
  expect_gt(length(unique(as.integer(table(p$unit)))), 1L)  # really unbalanced

  fit_ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = p, outcome = "Y", group = "unit", time = "period",
      treatment = "D", effects = 3, placebo = 2, graph_off = TRUE
    )
  ))
  e_ref <- fit_ref$results$Effects[, 1]
  p_ref <- fit_ref$results$Placebos[, 1]

  for (bk in c("reference", "r")) {
    fit_us <- didgpu(
      df = p, outcome = "Y", group = "unit", time = "period",
      treatment = "D", effects = 3L, placebo = 2L,
      bootstrap_reps = 0L, backend = bk, verbose = FALSE
    )
    expect_lt(max(abs(fit_us$results$Effects[, "Estimate"] - e_ref)), 1e-12)
    expect_lt(max(abs(fit_us$results$Placebos[, "Estimate"] - p_ref)), 1e-12)
  }
})
