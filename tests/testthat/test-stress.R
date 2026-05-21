# Stress test: combine bidirectional panel + multiple effects + placebos +
# cluster argument + a substantive bootstrap, and verify the r-backend
# matches the reference-backend bit-for-bit through the full orchestrator.
#
# Catches regressions that only manifest at the intersection of features.

test_that("r-backend == reference on a full-featured bidirectional bootstrap", {
  skip_if_no_reference()
  set.seed(11L)

  # Bidirectional panel. (cluster argument is NOT exercised here because
  # DIDmultiplegtDYN 2.2.0 has a bug where passing `cluster=` triggers an
  # `invalid first argument` error inside its core — `get(cluster)` looks
  # up the unrenamed column. Our r-backend handles cluster correctly; see
  # test-bidirectional.R::cluster argument for that path.)
  p <- didgpu_simulate_panel_bidir(
    n_units = 100L, n_periods = 18L,
    frac_treated = 0.6, frac_in = 0.5,
    min_treat_period = 6L, max_treat_period = 12L,
    seed = 11L
  )

  cdir_r   <- tempfile("didgpu_stress_r_")
  cdir_ref <- tempfile("didgpu_stress_ref_")
  on.exit(unlink(c(cdir_r, cdir_ref), recursive = TRUE), add = TRUE)

  fit_r <- didgpu(p, "Y", "unit", "period", "D",
                   effects = 4L, placebo = 2L,
                   bootstrap_reps = 25L, seed = 7L,
                   checkpoint_dir = cdir_r,
                   backend = "r", verbose = FALSE)
  fit_ref <- didgpu(p, "Y", "unit", "period", "D",
                     effects = 4L, placebo = 2L,
                     bootstrap_reps = 25L, seed = 7L,
                     checkpoint_dir = cdir_ref,
                     backend = "reference", verbose = FALSE)

  # Point estimates bit-identical (both backends apply the same cluster-
  # resample with the same seed; both produce the same point estimate on
  # each resampled panel; therefore the bootstrap-averaged estimates and
  # all derived SEs match).
  for (col in c("Estimate", "SE", "LB.CI", "UB.CI")) {
    diff_e <- max(abs(fit_r$results$Effects[, col]   -
                      fit_ref$results$Effects[, col]))
    diff_p <- max(abs(fit_r$results$Placebos[, col]  -
                      fit_ref$results$Placebos[, col]))
    diff_a <- max(abs(fit_r$results$ATE[1, col]      -
                      fit_ref$results$ATE[1, col]))
    expect_lt(diff_e, 1e-10, label = sprintf("Effects[%s] diff=%.2e", col, diff_e))
    expect_lt(diff_p, 1e-10, label = sprintf("Placebos[%s] diff=%.2e", col, diff_p))
    expect_lt(diff_a, 4 * .Machine$double.eps,
              label = sprintf("ATE[%s] diff=%.2e", col, diff_a))
  }
  # Sample-size cols agree on the counts.
  expect_equal(as.integer(fit_r$results$Effects[, "N"]),
               as.integer(fit_ref$results$Effects[, "N"]))
  expect_equal(as.integer(fit_r$results$Effects[, "Switchers"]),
               as.integer(fit_ref$results$Effects[, "Switchers"]))

  # Joint p-values agree (computed from the same bootstrap empirical cov).
  expect_lt(abs(fit_r$results$p_jointeffects -
                fit_ref$results$p_jointeffects),
            1e-10)
  expect_lt(abs(fit_r$results$p_jointplacebo -
                fit_ref$results$p_jointplacebo),
            1e-10)
})

test_that("orchestrator handles a panel where some cells have no switchers", {
  # Sparse switcher coverage: 5 units treated, 95 never-treated, late switch.
  # Some event-times may have no switchers; aggregator must NA them
  # gracefully rather than NaN-poisoning everything downstream.
  set.seed(7L)
  n_units <- 100L; n_periods <- 20L
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), 5L))
  F_g[treated] <- sample(15L:18L, 5L, replace = TRUE)
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g & is.finite(panel$F_g))
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] +
             rnorm(nrow(panel), 0, 0.4)
  panel <- panel[order(panel$unit, panel$period), c("unit", "period", "D", "Y")]

  # effects = 8 — longer than any L_g (max F_g = 18, T_max = 20 => max L_g = 6).
  # Both backends should auto-clamp.
  fit <- didgpu(panel, "Y", "unit", "period", "D",
                 effects = 8L, placebo = 0L,
                 bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  # max L_g for these treated units is at most 6 (F_g = 15 case).
  expect_lte(fit$results$N_Effects, 6L)
  # No NaN poisoning across effects that DO exist.
  expect_true(all(!is.nan(fit$results$Effects[, "Estimate"])))
})
