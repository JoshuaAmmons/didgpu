# Callaway-Sant'Anna (2021) staggered DiD estimator.
# v1: OR estimator, never-treated controls, no covariates.

test_that("didgpu_cs returns a result with the expected shape", {
  p <- didgpu_simulate_panel(n_units = 100L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              seed = 17L)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR", aggregation = "event",
                    bootstrap_reps = 0L,
                    backend = "r", verbose = FALSE)
  expect_s3_class(fit, "didgpu_cs_result")
  # ATT(g, t) table: one row per (cohort, post-treatment time).
  expect_true(all(c("g", "t", "event_time", "att",
                     "n_treated", "n_control") %in% names(fit$att_gt)))
  expect_true(nrow(fit$att_gt) > 0L)
  # Event-study aggregation present.
  expect_true("event_time" %in% names(fit$aggregation))
})

test_that("didgpu_cs OR recovers approximately the true ATT", {
  # Simulate a panel with constant ATT = 1.0; OR should recover ~1.0 at
  # every post-treatment event-time within sampling error.
  set.seed(7L)
  n_units <- 200L; n_periods <- 12L
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), 100L))
  F_g[treated] <- sample(seq(3L, n_periods - 2L), 100L, replace = TRUE)
  panel <- expand.grid(unit = seq_len(n_units), period = seq_len(n_periods))
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g)
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] +
             1.0 * panel$D + rnorm(nrow(panel), 0, 0.3)
  p <- panel[order(panel$unit, panel$period),
              c("unit", "period", "D", "Y")]
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR", aggregation = "event",
                    bootstrap_reps = 0L,
                    backend = "r", verbose = FALSE)
  # Average POST-TREATMENT ATT across all (g, t) cells should be near 1.0
  # (pre-treatment placebos with event_time < 0 are ~0 and would drag
  # the overall mean toward zero).
  post <- fit$att_gt[fit$att_gt$event_time >= 0L, ]
  mean_att <- mean(post$att, na.rm = TRUE)
  expect_lt(abs(mean_att - 1.0), 0.2)
})

test_that("didgpu_cs event-study aggregation produces one row per event-time", {
  p <- didgpu_simulate_panel(n_units = 100L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR", aggregation = "event",
                    bootstrap_reps = 0L,
                    backend = "r", verbose = FALSE)
  ev <- fit$aggregation
  expect_equal(length(unique(ev$event_time)), nrow(ev))
  # All event-times present in the att_gt table appear in the aggregation.
  expect_setequal(ev$event_time, unique(fit$att_gt$event_time))
})

test_that("didgpu_cs_aggregate re-aggregates without refitting", {
  p <- didgpu_simulate_panel(n_units = 100L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR", aggregation = "event",
                    bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  # Switch to group, calendar, overall — same att_gt, different summary.
  fit_g <- didgpu_cs_aggregate(fit, "group")
  expect_true("g" %in% names(fit_g$aggregation))
  fit_c <- didgpu_cs_aggregate(fit, "calendar")
  expect_true("t" %in% names(fit_c$aggregation))
  fit_o <- didgpu_cs_aggregate(fit, "overall")
  expect_true(nrow(fit_o$aggregation) == 1L)
  expect_true(is.finite(fit_o$aggregation$estimate))
})

test_that("didgpu_cs with bootstrap fills in SE and CI columns", {
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 10L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR",
                    bootstrap_reps = 5L, seed = 1L,
                    backend = "r", verbose = FALSE)
  expect_true("se" %in% names(fit$att_gt))
  expect_true("ci_low" %in% names(fit$att_gt))
  expect_true("ci_high" %in% names(fit$att_gt))
  # At least one (g, t) cell has a finite SE.
  expect_true(any(is.finite(fit$att_gt$se)))
})

test_that("didgpu_cs validates required args", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 8L, seed = 1L)
  expect_error(
    didgpu_cs(p, "nope", "unit", "period", "D"),
    "column not in df: nope"
  )
  expect_error(
    didgpu_cs(p, "Y", "unit", "period", "D", est_method = "BAD"),
    "should be one of"
  )
})

test_that("didgpu_cs IPW + DR + notyet + covariates now run successfully", {
  set.seed(1L)
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 10L, seed = 1L)
  # Add a unit-level covariate.
  p$x1 <- rnorm(length(unique(p$unit)))[p$unit]
  for (combo in list(
    list(est_method = "DR"),
    list(est_method = "IPW"),
    list(control_group = "notyet"),
    list(est_method = "DR", covariates = "x1"),
    list(est_method = "IPW", covariates = "x1"),
    list(est_method = "OR",  covariates = "x1")
  )) {
    fit <- do.call(didgpu_cs,
                    c(list(p, "Y", "unit", "period", "D",
                            bootstrap_reps = 0L, backend = "r",
                            verbose = FALSE), combo))
    expect_s3_class(fit, "didgpu_cs_result")
    expect_true(nrow(fit$att_gt) > 0L)
    expect_true(all(is.finite(fit$att_gt$att)) ||
                 sum(is.na(fit$att_gt$att)) < nrow(fit$att_gt))
  }
})

test_that("didgpu_cs reports pre-treatment placebos by default", {
  p <- didgpu_simulate_panel(n_units = 100L, n_periods = 12L,
                              tau_profile = c(0.5, 1.0),
                              seed = 17L)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR", aggregation = "event",
                    bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_true(any(fit$att_gt$event_time < 0L))
  expect_true(!is.null(fit$placebo))
  expect_true("per_event" %in% names(fit$placebo))
  # Pre-treatment placebos should be small (clean panel, no anticipation).
  pre_ev <- fit$placebo$per_event
  expect_true(all(abs(pre_ev$estimate) < 0.5))
})

test_that("didgpu_cs multiplier bootstrap runs and produces finite SEs", {
  p <- didgpu_simulate_panel(n_units = 80L, n_periods = 10L,
                              tau_profile = c(0.5, 1.0), seed = 17L)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR",
                    bootstrap_reps = 50L,
                    bootstrap_kind = "multiplier",
                    seed = 1L,
                    backend = "r", verbose = FALSE)
  expect_true(any(is.finite(fit$att_gt$se)))
  expect_true(all(fit$att_gt$se >= 0 | is.na(fit$att_gt$se)))
})

test_that("didgpu_cs notyet uses not-yet-treated units as controls", {
  # On a panel WITHOUT any never-treated unit, "never" errors out but
  # "notyet" should still work (uses not-yet-treated as controls).
  set.seed(1L)
  n_units <- 60L; n_periods <- 12L
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  # Every unit eventually gets treated; F_g varies in [4..10].
  F_g <- sample(4L:10L, n_units, replace = TRUE)
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g)
  panel$Y <- rnorm(n_units)[panel$unit] +
             rnorm(n_periods, 0, 0.3)[panel$period] +
             0.5 * panel$D + rnorm(nrow(panel), 0, 0.3)
  p <- panel[order(panel$unit, panel$period),
              c("unit", "period", "D", "Y")]
  expect_error(
    didgpu_cs(p, "Y", "unit", "period", "D",
               est_method = "OR", control_group = "never",
               bootstrap_reps = 0L, backend = "r", verbose = FALSE),
    "never-treated"
  )
  # notyet should succeed.
  fit_ny <- didgpu_cs(p, "Y", "unit", "period", "D",
                       est_method = "OR", control_group = "notyet",
                       bootstrap_reps = 0L, backend = "r",
                       verbose = FALSE)
  expect_s3_class(fit_ny, "didgpu_cs_result")
})

test_that("didgpu_cs needs at least one never-treated unit", {
  set.seed(1L)
  # Panel where every unit is treated by period 5 (no never-treated).
  panel <- expand.grid(unit = 1:30, period = 1:10)
  panel$F_g <- 5L
  panel$D <- as.integer(panel$period >= panel$F_g)
  panel$Y <- rnorm(nrow(panel))
  p <- panel[order(panel$unit, panel$period),
              c("unit", "period", "D", "Y")]
  expect_error(
    didgpu_cs(p, "Y", "unit", "period", "D",
               est_method = "OR", bootstrap_reps = 0L,
               backend = "r", verbose = FALSE),
    "never-treated"
  )
})
