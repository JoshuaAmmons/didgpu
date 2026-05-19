# trends_lin = TRUE: allow group-specific linear trends.
# Estimator runs on first-difference of outcome; Effect_k = cumulative
# sum of per-event-time FD DIDs from j=1..k. Reference: did_multiplegt_main.R
# orchestrator at line 831-857, which calls the core repeatedly with
# effects = i for i in 1..n_effects, each call same_switchers=TRUE.

build_trend_panel <- function(seed = 7L) {
  set.seed(seed)
  n_units <- 100L; n_periods <- 15L
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(1:n_units, 60L))
  F_g[treated] <- sample(5L:10L, 60L, replace = TRUE)
  unit_fe <- rnorm(n_units, 0, 1)
  # group-specific linear trend: trends_lin must "absorb" this
  unit_trend <- rnorm(n_units, 0, 0.3)
  time_fe <- rnorm(n_periods, 0, 0.3)
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g)
  panel$k_evt <- panel$period - panel$F_g
  panel$tau_k <- 0
  prof <- c(0.5, 1.0, 1.2, 1.0)
  post <- is.finite(panel$F_g) & panel$k_evt >= 0
  panel$tau_k[post] <- prof[pmin(panel$k_evt[post] + 1L, length(prof))]
  panel$Y <- unit_fe[panel$unit] +
             unit_trend[panel$unit] * panel$period +
             time_fe[panel$period] +
             panel$tau_k + rnorm(nrow(panel), 0, 0.4)
  as.data.frame(panel[order(panel$unit, panel$period),
                       c("unit", "period", "D", "Y")])
}

test_that("trends_lin = TRUE matches reference bit-for-bit (effects = 3)", {
  skip_if_no_reference()
  p <- build_trend_panel()
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = p, outcome = "Y", group = "unit", time = "period", treatment = "D",
      effects = as.double(3), placebo = 0, graph_off = TRUE,
      trends_lin = TRUE)))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, trends_lin = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
})

test_that("trends_lin = TRUE matches reference at effects = 1, 2, 4", {
  skip_if_no_reference()
  p <- build_trend_panel()
  for (k in c(1L, 2L, 4L)) {
    ref <- suppressMessages(suppressWarnings(
      DIDmultiplegtDYN::did_multiplegt_dyn(
        df = p, outcome = "Y", group = "unit", time = "period", treatment = "D",
        effects = as.double(k), placebo = 0, graph_off = TRUE,
        trends_lin = TRUE)))
    us <- didgpu(p, "Y", "unit", "period", "D",
                  effects = k, placebo = 0L, trends_lin = TRUE,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE)
    expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                      as.numeric(ref$results$Effects[, 1]))),
              8 * .Machine$double.eps,
              label = sprintf("effects = %d", k))
  }
})

test_that("trends_lin = TRUE suppresses ATE (matches reference behaviour)", {
  p <- build_trend_panel()
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, trends_lin = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  # Reference returns no ATE (or an NA) under trends_lin = TRUE; we mirror
  # by setting it to NA.
  expect_true(is.na(as.numeric(us$results$ATE[1, "Estimate"])))
})

test_that("trends_lin = TRUE absorbs group-specific linear trends (strong-trend DGP)", {
  # Identifying-assumption sanity check: with trends_lin = TRUE the
  # recovered effect at k = 1 should be much closer to the true
  # tau_0 = 0.5 than without trends_lin, when the DGP carries large
  # group-specific linear trends. The trend SD is dialled up to 2.0
  # here (vs. the default 0.3 in build_trend_panel) so the bias
  # signal dominates sampling noise. Mean-squared over 5 seeds.
  trend_panel_strong <- function(seed) {
    set.seed(seed)
    n_units <- 100L; n_periods <- 15L
    F_g <- rep(Inf, n_units)
    treated <- sort(sample(1:n_units, 60L))
    F_g[treated] <- sample(5L:10L, 60L, replace = TRUE)
    unit_fe <- rnorm(n_units, 0, 1)
    unit_trend <- rnorm(n_units, 0, 2.0)  # MUCH bigger trends
    # Crucially: correlate trend with treatment timing so the linear
    # trend leaks into the vanilla DID estimate.
    unit_trend[treated] <- unit_trend[treated] + 1.0
    time_fe <- rnorm(n_periods, 0, 0.3)
    panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
    panel$F_g <- F_g[panel$unit]
    panel$D <- as.integer(panel$period >= panel$F_g)
    panel$k_evt <- panel$period - panel$F_g
    panel$tau_k <- 0
    prof <- c(0.5, 1.0, 1.2, 1.0)
    post <- is.finite(panel$F_g) & panel$k_evt >= 0
    panel$tau_k[post] <- prof[pmin(panel$k_evt[post] + 1L, length(prof))]
    panel$Y <- unit_fe[panel$unit] +
               unit_trend[panel$unit] * panel$period +
               time_fe[panel$period] +
               panel$tau_k + rnorm(nrow(panel), 0, 0.4)
    as.data.frame(panel[order(panel$unit, panel$period),
                         c("unit", "period", "D", "Y")])
  }
  seeds <- 1L:5L
  bias_n <- numeric(length(seeds))
  bias_y <- numeric(length(seeds))
  for (i in seq_along(seeds)) {
    p <- trend_panel_strong(seeds[i])
    us_n <- didgpu(p, "Y", "unit", "period", "D",
                    effects = 1L, bootstrap_reps = 0L,
                    backend = "r", verbose = FALSE)
    us_y <- didgpu(p, "Y", "unit", "period", "D",
                    effects = 1L, trends_lin = TRUE,
                    bootstrap_reps = 0L, backend = "r", verbose = FALSE)
    bias_n[i] <- as.numeric(us_n$results$Effects[1, "Estimate"]) - 0.5
    bias_y[i] <- as.numeric(us_y$results$Effects[1, "Estimate"]) - 0.5
  }
  expect_lt(mean(bias_y^2), mean(bias_n^2))
})

test_that("trends_lin + weight matches reference bit-for-bit", {
  skip_if_no_reference()
  p <- build_trend_panel()
  set.seed(99L)
  p$w <- runif(nrow(p), 0.5, 2.5)
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = p, outcome = "Y", group = "unit", time = "period", treatment = "D",
      effects = as.double(3), placebo = 0, graph_off = TRUE,
      trends_lin = TRUE, weight = "w")))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, trends_lin = TRUE, weight = "w",
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
})

test_that("trends_lin + controls matches reference bit-for-bit", {
  skip_if_no_reference()
  p <- build_trend_panel()
  set.seed(13L)
  p$X <- rnorm(nrow(p))
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = p, outcome = "Y", group = "unit", time = "period", treatment = "D",
      effects = as.double(3), placebo = 0, graph_off = TRUE,
      trends_lin = TRUE, controls = "X")))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, trends_lin = TRUE,
                controls = "X",
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            1e-8)
})

test_that("trends_lin + switchers = 'in' matches reference", {
  skip_if_no_reference()
  p <- build_trend_panel()
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = p, outcome = "Y", group = "unit", time = "period", treatment = "D",
      effects = as.double(3), placebo = 0, graph_off = TRUE,
      trends_lin = TRUE, switchers = "in")))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, trends_lin = TRUE,
                switchers = "in",
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
})

test_that("trends_lin sample-size columns match reference", {
  skip_if_no_reference()
  p <- build_trend_panel()
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = p, outcome = "Y", group = "unit", time = "period", treatment = "D",
      effects = as.double(3), placebo = 0, graph_off = TRUE,
      trends_lin = TRUE)))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, trends_lin = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  # N, Switchers, N.w, Switchers.w columns should match.
  for (col_name in c("N", "Switchers", "N.w", "Switchers.w")) {
    ref_col <- as.numeric(ref$results$Effects[, col_name])
    us_col  <- as.numeric(us$results$Effects[, col_name])
    expect_equal(us_col, ref_col,
                 label = sprintf("trends_lin %s column", col_name))
  }
})

test_that("trends_lin + normalized matches reference bit-for-bit", {
  skip_if_no_reference()
  p <- build_trend_panel()
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = p, outcome = "Y", group = "unit", time = "period", treatment = "D",
      effects = as.double(3), placebo = 0, graph_off = TRUE,
      trends_lin = TRUE, normalized = TRUE)))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, trends_lin = TRUE,
                normalized = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
})

test_that("trends_lin + normalized + placebos all bit-identical to reference", {
  skip_if_no_reference()
  # Late-switching DGP so placebos at k = 1, 2 fit.
  set.seed(7L)
  n_units <- 100L; n_periods <- 15L
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(1:n_units, 60L))
  F_g[treated] <- sample(6L:10L, 60L, replace = TRUE)
  unit_fe <- rnorm(n_units, 0, 1)
  unit_trend <- rnorm(n_units, 0, 0.3)
  time_fe <- rnorm(n_periods, 0, 0.3)
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g)
  panel$k_evt <- panel$period - panel$F_g
  panel$tau_k <- 0
  prof <- c(0.5, 1.0, 1.2, 1.0)
  post <- is.finite(panel$F_g) & panel$k_evt >= 0
  panel$tau_k[post] <- prof[pmin(panel$k_evt[post] + 1L, length(prof))]
  panel$Y <- unit_fe[panel$unit] + unit_trend[panel$unit] * panel$period +
             time_fe[panel$period] + panel$tau_k + rnorm(nrow(panel), 0, 0.4)
  p <- as.data.frame(panel[order(panel$unit, panel$period),
                            c("unit", "period", "D", "Y")])
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = p, outcome = "Y", group = "unit", time = "period", treatment = "D",
      effects = as.double(2), placebo = as.double(2), graph_off = TRUE,
      trends_lin = TRUE, normalized = TRUE)))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 2L, placebo = 2L,
                trends_lin = TRUE, normalized = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            1e-10)
  expect_lt(max(abs(as.numeric(us$results$Placebos[, "Estimate"]) -
                    as.numeric(ref$results$Placebos[, 1]))),
            1e-10)
})

test_that("trends_lin + placebos match reference at k = 1, 2, 3", {
  skip_if_no_reference()
  # Need F_g >= placebo + 2 so placebos at lag k exist; use late switches.
  set.seed(7L)
  n_units <- 100L; n_periods <- 15L
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(1:n_units, 60L))
  F_g[treated] <- sample(6L:10L, 60L, replace = TRUE)
  unit_fe <- rnorm(n_units, 0, 1)
  unit_trend <- rnorm(n_units, 0, 0.3)
  time_fe <- rnorm(n_periods, 0, 0.3)
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g)
  panel$k_evt <- panel$period - panel$F_g
  panel$tau_k <- 0
  prof <- c(0.5, 1.0, 1.2, 1.0)
  post <- is.finite(panel$F_g) & panel$k_evt >= 0
  panel$tau_k[post] <- prof[pmin(panel$k_evt[post] + 1L, length(prof))]
  panel$Y <- unit_fe[panel$unit] + unit_trend[panel$unit] * panel$period +
             time_fe[panel$period] + panel$tau_k + rnorm(nrow(panel), 0, 0.4)
  p <- as.data.frame(panel[order(panel$unit, panel$period),
                            c("unit", "period", "D", "Y")])
  for (npl in c(1L, 2L, 3L)) {
    ref <- suppressMessages(suppressWarnings(
      DIDmultiplegtDYN::did_multiplegt_dyn(
        df = p, outcome = "Y", group = "unit", time = "period", treatment = "D",
        effects = as.double(2), placebo = as.double(npl), graph_off = TRUE,
        trends_lin = TRUE)))
    us <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 2L, placebo = npl, trends_lin = TRUE,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE)
    expect_lt(max(abs(as.numeric(us$results$Placebos[, "Estimate"]) -
                      as.numeric(ref$results$Placebos[, 1]))),
              1e-10,
              label = sprintf("placebo = %d", npl))
  }
})

test_that("trends_lin = FALSE (default) is unchanged", {
  # Regression guard: passing trends_lin = FALSE explicitly must NOT
  # change existing behaviour.
  p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
                              frac_treated = 0.6,
                              min_treat_period = 4L, max_treat_period = 8L,
                              tau_profile = c(0.5, 1.0, 1.2),
                              sigma = 0.4, seed = 17L)
  us_default <- didgpu(p, "Y", "unit", "period", "D",
                       effects = 2L, bootstrap_reps = 0L,
                       backend = "r", verbose = FALSE)
  us_explicit <- didgpu(p, "Y", "unit", "period", "D",
                        effects = 2L, trends_lin = FALSE,
                        bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_identical(as.numeric(us_default$results$Effects[, "Estimate"]),
                   as.numeric(us_explicit$results$Effects[, "Estimate"]))
})
