# Controls support (Phase 1: point estimates, no variance correction).
# The r-backend matches the reference within ~0.5% on the supported
# subset. Bit-identical match requires porting a subtle weighting
# detail that's still being investigated; the current implementation
# is documented in R/controls.R.

test_that("controls arg matches reference bit-for-bit", {
  skip_if_no_reference()
  set.seed(11L)
  n_units <- 80L; n_periods <- 15L
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), 48L))
  F_g[treated] <- sample(5L:10L, 48L, replace = TRUE)
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  X_unit  <- rnorm(n_units, 0, 1)
  X_time  <- rnorm(n_periods, 0, 0.5)

  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g & is.finite(panel$F_g))
  panel$X <- X_unit[panel$unit] + X_time[panel$period] +
             rnorm(nrow(panel), 0, 0.2)
  panel$k_evt <- panel$period - panel$F_g
  panel$tau_k <- 0
  post <- is.finite(panel$F_g) & panel$k_evt >= 0
  panel$tau_k[post] <- c(0.5, 1.0, 1.2, 1.0)[pmin(panel$k_evt[post] + 1L, 4L)]
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] +
             panel$tau_k + 0.7 * panel$X + rnorm(nrow(panel), 0, 0.4)
  panel <- panel[order(panel$unit, panel$period), c("unit", "period", "D", "Y", "X")]

  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(panel), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE, controls = "X"
    )
  ))
  us <- didgpu(panel, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, controls = "X",
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)

  ref_e <- as.numeric(ref$results$Effects[, 1])
  us_e  <- as.numeric(us$results$Effects[, "Estimate"])

  # Bit-identical match (within machine epsilon).
  expect_lt(max(abs(ref_e - us_e)), 1e-10,
            label = sprintf("controls max abs diff=%.2e",
                            max(abs(ref_e - us_e))))
})

test_that("multi-control panel matches reference bit-for-bit", {
  skip_if_no_reference()
  set.seed(11L)
  n_units <- 100L; n_periods <- 18L
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), 60L))
  F_g[treated] <- sample(5L:12L, 60L, replace = TRUE)
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  X1u <- rnorm(n_units, 0, 1); X1t <- rnorm(n_periods, 0, 0.5)
  X2u <- rnorm(n_units, 0, 1); X2t <- rnorm(n_periods, 0, 0.5)
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g & is.finite(panel$F_g))
  panel$X1 <- X1u[panel$unit] + X1t[panel$period] + rnorm(nrow(panel), 0, 0.2)
  panel$X2 <- X2u[panel$unit] + X2t[panel$period] + rnorm(nrow(panel), 0, 0.2)
  panel$k_evt <- panel$period - panel$F_g
  panel$tau_k <- 0
  post <- is.finite(panel$F_g) & panel$k_evt >= 0
  panel$tau_k[post] <- c(0.5, 1.0, 1.2, 1.0, 0.8)[pmin(panel$k_evt[post] + 1L, 5L)]
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] + panel$tau_k +
             0.7 * panel$X1 - 0.3 * panel$X2 + rnorm(nrow(panel), 0, 0.4)
  panel <- panel[order(panel$unit, panel$period),
                 c("unit", "period", "D", "Y", "X1", "X2")]

  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(panel), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 4, placebo = 2, graph_off = TRUE,
      controls = c("X1", "X2")
    )
  ))
  us <- didgpu(panel, "Y", "unit", "period", "D",
                effects = 4L, placebo = 2L, controls = c("X1", "X2"),
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)

  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))), 1e-10)
  expect_lt(max(abs(as.numeric(us$results$Placebos[, "Estimate"]) -
                    as.numeric(ref$results$Placebos[, 1]))), 1e-10)
  expect_lt(abs(as.numeric(us$results$ATE[1, "Estimate"]) -
                as.numeric(ref$results$ATE[1, 1])), 4 * .Machine$double.eps)
})

test_that("controls + placebos match reference bit-for-bit", {
  skip_if_no_reference()
  set.seed(11L)
  n_units <- 80L; n_periods <- 15L
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), 48L))
  F_g[treated] <- sample(5L:10L, 48L, replace = TRUE)
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  X_unit  <- rnorm(n_units, 0, 1)
  X_time  <- rnorm(n_periods, 0, 0.5)
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g & is.finite(panel$F_g))
  panel$X <- X_unit[panel$unit] + X_time[panel$period] + rnorm(nrow(panel), 0, 0.2)
  panel$k_evt <- panel$period - panel$F_g
  panel$tau_k <- 0
  post <- is.finite(panel$F_g) & panel$k_evt >= 0
  panel$tau_k[post] <- c(0.5, 1.0, 1.2, 1.0)[pmin(panel$k_evt[post] + 1L, 4L)]
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] + panel$tau_k +
             0.7 * panel$X + rnorm(nrow(panel), 0, 0.4)
  panel <- panel[order(panel$unit, panel$period),
                 c("unit", "period", "D", "Y", "X")]

  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(panel), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 1, graph_off = TRUE, controls = "X"
    )
  ))
  us <- didgpu(panel, "Y", "unit", "period", "D",
                effects = 3L, placebo = 1L, controls = "X",
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)

  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))), 1e-10)
  expect_lt(max(abs(as.numeric(us$results$Placebos[, "Estimate"]) -
                    as.numeric(ref$results$Placebos[, 1]))), 1e-10)
})

test_that("controls arg degrades gracefully to no-controls when controls=NULL", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 5L,
                              min_treat_period = 3L, max_treat_period = 5L)
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 2L, placebo = 0L, controls = NULL,
                 bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_equal(fit$results$N_Effects, 2L)
})

test_that("controls arg errors clearly when column not in df", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 5L)
  expect_error(
    didgpu(p, "Y", "unit", "period", "D",
            effects = 1L, controls = "nonexistent_col",
            bootstrap_reps = 0L, backend = "r", verbose = FALSE),
    "not in df"
  )
})
