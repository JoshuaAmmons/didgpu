# Continuous treatment: D is genuinely continuous (e.g., dosage). The
# reference handles this by binarizing treatment to a +/-1 indicator,
# collapsing d_sq cohorts, and adding (time-onset >= j) * baseline^k
# interaction controls for j = 2..T_max and k = 1..continuous.
# didgpu's r-backend replicates this.

build_continuous_panel <- function(seed = 11L, n_units = 80L, n_periods = 15L) {
  set.seed(seed)
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), as.integer(n_units * 0.6)))
  F_g[treated] <- sample(5L:10L, length(treated), replace = TRUE)
  baseline_d <- runif(n_units, 0, 1)
  post_d <- baseline_d + rnorm(n_units, 0, 0.5)
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$F_g <- F_g[panel$unit]
  panel$D <- ifelse(panel$period >= panel$F_g,
                    post_d[panel$unit], baseline_d[panel$unit])
  panel$k_evt <- panel$period - panel$F_g
  panel$tau_k <- 0
  post <- is.finite(panel$F_g) & panel$k_evt >= 0
  prof <- c(0.3, 0.5, 0.6, 0.5)
  panel$tau_k[post] <- (post_d[panel$unit[post]] - baseline_d[panel$unit[post]]) *
                       prof[pmin(panel$k_evt[post] + 1L, length(prof))]
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] +
             panel$tau_k + rnorm(nrow(panel), 0, 0.4)
  panel[order(panel$unit, panel$period), c("unit", "period", "D", "Y")]
}

test_that("continuous = 1 matches reference bit-for-bit", {
  skip_if_no_reference()
  p <- build_continuous_panel(seed = 11L)
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE, continuous = 1
    )
  ))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, continuous = 1L,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
})

test_that("continuous = 2 (quadratic) matches reference bit-for-bit", {
  skip_if_no_reference()
  p <- build_continuous_panel(seed = 11L)
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE, continuous = 2
    )
  ))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, continuous = 2L,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
})

test_that("continuous returns NA on a panel without the flag (degenerate cohorts)", {
  # Without the continuous= flag, a continuous-D panel has each unit in
  # its own d_sq cohort, leaving no controls. didgpu returns NA cleanly.
  p <- build_continuous_panel(seed = 11L)
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 1L, bootstrap_reps = 0L,
                backend = "r", verbose = FALSE)
  expect_true(all(is.na(us$results$Effects[, "Estimate"])))
})
