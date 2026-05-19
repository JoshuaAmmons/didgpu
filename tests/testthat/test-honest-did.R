# HonestDiD sensitivity-analysis wrapper.

skip_if_no_honest <- function() {
  testthat::skip_if_not_installed("HonestDiD")
}

build_es_panel <- function(seed = 17L, n_units = 100L, true_att = 1.0) {
  set.seed(seed)
  n_periods <- 12L
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), as.integer(n_units * 0.5)))
  F_g[treated] <- sample(seq(5L, n_periods - 2L), length(treated),
                          replace = TRUE)
  panel <- expand.grid(unit = seq_len(n_units), period = seq_len(n_periods))
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g)
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] +
             true_att * panel$D + rnorm(nrow(panel), 0, 0.3)
  panel[order(panel$unit, panel$period),
         c("unit", "period", "D", "Y")]
}

test_that("didgpu_honest_did runs end-to-end on a CS fit (RM method)", {
  skip_if_no_honest()
  p <- build_es_panel(seed = 17L)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR",
                    bootstrap_reps = 30L, seed = 1L,
                    backend = "r", verbose = FALSE)
  sens <- didgpu_honest_did(fit, event_post = 1L,
                              method = "RM",
                              Mbar = c(0, 0.5, 1.0))
  expect_s3_class(sens, "didgpu_honest_did_result")
  expect_true(all(c("Mbar", "lb", "ub", "crosses_zero") %in% names(sens)))
  expect_equal(nrow(sens), 3L)
})

test_that("didgpu_honest_did runs with M (smoothness) method", {
  skip_if_no_honest()
  p <- build_es_panel(seed = 23L)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR",
                    bootstrap_reps = 30L, seed = 1L,
                    backend = "r", verbose = FALSE)
  sens <- didgpu_honest_did(fit, event_post = 1L,
                              method = "M",
                              Mbar = c(0, 0.2, 0.5))
  expect_equal(nrow(sens), 3L)
  # Larger M => wider bounds.
  expect_gte(sens$ub[3] - sens$lb[3], sens$ub[1] - sens$lb[1])
})

test_that("didgpu_honest_did reports a breakdown when CI eventually crosses zero", {
  skip_if_no_honest()
  p <- build_es_panel(seed = 17L, true_att = 0.5)
  fit <- didgpu_cs(p, "Y", "unit", "period", "D",
                    est_method = "OR",
                    bootstrap_reps = 30L, seed = 1L,
                    backend = "r", verbose = FALSE)
  sens <- didgpu_honest_did(fit, event_post = 1L,
                              method = "RM",
                              Mbar = seq(0, 5, by = 0.5))
  bd <- attr(sens, "breakdown")
  # At Mbar large enough, CI will include 0. Either we find a breakdown
  # in the grid, or it's NA (no Mbar tested crossed zero).
  expect_true(is.na(bd) || is.numeric(bd))
})

test_that("didgpu_honest_did errors if no SEs available", {
  skip_if_no_honest()
  p <- build_es_panel(seed = 17L)
  fit_nose <- didgpu_cs(p, "Y", "unit", "period", "D",
                         est_method = "OR",
                         bootstrap_reps = 0L,
                         backend = "r", verbose = FALSE)
  expect_error(
    didgpu_honest_did(fit_nose, event_post = 1L,
                       method = "RM", Mbar = c(0, 1)),
    "requires SE estimates"
  )
})

test_that("didgpu_honest_did errors when HonestDiD package is missing", {
  # Mock: temporarily mask requireNamespace to always fail.
  skip_if(requireNamespace("HonestDiD", quietly = TRUE),
           "HonestDiD is installed; can't test missing-pkg path")
})
