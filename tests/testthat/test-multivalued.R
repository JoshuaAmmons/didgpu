# Multivalued discrete treatment: D takes values in {0, 1, 2, ...} rather
# than just {0, 1}. The reference handles this naturally via per-d_sq
# cohort grouping. didgpu's r-backend uses the same grouping and so
# inherits the support without code changes.

build_multivalued_panel <- function(seed = 11L, n_units = 80L, n_periods = 15L,
                                      n_levels = 3L) {
  set.seed(seed)
  baseline <- sample(0L:(n_levels - 1L), n_units, replace = TRUE)
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), as.integer(n_units * 0.6)))
  F_g[treated] <- sample(5L:10L, length(treated), replace = TRUE)
  post_level <- baseline
  for (u in treated) {
    others <- setdiff(0L:(n_levels - 1L), baseline[u])
    post_level[u] <- sample(others, 1L)
  }
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$baseline   <- baseline[panel$unit]
  panel$F_g        <- F_g[panel$unit]
  panel$post_level <- post_level[panel$unit]
  panel$D <- ifelse(panel$period >= panel$F_g,
                    panel$post_level, panel$baseline)
  panel$k_evt <- panel$period - panel$F_g
  panel$tau_k <- 0
  post <- is.finite(panel$F_g) & panel$k_evt >= 0
  prof <- c(0.3, 0.5, 0.6, 0.5)
  panel$tau_k[post] <- (panel$post_level[post] - panel$baseline[post]) *
                       prof[pmin(panel$k_evt[post] + 1L, length(prof))]
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] +
             panel$tau_k + rnorm(nrow(panel), 0, 0.4)
  panel[order(panel$unit, panel$period), c("unit", "period", "D", "Y")]
}

test_that("multivalued (3 levels) matches reference bit-for-bit", {
  skip_if_no_reference()
  p <- build_multivalued_panel(seed = 11L, n_levels = 3L)
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 1, graph_off = TRUE
    )
  ))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 1L,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
  expect_lt(max(abs(as.numeric(us$results$Placebos[, "Estimate"]) -
                    as.numeric(ref$results$Placebos[, 1]))),
            8 * .Machine$double.eps)
})

test_that("multivalued (4 levels) matches reference bit-for-bit", {
  skip_if_no_reference()
  p <- build_multivalued_panel(seed = 42L, n_units = 120L, n_periods = 16L,
                                n_levels = 4L)
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE
    )
  ))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
})

test_that("multivalued with controls matches reference", {
  skip_if_no_reference()
  set.seed(11L)
  p <- build_multivalued_panel(seed = 11L, n_levels = 3L)
  # Add a control variable
  p$X <- rnorm(nrow(p), 0, 1) + 0.5 * p$D
  p$Y <- p$Y + 0.5 * p$X  # X affects Y too
  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(p), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE, controls = "X"
    )
  ))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, controls = "X",
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
})
