test_that("weight arg matches reference bit-for-bit", {
  skip_if_no_reference()
  set.seed(11L)
  n_units <- 80L; n_periods <- 15L
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), 48L))
  F_g[treated] <- sample(5L:10L, 48L, replace = TRUE)
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g & is.finite(panel$F_g))
  panel$k_evt <- panel$period - panel$F_g
  panel$tau_k <- 0
  post <- is.finite(panel$F_g) & panel$k_evt >= 0
  panel$tau_k[post] <- c(0.5, 1.0, 1.2)[pmin(panel$k_evt[post] + 1L, 3L)]
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] + panel$tau_k +
             rnorm(nrow(panel), 0, 0.4)
  panel$w <- runif(nrow(panel), 0.5, 2.0)
  panel <- panel[order(panel$unit, panel$period),
                 c("unit", "period", "D", "Y", "w")]

  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(panel), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE, weight = "w"
    )
  ))
  us <- didgpu(panel, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, weight = "w",
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)

  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))), 1e-10)
})

test_that("weight arg errors clearly when column not in df", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 5L)
  expect_error(
    didgpu(p, "Y", "unit", "period", "D",
            effects = 1L, weight = "nonexistent",
            bootstrap_reps = 0L, backend = "r", verbose = FALSE),
    "not in df"
  )
})

test_that("weight=NULL produces same result as default (all 1s)", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 5L,
                              min_treat_period = 3L, max_treat_period = 7L)
  fit_no <- didgpu(p, "Y", "unit", "period", "D",
                    effects = 2L, bootstrap_reps = 0L,
                    backend = "r", verbose = FALSE)
  p$w <- 1
  fit_w <- didgpu(p, "Y", "unit", "period", "D",
                   effects = 2L, weight = "w", bootstrap_reps = 0L,
                   backend = "r", verbose = FALSE)
  expect_equal(fit_no$results$Effects[, "Estimate"],
               fit_w$results$Effects[, "Estimate"])
})
