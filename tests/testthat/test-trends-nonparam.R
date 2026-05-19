test_that("trends_nonparam matches reference bit-for-bit", {
  skip_if_no_reference()
  set.seed(11L)
  n_units <- 100L; n_periods <- 15L
  F_g <- rep(Inf, n_units)
  treated <- sort(sample(seq_len(n_units), 60L))
  F_g[treated] <- sample(5L:10L, 60L, replace = TRUE)
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  industry <- sample.int(3L, n_units, replace = TRUE)
  industry_trend <- matrix(rnorm(3L * n_periods, 0, 0.4),
                            nrow = 3L, ncol = n_periods)
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  panel$industry <- industry[panel$unit]
  panel$F_g <- F_g[panel$unit]
  panel$D <- as.integer(panel$period >= panel$F_g & is.finite(panel$F_g))
  panel$k_evt <- panel$period - panel$F_g
  panel$tau_k <- 0
  post <- is.finite(panel$F_g) & panel$k_evt >= 0
  panel$tau_k[post] <- c(0.5, 1.0, 1.2, 1.0)[pmin(panel$k_evt[post] + 1L, 4L)]
  panel$ind_time <- industry_trend[cbind(panel$industry, panel$period)]
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] + panel$tau_k +
             panel$ind_time + rnorm(nrow(panel), 0, 0.4)
  panel <- panel[order(panel$unit, panel$period),
                 c("unit", "period", "D", "Y", "industry")]

  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(panel), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 3, placebo = 0, graph_off = TRUE,
      trends_nonparam = "industry"
    )
  ))
  us <- didgpu(panel, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, trends_nonparam = "industry",
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)

  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))), 1e-10)
})

test_that("trends_nonparam errors when column not in df", {
  p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L, seed = 5L)
  expect_error(
    didgpu(p, "Y", "unit", "period", "D",
            effects = 1L, trends_nonparam = "nonexistent",
            bootstrap_reps = 0L, backend = "r", verbose = FALSE),
    "not in df"
  )
})
