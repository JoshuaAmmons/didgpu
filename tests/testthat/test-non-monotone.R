# Non-monotone treatment paths: some units have treatment that BOTH
# increases above and decreases below the baseline. By default the
# reference drops the rows AFTER the second-direction switch. didgpu's
# .prep_panel now replicates this behavior.

test_that("non-monotone panel matches reference (default drop)", {
  skip_if_no_reference()
  set.seed(11L)
  # 50 units, 12 periods. Some have monotone increases, some have
  # monotone decreases, some have 0->1->0 (non-monotone — gets dropped
  # post-second-switch by both backends).
  n_units <- 50L; n_periods <- 12L
  panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
  unit_fe <- rnorm(n_units, 0, 1)
  time_fe <- rnorm(n_periods, 0, 0.3)
  panel$D <- 0L
  for (u in 1:n_units) {
    # 1/3 monotone increase, 1/3 monotone decrease, 1/6 non-monotone, 1/6 always 0
    typ <- u %% 6L
    if (typ == 0L) {
      panel$D[panel$unit == u] <- 0L
    } else if (typ == 1L || typ == 2L) {
      # Increase at period 5
      panel$D[panel$unit == u] <- as.integer(panel$period[panel$unit == u] >= 5L)
    } else if (typ == 3L || typ == 4L) {
      # Decrease at period 5 (starts at 1)
      panel$D[panel$unit == u] <- as.integer(panel$period[panel$unit == u] < 5L)
    } else {
      # Non-monotone: 0,0,1,1,1,0,0,...
      d <- ifelse(panel$period[panel$unit == u] >= 3L &
                  panel$period[panel$unit == u] <= 5L, 1L, 0L)
      panel$D[panel$unit == u] <- d
    }
  }
  panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] +
             0.5 * panel$D + rnorm(nrow(panel), 0, 0.4)
  panel <- panel[order(panel$unit, panel$period), c("unit", "period", "D", "Y")]

  ref <- suppressMessages(suppressWarnings(
    DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(panel), outcome = "Y", group = "unit",
      time = "period", treatment = "D",
      effects = 2, placebo = 0, graph_off = TRUE
    )
  ))
  us <- didgpu(panel, "Y", "unit", "period", "D",
                effects = 2L, placebo = 0L,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)

  expect_lt(max(abs(as.numeric(us$results$Effects[, "Estimate"]) -
                    as.numeric(ref$results$Effects[, 1]))),
            8 * .Machine$double.eps)
})
