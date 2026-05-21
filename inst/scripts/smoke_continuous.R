library(didgpu); library(DIDmultiplegtDYN); library(data.table)

set.seed(11)
n_units <- 80L; n_periods <- 15L
F_g <- rep(Inf, n_units)
treated <- sort(sample(seq_len(n_units), 48L))
F_g[treated] <- sample(5L:10L, 48L, replace = TRUE)
baseline_d <- runif(n_units, 0, 1)
post_d <- baseline_d + rnorm(n_units, 0, 0.5)
unit_fe <- rnorm(n_units, 0, 1); time_fe <- rnorm(n_periods, 0, 0.3)

panel <- data.table(unit = rep(1:n_units, each = n_periods),
                     period = rep(1:n_periods, n_units))
panel[, F_g := F_g[unit]]
panel[, D := ifelse(period >= F_g, post_d[unit], baseline_d[unit])]
panel[, k_evt := period - F_g]
panel[, tau_k := 0]
panel[is.finite(F_g) & k_evt >= 0,
      tau_k := (post_d[unit] - baseline_d[unit]) *
        c(0.3, 0.5, 0.6, 0.5)[pmin(k_evt + 1L, 4L)]]
panel[, Y := unit_fe[unit] + time_fe[period] + tau_k + rnorm(.N, 0, 0.4)]
panel_df <- as.data.frame(panel[order(unit, period), .(unit, period, D, Y)])

cat("--- ref WITH continuous = 1 ---\n")
ref <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = panel_df, outcome = "Y", group = "unit", time = "period",
  treatment = "D", effects = 3, placebo = 0, graph_off = TRUE,
  continuous = 1)))
print(as.numeric(ref$results$Effects[, 1]))

cat("\n--- didgpu WITH continuous = 1 ---\n")
us <- didgpu(panel_df, "Y", "unit", "period", "D",
              effects = 3L, placebo = 0L, continuous = 1L,
              bootstrap_reps = 0L, backend = "r", verbose = FALSE)
print(as.numeric(us$results$Effects[, "Estimate"]))

cat("\nmax diff:",
    max(abs(as.numeric(ref$results$Effects[, 1]) -
            as.numeric(us$results$Effects[, "Estimate"]))), "\n")

cat("\n--- ref WITH continuous = 2 (quadratic) ---\n")
ref2 <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = panel_df, outcome = "Y", group = "unit", time = "period",
  treatment = "D", effects = 3, placebo = 0, graph_off = TRUE,
  continuous = 2)))
print(as.numeric(ref2$results$Effects[, 1]))

cat("\n--- didgpu WITH continuous = 2 ---\n")
us2 <- didgpu(panel_df, "Y", "unit", "period", "D",
               effects = 3L, placebo = 0L, continuous = 2L,
               bootstrap_reps = 0L, backend = "r", verbose = FALSE)
print(as.numeric(us2$results$Effects[, "Estimate"]))

cat("\nmax diff (continuous=2):",
    max(abs(as.numeric(ref2$results$Effects[, 1]) -
            as.numeric(us2$results$Effects[, "Estimate"]))), "\n")
