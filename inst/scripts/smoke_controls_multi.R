# Stress-test controls: multiple controls, bidirectional, ATE.
library(didgpu); library(DIDmultiplegtDYN); library(data.table)

set.seed(11)
n_units <- 100L; n_periods <- 18L
F_g <- rep(Inf, n_units)
treated <- sort(sample(seq_len(n_units), 60L))
F_g[treated] <- sample(5L:12L, 60L, replace = TRUE)
unit_fe <- rnorm(n_units, 0, 1); time_fe <- rnorm(n_periods, 0, 0.3)

# Two controls
X1u <- rnorm(n_units, 0, 1); X1t <- rnorm(n_periods, 0, 0.5)
X2u <- rnorm(n_units, 0, 1); X2t <- rnorm(n_periods, 0, 0.5)

panel <- data.table(unit = rep(1:n_units, each = n_periods),
                     period = rep(1:n_periods, n_units))
panel[, F_g := F_g[unit]]
panel[, D := as.integer(period >= F_g & is.finite(F_g))]
panel[, X1 := X1u[unit] + X1t[period] + rnorm(.N, 0, 0.2)]
panel[, X2 := X2u[unit] + X2t[period] + rnorm(.N, 0, 0.2)]
panel[, k_evt := period - F_g]
tau <- c(0.5, 1.0, 1.2, 1.0, 0.8)
panel[, tau_k := 0]
panel[is.finite(F_g) & k_evt >= 0,
      tau_k := tau[pmin(k_evt + 1L, length(tau))]]
panel[, Y := unit_fe[unit] + time_fe[period] + tau_k +
            0.7 * X1 - 0.3 * X2 + rnorm(.N, 0, 0.4)]
panel_df <- as.data.frame(panel[order(unit, period),
                                  .(unit, period, D, Y, X1, X2)])

for (ctrl in list(c("X1"), c("X2"), c("X1", "X2"))) {
  cat(sprintf("\n=== controls = c(%s) ===\n",
              paste(shQuote(ctrl), collapse = ", ")))
  ref <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
    df = panel_df, outcome = "Y", group = "unit", time = "period",
    treatment = "D", effects = 4, placebo = 2, graph_off = TRUE,
    controls = ctrl)))
  us <- didgpu(panel_df, "Y", "unit", "period", "D",
                effects = 4L, placebo = 2L, controls = ctrl,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  cat("  effects max diff:",
      max(abs(as.numeric(ref$results$Effects[, 1]) -
              as.numeric(us$results$Effects[, "Estimate"]))), "\n")
  cat("  placebos max diff:",
      max(abs(as.numeric(ref$results$Placebos[, 1]) -
              as.numeric(us$results$Placebos[, "Estimate"]))), "\n")
  cat("  ATE diff:",
      abs(as.numeric(ref$results$ATE[1, 1]) -
          as.numeric(us$results$ATE[1, "Estimate"])), "\n")
}
