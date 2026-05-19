library(didgpu); library(DIDmultiplegtDYN); library(data.table)

set.seed(11)
n_units <- 80L; n_periods <- 15L
F_g <- rep(Inf, n_units)
treated <- sort(sample(seq_len(n_units), 48L))
F_g[treated] <- sample(5L:10L, 48L, replace = TRUE)
unit_fe <- rnorm(n_units, 0, 1); time_fe <- rnorm(n_periods, 0, 0.3)
panel <- data.table(unit = rep(1:n_units, each = n_periods),
                     period = rep(1:n_periods, n_units))
panel[, F_g := F_g[unit]]
panel[, D := as.integer(period >= F_g & is.finite(F_g))]
panel[, k_evt := period - F_g]
tau <- c(0.5, 1.0, 1.2)
panel[, tau_k := 0]
panel[is.finite(F_g) & k_evt >= 0, tau_k := tau[pmin(k_evt + 1L, length(tau))]]
panel[, Y := unit_fe[unit] + time_fe[period] + tau_k + rnorm(.N, 0, 0.4)]
# Add a non-trivial weight: random positive numbers.
panel[, w := runif(.N, 0.5, 2.0)]
panel_df <- as.data.frame(panel[order(unit, period), .(unit, period, D, Y, w)])

cat("--- ref without weight ---\n")
ref_nw <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = panel_df, outcome = "Y", group = "unit", time = "period",
  treatment = "D", effects = 3, placebo = 0, graph_off = TRUE)))
print(as.numeric(ref_nw$results$Effects[, 1]))

cat("\n--- ref with weight = 'w' ---\n")
ref_w <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = panel_df, outcome = "Y", group = "unit", time = "period",
  treatment = "D", effects = 3, placebo = 0, graph_off = TRUE,
  weight = "w")))
print(as.numeric(ref_w$results$Effects[, 1]))

cat("\n--- didgpu with weight = 'w' ---\n")
us <- didgpu(panel_df, "Y", "unit", "period", "D",
              effects = 3L, placebo = 0L, weight = "w",
              bootstrap_reps = 0L, backend = "r", verbose = FALSE)
print(as.numeric(us$results$Effects[, "Estimate"]))

cat("\nmax diff:",
    max(abs(as.numeric(ref_w$results$Effects[, 1]) -
            as.numeric(us$results$Effects[, "Estimate"]))), "\n")
