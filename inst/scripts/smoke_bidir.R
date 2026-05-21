# Test r-backend on a panel with BOTH switcher directions.
library(didgpu)
library(DIDmultiplegtDYN)

p <- didgpu_simulate_panel_bidir(
  n_units = 80L, n_periods = 15L, frac_treated = 0.6, frac_in = 0.5,
  min_treat_period = 5L, max_treat_period = 10L, seed = 11L
)
cat("D distribution:\n"); print(table(p$D))
cat("\nFirst-period treatment d_sq:\n")
print(table(p$D[p$period == 1L]))
cat("\nDirections of treated units (from truth):\n")
print(table(attr(p, "truth")$direction, useNA = "ifany"))

ref <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(p), outcome = "Y", group = "unit",
  time = "period", treatment = "D",
  effects = 3, placebo = 1, graph_off = TRUE
)))
us <- didgpu(p, "Y", "unit", "period", "D",
              effects = 3L, placebo = 1L,
              bootstrap_reps = 0L, backend = "r", verbose = FALSE)

cat("\nref effects:    ", sprintf("%.6f", as.numeric(ref$results$Effects[, 1])), "\n")
cat("didgpu effects: ", sprintf("%.6f", as.numeric(us$results$Effects[, "Estimate"])), "\n")
cat("max diff:       ", sprintf("%.2e",
      max(abs(as.numeric(ref$results$Effects[, 1]) -
              as.numeric(us$results$Effects[, "Estimate"])))), "\n")

cat("\nref placebos:    ", sprintf("%.6f", as.numeric(ref$results$Placebos[, 1])), "\n")
cat("didgpu placebos: ", sprintf("%.6f", as.numeric(us$results$Placebos[, "Estimate"])), "\n")

cat("\nref ATE:    ", as.numeric(ref$results$ATE[1, 1]), "\n")
cat("didgpu ATE: ", as.numeric(us$results$ATE[1, "Estimate"]), "\n")
