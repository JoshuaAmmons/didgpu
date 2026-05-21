library(didgpu); library(DIDmultiplegtDYN); library(data.table)

p_full <- as.data.table(didgpu_simulate_panel_bidir(
  n_units = 80L, n_periods = 15L, frac_treated = 0.6, frac_in = 0.5,
  min_treat_period = 5L, max_treat_period = 10L, seed = 11L
))
truth <- attr(didgpu_simulate_panel_bidir(
  n_units = 80L, n_periods = 15L, frac_treated = 0.6, frac_in = 0.5,
  min_treat_period = 5L, max_treat_period = 10L, seed = 11L), "truth")
p_full[, direction := truth$direction[as.character(unit)]]

# Subset to ONLY out-switchers (drop in-switchers and never-treated entirely)
p_out_only <- p_full[!is.na(direction) & direction == 0L, list(unit, period, D, Y)]
# Re-number unit IDs to be consecutive
p_out_only[, unit := as.integer(factor(unit))]
cat("p_out_only: ", nrow(p_out_only), "rows,", length(unique(p_out_only$unit)), "units\n")
cat("D at t=1:\n"); print(table(p_out_only$D[p_out_only$period == 1L]))

cat("\n--- ref on out-only subset ---\n")
ref <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(p_out_only), outcome = "Y", group = "unit",
  time = "period", treatment = "D",
  effects = 3, placebo = 0, graph_off = TRUE)))
print(as.numeric(ref$results$Effects[, 1]))

cat("\n--- didgpu on out-only subset ---\n")
us <- didgpu(as.data.frame(p_out_only), "Y", "unit", "period", "D",
              effects = 3L, placebo = 0L, bootstrap_reps = 0L,
              backend = "r", verbose = FALSE)
print(as.numeric(us$results$Effects[, "Estimate"]))

# But — does the reference even run on a panel with no never-changers in d_sq=1 cohort?
# (Since the only d_sq=1 units here all switch eventually.)
cat("\n--- reference behavior summary ---\n")
cat("If the reference complained, that explains the bidir difference.\n")
