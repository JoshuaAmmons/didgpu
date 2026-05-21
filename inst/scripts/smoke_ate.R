# Compare r-backend ATE vs reference ATE for multiple effects.
library(didgpu)
library(DIDmultiplegtDYN)

p <- didgpu_simulate_panel(
  n_units = 80L, n_periods = 18L, frac_treated = 0.6,
  min_treat_period = 4L, max_treat_period = 9L,
  tau_profile = c(0.5, 1.0, 1.2, 1.0, 0.8),
  sigma = 0.4, seed = 17L
)

for (eff in c(1L, 2L, 3L, 5L)) {
  cat(sprintf("\n=== effects = %d ===\n", eff))
  ref <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
    df = as.data.frame(p), outcome = "Y", group = "unit",
    time = "period", treatment = "D",
    effects = as.double(eff), placebo = 0, graph_off = TRUE
  )))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = eff, placebo = 0L,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  ref_ate <- as.numeric(ref$results$ATE[1, 1])
  us_ate  <- as.numeric(us$results$ATE[1, "Estimate"])
  cat("  reference ATE: ", sprintf("%.6f", ref_ate), "\n")
  cat("  didgpu ATE:    ", sprintf("%.6f", us_ate), "\n")
  cat("  diff:          ", sprintf("%.2e", abs(ref_ate - us_ate)), "\n")
}
