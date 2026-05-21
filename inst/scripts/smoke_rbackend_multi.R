# Compare r-backend vs reference at effects = 1, 2, 3, 5.
library(didgpu)
library(DIDmultiplegtDYN)

for (effects in c(1L, 2L, 3L, 5L)) {
  p <- didgpu_simulate_panel(
    n_units = 80L, n_periods = 18L, frac_treated = 0.6,
    min_treat_period = 4L, max_treat_period = 9L,
    tau_profile = c(0.5, 1.0, 1.2, 1.0, 0.8, 0.6),
    sigma = 0.4, seed = 17L
  )
  cat(sprintf("\n=== effects = %d ===\n", effects))
  ref <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
    df = as.data.frame(p), outcome = "Y", group = "unit",
    time = "period", treatment = "D",
    effects = as.double(effects), placebo = 0, graph_off = TRUE
  )))
  us <- didgpu(p, "Y", "unit", "period", "D",
                effects = effects, placebo = 0L,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  ref_e <- as.numeric(ref$results$Effects[, 1])
  us_e  <- as.numeric(us$results$Effects[, "Estimate"])
  cat("  reference: ", paste(sprintf("%.6f", ref_e), collapse = "  "), "\n")
  cat("  didgpu:    ", paste(sprintf("%.6f", us_e),  collapse = "  "), "\n")
  cat("  diff:      ", paste(sprintf("%.2e", abs(ref_e - us_e)), collapse = "  "), "\n")
}
