# Compare the partial r-backend against the reference for the simplest
# case it supports (effects=1, placebo=0, no controls).
library(didgpu)
library(DIDmultiplegtDYN)

p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L, frac_treated = 0.6,
                            min_treat_period = 4L, max_treat_period = 9L,
                            tau_profile = c(0.5, 1.0), sigma = 0.4, seed = 17L)

cat("--- reference (effects=1, placebo=0) ---\n")
ref <- did_multiplegt_dyn(df = as.data.frame(p),
                           outcome = "Y", group = "unit",
                           time = "period", treatment = "D",
                           effects = 1, placebo = 0, graph_off = TRUE)
cat("  ref effect_1: ", ref$results$Effects[1, 1], "\n")

cat("--- didgpu, r backend ---\n")
us <- didgpu(p, "Y", "unit", "period", "D",
              effects = 1L, placebo = 0L,
              bootstrap_reps = 0L, backend = "r", verbose = FALSE)
cat("  didgpu effect_1: ", us$results$Effects[1, "Estimate"], "\n")
cat("  diff:            ", us$results$Effects[1, "Estimate"] - ref$results$Effects[1, 1], "\n")
