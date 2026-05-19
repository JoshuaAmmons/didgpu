library(didgpu)
p <- didgpu_simulate_panel(n_units = 200L, n_periods = 50L,
                            frac_treated = 0.6,
                            tau_profile = c(0.5, 1.0, 1.2, 1.0, 0.8),
                            sigma = 0.4, seed = 17L)

cat("--- estimate sequential 100-rep bootstrap ---\n")
e1 <- didgpu_estimate_runtime(p, "Y", "unit", "period", "D",
                                effects = 5L, placebo = 2L,
                                bootstrap_reps = 100L,
                                n_workers = 1L)
str(e1, max.level = 1)

cat("\n--- estimate with 8 workers ---\n")
e2 <- didgpu_estimate_runtime(p, "Y", "unit", "period", "D",
                                effects = 5L, placebo = 2L,
                                bootstrap_reps = 100L,
                                n_workers = 8L)
str(e2, max.level = 1)
