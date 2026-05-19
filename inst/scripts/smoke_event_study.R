library(didgpu)
p <- didgpu_simulate_panel(n_units = 80L, n_periods = 18L,
                            frac_treated = 0.6,
                            min_treat_period = 7L, max_treat_period = 12L,
                            tau_profile = c(0.5, 1.0, 1.2, 1.0, 0.8),
                            sigma = 0.4, seed = 17L)
fit <- didgpu(p, "Y", "unit", "period", "D",
               effects = 4L, placebo = 3L,
               bootstrap_reps = 20L, seed = 1L,
               backend = "r", verbose = FALSE)
cat("--- event study data ---\n")
print(didgpu_event_study_data(fit))
