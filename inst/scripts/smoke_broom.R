library(didgpu)
p <- didgpu_simulate_panel(n_units = 60L, n_periods = 12L,
                            frac_treated = 0.6,
                            min_treat_period = 4L, max_treat_period = 9L,
                            tau_profile = c(0.5, 1.0, 1.2, 1.0),
                            sigma = 0.4, seed = 17L)
fit <- didgpu(p, "Y", "unit", "period", "D",
               effects = 3L, placebo = 1L,
               bootstrap_reps = 20L, seed = 1L,
               backend = "r", verbose = FALSE)
cat("--- didgpu_tidy(fit) ---\n")
print(didgpu_tidy(fit))
cat("\n--- didgpu_glance(fit) ---\n")
print(didgpu_glance(fit))

# Dispatch via tidy() / glance() also works (S3 methods registered).
cat("\n--- tidy.didgpu_result direct call ---\n")
print(tidy.didgpu_result(fit))
