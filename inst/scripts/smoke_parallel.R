# Test parallel bootstrap workers.
library(didgpu)

p <- didgpu_simulate_panel(n_units = 60L, n_periods = 14L,
                            frac_treated = 0.6,
                            min_treat_period = 4L, max_treat_period = 9L,
                            tau_profile = c(0.5, 1.0, 1.2),
                            sigma = 0.4, seed = 17L)

cat("=== sequential, 20 reps ===\n")
t0 <- Sys.time()
fit_seq <- didgpu(p, "Y", "unit", "period", "D",
                   effects = 3L, placebo = 1L,
                   bootstrap_reps = 20L, seed = 1L,
                   backend = "r", verbose = FALSE, n_workers = 1L)
cat(sprintf("  wall: %.2fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

cat("\n=== parallel (4 workers), 20 reps ===\n")
t0 <- Sys.time()
fit_par <- didgpu(p, "Y", "unit", "period", "D",
                   effects = 3L, placebo = 1L,
                   bootstrap_reps = 20L, seed = 1L,
                   backend = "r", verbose = FALSE, n_workers = 4L)
cat(sprintf("  wall: %.2fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

cat("\n--- point estimates (should be identical) ---\n")
cat("  seq: ", fit_seq$results$Effects[, "Estimate"], "\n")
cat("  par: ", fit_par$results$Effects[, "Estimate"], "\n")
cat("  diff:", max(abs(fit_seq$results$Effects[, "Estimate"] -
                        fit_par$results$Effects[, "Estimate"])), "\n")

cat("\n--- SEs ---\n")
cat("  seq: ", fit_seq$results$Effects[, "SE"], "\n")
cat("  par: ", fit_par$results$Effects[, "SE"], "\n")
cat("  diff:", max(abs(fit_seq$results$Effects[, "SE"] -
                        fit_par$results$Effects[, "SE"])), "\n")
