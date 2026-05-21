# Benchmark at the project's actual config: effects=90, placebo=90,
# ~280K row panel (close to the DAO daily-frequency Stage 1 dimensions).
# Single point-estimate fit only (no bootstrap) to bound the per-iter
# cost, since a 100-rep bootstrap is just 100x this.
library(didgpu)

# 1000 units × 280 periods = 280,000 rows.
p <- didgpu_simulate_panel(
  n_units = 1000L, n_periods = 280L,
  frac_treated = 0.6,
  min_treat_period = 50L, max_treat_period = 200L,
  tau_profile = c(seq(0.1, 1.0, length.out = 10), rep(1.0, 100)),
  sigma = 0.4, seed = 17L
)
cat(sprintf("panel: %d rows, %d unique units, %d unique periods\n",
            nrow(p), length(unique(p$unit)), length(unique(p$period))))

# Effects = 90, placebos = 30 (won't go to 90 because of horizon clamping
# given our treatment window, but the panel will accept whatever fits).
for (cfg in list(c(eff = 30L, pl = 10L),
                  c(eff = 90L, pl = 30L))) {
  cat(sprintf("\n=== effects = %d, placebos = %d ===\n", cfg["eff"], cfg["pl"]))
  # Warmup
  invisible(didgpu(p, "Y", "unit", "period", "D",
                    effects = as.integer(cfg["eff"]),
                    placebo = as.integer(cfg["pl"]),
                    bootstrap_reps = 0L, backend = "r", verbose = FALSE))
  t0 <- Sys.time()
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = as.integer(cfg["eff"]),
                 placebo = as.integer(cfg["pl"]),
                 bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  r-backend wall: %.2fs (effective effects=%d, placebos=%d)\n",
              wall, fit$results$N_Effects, fit$results$N_Placebos))
  cat(sprintf("  Implied 100-rep bootstrap: %.1f min\n", wall * 100 / 60))
}
