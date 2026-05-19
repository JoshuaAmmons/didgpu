# Scaling benchmark: r-backend vs reference on increasing panel size.
# Measures wall time PER POINT-ESTIMATE FIT (no bootstrap; the bootstrap
# loop is the same outer code for both backends).
library(didgpu)
library(DIDmultiplegtDYN)

bench_one <- function(n_units, n_periods, effects, placebo, seed = 1L,
                       runs = 3L) {
  p <- didgpu_simulate_panel(
    n_units = n_units, n_periods = n_periods, frac_treated = 0.6,
    tau_profile = c(0.5, 1.0, 1.2, 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2),
    sigma = 0.4, seed = seed
  )
  cat(sprintf("\n--- %d units x %d periods (%d rows), effects=%d, placebo=%d ---\n",
              n_units, n_periods, nrow(p), effects, placebo))

  # Warmup (R caches some compile + datatable infra on first call)
  invisible(didgpu(p, "Y", "unit", "period", "D",
                    effects = effects, placebo = placebo,
                    bootstrap_reps = 0L, backend = "r", verbose = FALSE))

  t_r <- replicate(runs, {
    t0 <- Sys.time()
    invisible(didgpu(p, "Y", "unit", "period", "D",
                      effects = effects, placebo = placebo,
                      bootstrap_reps = 0L, backend = "r", verbose = FALSE))
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  })
  t_ref <- replicate(runs, {
    t0 <- Sys.time()
    invisible(didgpu(p, "Y", "unit", "period", "D",
                      effects = effects, placebo = placebo,
                      bootstrap_reps = 0L, backend = "reference", verbose = FALSE))
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  })

  cat(sprintf("  r-backend:  median %.3fs  (runs: %s)\n",
              median(t_r),  paste(sprintf("%.3f", t_r),  collapse = ", ")))
  cat(sprintf("  reference:  median %.3fs  (runs: %s)\n",
              median(t_ref), paste(sprintf("%.3f", t_ref), collapse = ", ")))
  cat(sprintf("  speedup:    %.1fx\n", median(t_ref) / median(t_r)))
  invisible(list(units = n_units, periods = n_periods,
                  r = median(t_r), ref = median(t_ref)))
}

# Progressively larger panels. Keep effects + placebo modest so the
# comparison is apples-to-apples (and so the reference doesn't time out).
results <- list()
results[[1]] <- bench_one(100L,    20L,  3L, 1L)
results[[2]] <- bench_one(200L,    50L,  3L, 1L)
results[[3]] <- bench_one(500L,   100L,  5L, 2L)
results[[4]] <- bench_one(1000L,  200L,  5L, 2L)

cat("\n\n=== summary ===\n")
cat(sprintf("%-20s  %-8s  %-12s  %-12s  %-8s\n",
            "panel", "rows", "r (s)", "ref (s)", "speedup"))
for (r in results) {
  cat(sprintf("%-20s  %-8d  %-12.3f  %-12.3f  %-8.1f\n",
              sprintf("%dx%d", r$units, r$periods),
              r$units * r$periods, r$r, r$ref, r$ref / r$r))
}
