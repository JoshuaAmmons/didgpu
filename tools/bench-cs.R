# ============================================================================
# Microbenchmark: didgpu_cs(backend = "r") vs backend = "cuda".
#
# Measures end-to-end wall-clock time per call across a small grid of
# panel sizes. Runs each configuration `reps` times and reports
# median, min, max so transient noise (JIT, OS jitter) is visible.
#
# Usage (from WSL):
#   Rscript tools/bench-cs.R
# ============================================================================

suppressPackageStartupMessages({
  library(didgpu)
})

cat("=== didgpu_cs CUDA vs R benchmark ===\n")
cat("CUDA compiled in: ", didgpu_has_cuda_support(), "\n\n", sep = "")

if (!didgpu_has_cuda_support()) {
  cat("CUDA not available; nothing to compare.\n"); quit(status = 0)
}

# Time one didgpu_cs call (median of `reps` runs).
time_call <- function(p, backend, reps = 3L,
                       bootstrap_reps = 0L,
                       bootstrap_kind = "cluster") {
  # Warm up.
  invisible(didgpu_cs(p, "Y", "unit", "period", "D",
                       est_method = "OR", aggregation = "event",
                       bootstrap_reps = 0L,
                       backend = backend, verbose = FALSE))
  ts <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- Sys.time()
    invisible(didgpu_cs(p, "Y", "unit", "period", "D",
                         est_method = "OR", aggregation = "event",
                         bootstrap_reps = bootstrap_reps,
                         bootstrap_kind = bootstrap_kind,
                         backend = backend, verbose = FALSE))
    ts[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  ts
}

# Build a CS-friendly panel: staggered adoption, integer D.
make_panel <- function(n_units, n_periods, seed = 17L) {
  p <- didgpu_simulate_panel(n_units = n_units, n_periods = n_periods,
                              tau_profile = c(0.5, 1.0), seed = seed)
  p$D <- as.integer(p$D >= 0.5)
  p
}

rows <- list()
for (n_units in c(40L, 100L, 200L)) {
  for (n_periods in c(8L, 12L)) {
    p <- make_panel(n_units, n_periods)

    # Point-estimate only.
    t_r <- time_call(p, "r",    reps = 3L, bootstrap_reps = 0L)
    t_g <- time_call(p, "cuda", reps = 3L, bootstrap_reps = 0L)
    rows[[length(rows) + 1L]] <- data.frame(
      scenario  = "point_estimate",
      n_units   = n_units, n_periods = n_periods,
      r_median  = median(t_r), cuda_median = median(t_g),
      speedup   = median(t_r) / median(t_g))

    # Cluster bootstrap (uses .cs_cluster_bootstrap_cuda fast path).
    t_r_cb <- time_call(p, "r",    reps = 3L,
                         bootstrap_reps = 200L, bootstrap_kind = "cluster")
    t_g_cb <- time_call(p, "cuda", reps = 3L,
                         bootstrap_reps = 200L, bootstrap_kind = "cluster")
    rows[[length(rows) + 1L]] <- data.frame(
      scenario  = "cluster_bootstrap_B200",
      n_units   = n_units, n_periods = n_periods,
      r_median  = median(t_r_cb), cuda_median = median(t_g_cb),
      speedup   = median(t_r_cb) / median(t_g_cb))

    # Multiplier bootstrap.
    t_r_mb <- time_call(p, "r",    reps = 3L,
                         bootstrap_reps = 200L, bootstrap_kind = "multiplier")
    t_g_mb <- time_call(p, "cuda", reps = 3L,
                         bootstrap_reps = 200L, bootstrap_kind = "multiplier")
    rows[[length(rows) + 1L]] <- data.frame(
      scenario  = "multiplier_bootstrap_B200",
      n_units   = n_units, n_periods = n_periods,
      r_median  = median(t_r_mb), cuda_median = median(t_g_mb),
      speedup   = median(t_r_mb) / median(t_g_mb))
  }
}

out <- do.call(rbind, rows)
out$r_median    <- sprintf("%6.3fs", out$r_median)
out$cuda_median <- sprintf("%6.3fs", out$cuda_median)
out$speedup     <- sprintf("%5.2fx", out$speedup)
print(out, row.names = FALSE)
