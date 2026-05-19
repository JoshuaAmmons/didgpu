# ============================================================================
# Microbenchmark: didgpu_fect(method = "ife"|"mc") CUDA vs R.
#
# Measures the cost of one fect_ife / fect_mc fit (point estimate
# only — bootstrap is a separate question) across a few panel sizes.
# fect_ife alternates FE step + truncated-SVD step until convergence.
# The per-iter SVD is currently on the GPU (Phase 1 #80) but the
# alternation loop is on the host; Phase 2 #86 would fuse them.
#
# Usage:
#   Rscript tools/bench-fect.R
# ============================================================================

suppressPackageStartupMessages({ library(didgpu) })

cat("=== didgpu_fect CUDA vs R benchmark ===\n")
cat("CUDA compiled in: ", didgpu_has_cuda_support(), "\n\n", sep = "")
if (!didgpu_has_cuda_support()) {
  cat("CUDA not available; nothing to compare.\n"); quit(status = 0)
}

time_call <- function(p, method, backend, reps = 3L) {
  invisible(didgpu_fect(p, "Y", "unit", "period", "D",
                         method = method, backend = backend,
                         bootstrap_reps = 0L, verbose = FALSE))
  ts <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- Sys.time()
    invisible(didgpu_fect(p, "Y", "unit", "period", "D",
                           method = method, backend = backend,
                           bootstrap_reps = 0L, verbose = FALSE))
    ts[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  ts
}

make_panel <- function(n_units, n_periods, seed = 17L) {
  p <- didgpu_simulate_panel(n_units = n_units, n_periods = n_periods,
                              tau_profile = c(0.5, 1.0), seed = seed)
  p$D <- as.integer(p$D >= 0.5)
  p
}

rows <- list()
for (method in c("ife", "mc")) {
  for (n_units in c(40L, 100L, 200L)) {
    for (n_periods in c(10L, 20L)) {
      p <- make_panel(n_units, n_periods)
      t_r <- time_call(p, method, "r",    reps = 3L)
      t_g <- time_call(p, method, "cuda", reps = 3L)
      rows[[length(rows) + 1L]] <- data.frame(
        method    = method,
        n_units   = n_units, n_periods = n_periods,
        r_median  = median(t_r), cuda_median = median(t_g),
        speedup   = median(t_r) / median(t_g))
    }
  }
}

out <- do.call(rbind, rows)
out$r_median    <- sprintf("%7.3fs", out$r_median)
out$cuda_median <- sprintf("%7.3fs", out$cuda_median)
out$speedup     <- sprintf("%6.2fx", out$speedup)
print(out, row.names = FALSE)
