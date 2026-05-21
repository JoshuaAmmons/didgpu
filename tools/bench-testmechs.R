# ============================================================================
# Microbenchmark: TestMechs nonparametric bootstrap, R vs CUDA.
#
# Times the bootstrap primitive directly (.testmechs_bootstrap_r vs
# .testmechs_bootstrap_cuda) across a grid of (n, B). The bootstrap is
# the dominant cost of the analytic-variance-off path in
# didgpu_test_sharp_null().
#
# Usage:
#   Rscript tools/bench-testmechs.R
# ============================================================================

suppressPackageStartupMessages({ library(didgpu) })

cat("=== TestMechs bootstrap CUDA vs R benchmark ===\n")
cat("CUDA compiled in: ", didgpu_has_cuda_support(), "\n\n", sep = "")
if (!didgpu_has_cuda_support()) {
  cat("CUDA not available; nothing to compare.\n"); quit(status = 0)
}

make_dmy <- function(n, K = 2L, dy = 3L, seed = 17L) {
  set.seed(seed)
  list(d = sample.int(2L, n, replace = TRUE) - 1L,
       m = sample.int(K,  n, replace = TRUE),
       y = sample.int(dy, n, replace = TRUE))
}

time_fn <- function(fn, reps = 3L) {
  fn()  # warm up
  ts <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- Sys.time()
    fn()
    ts[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  median(ts)
}

rows <- list()
for (n in c(500L, 2000L, 10000L)) {
  for (B in c(200L, 1000L)) {
    dat <- make_dmy(n)
    t_r <- time_fn(function()
      didgpu:::.testmechs_bootstrap_r(dat$d, dat$m, dat$y, B,
                                       method = "nonparametric", seed = 17L))
    t_g <- time_fn(function()
      didgpu:::.testmechs_bootstrap_cuda(dat$d, dat$m, dat$y, B,
                                          method = "nonparametric", seed = 17L))
    rows[[length(rows) + 1L]] <- data.frame(
      n = n, B = B,
      r_median = t_r, cuda_median = t_g,
      speedup = t_r / t_g)
  }
}

out <- do.call(rbind, rows)
out$r_median    <- sprintf("%7.4fs", out$r_median)
out$cuda_median <- sprintf("%7.4fs", out$cuda_median)
out$speedup     <- sprintf("%6.2fx", out$speedup)
print(out, row.names = FALSE)
