# Explore #88: where does didgpu_loo (CS family) spend its time, and
# does backend="cuda" already help? LOO refits are POINT estimates
# (bootstrap_reps=0), which are ~1x on GPU (copy-bound). This measures
# whether a batched LOO kernel is worth building or whether LOO is
# dominated by R-side orchestration (cell construction) -- in which
# case a kernel won't move the needle (cf. the fect #86/#87 finding).
suppressPackageStartupMessages({ library(didgpu) })
cat("CUDA: ", didgpu_has_cuda_support(), "\n", sep = "")

make_panel <- function(n_units, n_periods, seed = 17L) {
  p <- didgpu_simulate_panel(n_units = n_units, n_periods = n_periods,
                              tau_profile = c(0.5, 1.0), seed = seed)
  p$D <- as.integer(p$D >= 0.5); p
}

time_loo <- function(p, fit, by, backend) {
  fit$args$backend <- backend
  t0 <- Sys.time()
  invisible(didgpu_loo(fit, by = by, df = p, verbose = FALSE))
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

for (nu in c(80L, 200L)) {
  p <- make_panel(nu, 12L)
  fit_r <- didgpu_cs(p, "Y", "unit", "period", "D", est_method = "OR",
                      aggregation = "overall", bootstrap_reps = 0L,
                      backend = "r", verbose = FALSE)
  # how many cohorts / units?
  d <- data.table::as.data.table(p)
  ncoh <- length(unique(d[D == 1L, .(first = min(period)), by = unit]$first))
  cat(sprintf("\n--- n_units=%d, n_periods=12, ~%d cohorts ---\n", nu, ncoh))
  for (by in c("cohort", "unit")) {
    tr <- time_loo(p, fit_r, by, "r")
    tc <- if (didgpu_has_cuda_support()) time_loo(p, fit_r, by, "cuda") else NA
    cat(sprintf("  by=%-7s  R=%7.3fs  CUDA=%7.3fs  speedup=%.2fx\n",
                by, tr, tc, tr / tc))
  }
  # Single point-estimate refit cost (to see the per-refit baseline).
  t0 <- Sys.time()
  invisible(didgpu_cs(p, "Y", "unit", "period", "D", est_method = "OR",
                      aggregation = "overall", bootstrap_reps = 0L,
                      backend = "r", verbose = FALSE))
  cat(sprintf("  [one CS point-estimate refit: %.3fs]\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
