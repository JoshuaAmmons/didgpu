# Benchmark: IPW/DR cluster bootstrap, R vs CUDA. Now that IPW/DR run
# on the GPU (with per-cell influence functions), the cluster bootstrap
# gets the same IF-shortcut treatment as OR. This measures whether the
# ~200x win extends to the doubly-robust estimator.
suppressPackageStartupMessages({ library(didgpu) })
if (!didgpu_has_cuda_support()) { cat("no CUDA\n"); quit(status = 0) }

make_panel <- function(n_units, n_periods, seed = 17L) {
  p <- didgpu_simulate_panel(n_units = n_units, n_periods = n_periods,
                              tau_profile = c(0.5, 1.0), seed = seed)
  p$D <- as.integer(p$D >= 0.5)
  set.seed(seed + 1L)
  uv <- stats::rnorm(length(unique(p$unit)))
  names(uv) <- as.character(sort(unique(p$unit)))
  p$x1 <- uv[as.character(p$unit)]
  p
}
t1 <- function(p, method, backend, B) {
  t0 <- Sys.time()
  invisible(didgpu_cs(p, "Y", "unit", "period", "D", covariates = "x1",
                      est_method = method, bootstrap_reps = B,
                      bootstrap_kind = "cluster", backend = backend,
                      verbose = FALSE))
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

cat(sprintf("%-4s %-4s %-4s  %9s  %9s  %8s\n", "meth","nu","np","R","CUDA","speedup"))
for (method in c("IPW", "DR")) {
  for (nu in c(60L, 120L)) {
    p <- make_panel(nu, 10L)
    # warm up
    invisible(t1(p, method, "cuda", 0L))
    tr <- t1(p, method, "r",    200L)
    tc <- t1(p, method, "cuda", 200L)
    cat(sprintf("%-4s %-4d %-4d  %8.3fs  %8.3fs  %7.1fx\n",
                method, nu, 10L, tr, tc, tr / tc))
  }
}
