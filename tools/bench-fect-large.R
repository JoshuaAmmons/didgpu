# One-off: does the GPU fect SVD path actually win ABOVE the size gate?
# The gate is n_units >= 2000 AND n_units*n_periods >= 2e5. If the GPU
# loses even here, the CUDA fect path is effectively dead code and the
# gate threshold should be raised (or the path removed).
suppressPackageStartupMessages({ library(didgpu) })
if (!didgpu_has_cuda_support()) { cat("no CUDA\n"); quit(status = 0) }

make_panel <- function(n_units, n_periods, seed = 17L) {
  p <- didgpu_simulate_panel(n_units = n_units, n_periods = n_periods,
                              tau_profile = c(0.5, 1.0), seed = seed)
  p$D <- as.integer(p$D >= 0.5); p
}
time1 <- function(p, method, backend) {
  t0 <- Sys.time()
  invisible(didgpu_fect(p, "Y", "unit", "period", "D", method = method,
                         backend = backend, bootstrap_reps = 0L, verbose = FALSE))
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}
# Above-threshold sizes (n_units >= 2000, n_units*n_periods >= 2e5).
for (cfg in list(c(2000, 100), c(4000, 60), c(8000, 50))) {
  nu <- cfg[1]; np <- cfg[2]
  p <- make_panel(nu, np)
  for (method in c("ife", "mc")) {
    t_r <- time1(p, method, "r")
    t_g <- time1(p, method, "cuda")
    cat(sprintf("%s  n_units=%5d n_periods=%3d  R=%7.3fs  CUDA=%7.3fs  speedup=%5.2fx\n",
                method, nu, np, t_r, t_g, t_r / t_g))
  }
}
