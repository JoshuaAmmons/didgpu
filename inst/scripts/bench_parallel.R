library(didgpu)

p <- didgpu_simulate_panel(
  n_units = 500L, n_periods = 100L, frac_treated = 0.6,
  tau_profile = c(0.5, 1.0, 1.2, 1.0, 0.8, 0.6, 0.5),
  sigma = 0.4, seed = 17L
)
n_cores <- parallel::detectCores()
cat(sprintf("panel: %d rows. n_cores available: %d\n", nrow(p), n_cores))

for (nw in c(1L, 2L, 4L, 8L)) {
  if (nw > n_cores) next
  cdir <- tempfile(sprintf("didgpu_par_%d_", nw))
  on.exit(unlink(cdir, recursive = TRUE), add = TRUE)
  t0 <- Sys.time()
  fit <- didgpu(p, "Y", "unit", "period", "D",
                 effects = 5L, placebo = 2L,
                 bootstrap_reps = 32L, seed = 1L,
                 checkpoint_dir = cdir,
                 backend = "r", verbose = FALSE,
                 n_workers = nw)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  n_workers=%d:  wall=%.2fs   effects[1]=%.4f  SE[1]=%.4f\n",
              nw, wall, fit$results$Effects[1, "Estimate"],
              fit$results$Effects[1, "SE"]))
}
