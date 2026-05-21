# One-off verification: does the GPU fect_mc fit agree with the R fit
# at a panel size ABOVE the size gate (where the GPU path actually
# runs)? The softthreshold kernel is unit-tested at small sizes, but
# the full iterated mc fit accumulates many kernel calls; this checks
# the end-to-end agreement on the regime users actually benefit from.
#
# Slow (R mc is ~minutes at this size) -> a tools/ script, not CI.
suppressPackageStartupMessages({ library(didgpu) })
if (!didgpu_has_cuda_support()) { cat("no CUDA\n"); quit(status = 0) }

make_panel <- function(n_units, n_periods, seed = 17L) {
  p <- didgpu_simulate_panel(n_units = n_units, n_periods = n_periods,
                              tau_profile = c(0.5, 1.0), seed = seed)
  p$D <- as.integer(p$D >= 0.5); p
}

# 2500 x 80 = 2e5 -> above the .fect_cuda_svd_worthwhile gate.
p <- make_panel(2500L, 80L)
cat("Panel: 2500 units x 80 periods (above mc GPU size gate)\n\n")

t0 <- Sys.time()
fit_r <- didgpu_fect(p, "Y", "unit", "period", "D", method = "mc",
                     backend = "r", bootstrap_reps = 0L, verbose = FALSE)
t_r <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

t0 <- Sys.time()
fit_c <- didgpu_fect(p, "Y", "unit", "period", "D", method = "mc",
                     backend = "cuda", bootstrap_reps = 0L, verbose = FALSE)
t_c <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

ate_r <- fit_r$results$ATE[1, "Estimate"]
ate_c <- fit_c$results$ATE[1, "Estimate"]
cat(sprintf("R    backend: ATE = %.10f   (%.1fs)\n", ate_r, t_r))
cat(sprintf("CUDA backend: ATE = %.10f   (%.1fs)\n", ate_c, t_c))
cat(sprintf("abs diff = %.3e   rel diff = %.3e\n",
            abs(ate_r - ate_c), abs(ate_r - ate_c) / abs(ate_r)))
cat(sprintf("speedup  = %.2fx\n", t_r / t_c))

# Event-study estimates too.
es_r <- fit_r$results$Effects[, "Estimate"]
es_c <- fit_c$results$Effects[, "Estimate"]
cat(sprintf("\nEvent-study max abs diff = %.3e\n",
            max(abs(es_r - es_c), na.rm = TRUE)))

if (abs(ate_r - ate_c) < 1e-4) {
  cat("\nVERDICT: GPU mc fit AGREES with R mc fit (< 1e-4). Path verified.\n")
} else {
  cat("\nVERDICT: DIVERGENCE > 1e-4 -- investigate before trusting GPU mc.\n")
}
