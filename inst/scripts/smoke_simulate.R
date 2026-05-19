# Quick smoke test for the simulator. Run via:
#   Rscript inst/scripts/smoke_simulate.R
suppressPackageStartupMessages({
  library(data.table)
})
source("R/simulate.R")
p <- didgpu_simulate_panel(n_units = 40L, n_periods = 12L, frac_treated = 0.6,
                            min_treat_period = 4L, max_treat_period = 9L,
                            tau_profile = c(0.5, 1.0, 1.2, 1.1),
                            sigma = 0.3, seed = 42L)
cat("rows:", nrow(p), "  cols:", ncol(p), "\n")
cat("head:\n"); print(head(p))
cat("D distribution:\n"); print(table(p$D))
cat("First-treatment periods F_g (table):\n")
print(table(attr(p, "truth")$F_g))
cat("Mean Y by D:\n")
print(aggregate(Y ~ D, p, mean))
