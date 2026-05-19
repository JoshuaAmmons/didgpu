library(didgpu)
p <- didgpu_simulate_panel(n_units = 500L, n_periods = 100L,
                            frac_treated = 0.6,
                            tau_profile = c(0.5, 1.0, 1.2, 1.0, 0.8, 0.6,
                                            0.5, 0.4, 0.3, 0.2),
                            sigma = 0.4, seed = 17L)
invisible(didgpu(p, "Y", "unit", "period", "D",
                  effects = 5L, placebo = 2L,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE))

tmp <- tempfile("rprof_", fileext = ".out")
Rprof(tmp, interval = 0.005)
for (rep in 1:8) {
  didgpu(p, "Y", "unit", "period", "D",
          effects = 5L, placebo = 2L,
          bootstrap_reps = 0L, backend = "r", verbose = FALSE)
}
Rprof(NULL)
s <- summaryRprof(tmp)
cat("--- by self ---\n")
print(head(s$by.self, 30))
cat("\n--- by total ---\n")
print(head(s$by.total, 20))
cat("\nsampling time total:", s$sampling.time, "\n")
