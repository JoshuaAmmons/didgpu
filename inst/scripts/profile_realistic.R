library(didgpu)
p <- didgpu_simulate_panel(
  n_units = 1000L, n_periods = 280L,
  frac_treated = 0.6,
  min_treat_period = 50L, max_treat_period = 200L,
  tau_profile = c(seq(0.1, 1.0, length.out = 10), rep(1.0, 100)),
  sigma = 0.4, seed = 17L
)
# warmup
invisible(didgpu(p, "Y", "unit", "period", "D",
                  effects = 30L, placebo = 10L,
                  bootstrap_reps = 0L, backend = "r", verbose = FALSE))

tmp <- tempfile("rprof_", fileext = ".out")
Rprof(tmp, interval = 0.02)
fit <- didgpu(p, "Y", "unit", "period", "D",
               effects = 90L, placebo = 30L,
               bootstrap_reps = 0L, backend = "r", verbose = FALSE)
Rprof(NULL)
cat("Total fit time recorded.\n")
s <- summaryRprof(tmp)
cat("--- top by total time ---\n")
print(head(s$by.total, 20))
cat("\n--- top by self time ---\n")
print(head(s$by.self, 20))
cat("\nsampling time total:", s$sampling.time, "\n")
