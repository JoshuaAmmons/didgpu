library(didgpu)
p <- didgpu_simulate_panel_bidir(n_units = 80L, n_periods = 15L,
                                  frac_treated = 0.6, frac_in = 0.5,
                                  min_treat_period = 5L, max_treat_period = 10L,
                                  seed = 11L)
res <- didgpu_compare(p, "Y", "unit", "period", "D",
                       effects = 3L, placebo = 1L)
cat("\npass =", res$pass, "\n")
