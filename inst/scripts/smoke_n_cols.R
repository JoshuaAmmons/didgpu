library(didgpu); library(DIDmultiplegtDYN)
p <- didgpu_simulate_panel_bidir(n_units = 80L, n_periods = 15L,
                                  frac_treated = 0.6, frac_in = 0.5,
                                  min_treat_period = 5L, max_treat_period = 10L,
                                  seed = 11L)
ref <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(p), outcome = "Y", group = "unit", time = "period",
  treatment = "D", effects = 3, placebo = 1, graph_off = TRUE)))
us <- didgpu(p, "Y", "unit", "period", "D",
              effects = 3L, placebo = 1L,
              bootstrap_reps = 0L, backend = "r", verbose = FALSE)
cat("reference Effects:\n"); print(ref$results$Effects)
cat("\ndidgpu Effects:\n"); print(us$results$Effects)
cat("\nreference Placebos:\n"); print(ref$results$Placebos)
cat("\ndidgpu Placebos:\n"); print(us$results$Placebos)
