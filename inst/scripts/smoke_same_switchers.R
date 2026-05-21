library(didgpu); library(DIDmultiplegtDYN)

p <- didgpu_simulate_panel(n_units = 80L, n_periods = 15L,
                            frac_treated = 0.6,
                            min_treat_period = 5L, max_treat_period = 10L,
                            tau_profile = c(0.5, 1.0, 1.2),
                            sigma = 0.4, seed = 17L)

cat("--- ref WITH same_switchers = TRUE ---\n")
ref <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(p), outcome = "Y", group = "unit",
  time = "period", treatment = "D",
  effects = 3, placebo = 0, graph_off = TRUE,
  same_switchers = TRUE)))
print(as.numeric(ref$results$Effects[, 1]))

cat("\n--- didgpu WITH same_switchers = TRUE ---\n")
us <- didgpu(p, "Y", "unit", "period", "D",
              effects = 3L, placebo = 0L, same_switchers = TRUE,
              bootstrap_reps = 0L, backend = "r", verbose = FALSE)
print(as.numeric(us$results$Effects[, "Estimate"]))

cat("\nmax diff:",
    max(abs(as.numeric(ref$results$Effects[, 1]) -
            as.numeric(us$results$Effects[, "Estimate"]))), "\n")
