library(didgpu); library(DIDmultiplegtDYN)

p <- didgpu_simulate_panel(n_units = 80L, n_periods = 15L,
                            frac_treated = 0.6,
                            min_treat_period = 5L, max_treat_period = 10L,
                            tau_profile = c(0.5, 1.0, 1.2),
                            sigma = 0.4, seed = 17L)

cat("--- ref WITHOUT normalized ---\n")
ref_n <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(p), outcome = "Y", group = "unit",
  time = "period", treatment = "D",
  effects = 3, placebo = 0, graph_off = TRUE)))
print(as.numeric(ref_n$results$Effects[, 1]))

cat("\n--- ref WITH normalized = TRUE ---\n")
ref_y <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(p), outcome = "Y", group = "unit",
  time = "period", treatment = "D",
  effects = 3, placebo = 0, graph_off = TRUE, normalized = TRUE)))
print(as.numeric(ref_y$results$Effects[, 1]))

cat("\n--- diff ---\n")
print(as.numeric(ref_n$results$Effects[, 1]) - as.numeric(ref_y$results$Effects[, 1]))
