library(didgpu); library(DIDmultiplegtDYN)
p <- didgpu_simulate_panel_bidir(n_units = 80L, n_periods = 15L,
                                  frac_treated = 0.6, frac_in = 0.5,
                                  min_treat_period = 5L, max_treat_period = 10L,
                                  seed = 11L)
cat("--- ref switchers=in ---\n")
ref_in <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(p), outcome = "Y", group = "unit", time = "period",
  treatment = "D", effects = 3, placebo = 0, graph_off = TRUE, switchers = "in")))
print(ref_in$results$Effects[, 1])

cat("\n--- ref switchers=out ---\n")
ref_out <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(p), outcome = "Y", group = "unit", time = "period",
  treatment = "D", effects = 3, placebo = 0, graph_off = TRUE, switchers = "out")))
print(ref_out$results$Effects[, 1])

cat("\n--- ref switchers= (both) ---\n")
ref_both <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(p), outcome = "Y", group = "unit", time = "period",
  treatment = "D", effects = 3, placebo = 0, graph_off = TRUE)))
print(ref_both$results$Effects[, 1])
