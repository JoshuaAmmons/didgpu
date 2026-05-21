library(didgpu); library(DIDmultiplegtDYN)

# Build the prepped panel exactly as the r-backend does, then call
# .core_one_event_time per direction. Compare to reference's
# switchers="in"/"out" per-direction estimates.
p <- didgpu_simulate_panel_bidir(n_units = 80L, n_periods = 15L,
                                  frac_treated = 0.6, frac_in = 0.5,
                                  min_treat_period = 5L, max_treat_period = 10L,
                                  seed = 11L)

# Reference per direction
cat("--- ref switchers=in ---\n")
ref_in <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(p), outcome = "Y", group = "unit", time = "period",
  treatment = "D", effects = 3, placebo = 0, graph_off = TRUE, switchers = "in")))
print(as.numeric(ref_in$results$Effects[, 1]))

cat("\n--- ref switchers=out ---\n")
ref_out <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(p), outcome = "Y", group = "unit", time = "period",
  treatment = "D", effects = 3, placebo = 0, graph_off = TRUE, switchers = "out")))
print(as.numeric(ref_out$results$Effects[, 1]))

# Now: run .core_one_event_time per direction directly, see what we get
prep    <- didgpu:::.prep_panel(p, "Y", "unit", "period", "D")
cat("\n--- didgpu .core_one_event_time, direction=1 (in only) ---\n")
for (k in 1:3) {
  res <- didgpu:::.core_one_event_time(prep, k = k, direction = 1L)
  cat(sprintf("  k=%d: att=%.6f  N_inc=%d\n", k, res$att, res$N_inc))
}
cat("\n--- didgpu .core_one_event_time, direction=0 (out only) ---\n")
for (k in 1:3) {
  res <- didgpu:::.core_one_event_time(prep, k = k, direction = 0L)
  # The reference flips this sign before publishing.
  cat(sprintf("  k=%d: att=%.6f  -att=%.6f  N_inc=%d\n",
              k, res$att, -res$att, res$N_inc))
}
