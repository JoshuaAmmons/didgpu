# Smoke-test: pipe the simulated panel through DIDmultiplegtDYN to (a) confirm
# the reference can fit it cleanly, (b) inspect the return object so we know
# what shape to target for output compatibility.
suppressPackageStartupMessages({
  library(data.table)
  library(DIDmultiplegtDYN)
})
source("R/simulate.R")

p <- didgpu_simulate_panel(n_units = 80L, n_periods = 16L, frac_treated = 0.6,
                            min_treat_period = 4L, max_treat_period = 12L,
                            tau_profile = c(0.5, 1.0, 1.2, 1.0, 0.8),
                            sigma = 0.4, seed = 17L)

cat("=== fitting reference (effects=5, placebo=3, no bootstrap) ===\n")
fit <- did_multiplegt_dyn(
  df = as.data.frame(p),
  outcome = "Y", group = "unit", time = "period", treatment = "D",
  effects = 5L, placebo = 3L, graph_off = TRUE
)
cat("\n--- str(fit) ---\n")
str(fit, max.level = 2)
cat("\n--- top-level names ---\n")
print(names(fit))
cat("\n--- class ---\n")
print(class(fit))

# Save the return object so we can inspect shape later.
saveRDS(fit, "inst/scripts/reference_fit_smoke.rds")
cat("\nsaved -> inst/scripts/reference_fit_smoke.rds\n")
