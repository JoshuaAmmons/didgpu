# Probe: does my existing code handle multivalued treatment?
library(didgpu); library(DIDmultiplegtDYN); library(data.table)

set.seed(11)
n_units <- 80L; n_periods <- 15L
# Three treatment levels: 0, 1, 2. Each unit's baseline is one of these.
baseline <- sample(0:2, n_units, replace = TRUE)
# Switchers go to a different level (uniform among the other 2).
F_g <- rep(Inf, n_units)
treated <- sort(sample(seq_len(n_units), 48L))
F_g[treated] <- sample(5L:10L, 48L, replace = TRUE)
post_level <- baseline
for (u in treated) {
  others <- setdiff(0:2, baseline[u])
  post_level[u] <- sample(others, 1L)
}
unit_fe <- rnorm(n_units, 0, 1); time_fe <- rnorm(n_periods, 0, 0.3)

panel <- data.table(unit = rep(1:n_units, each = n_periods),
                     period = rep(1:n_periods, n_units))
panel[, baseline := baseline[unit]]
panel[, F_g := F_g[unit]]
panel[, post_level := post_level[unit]]
panel[, D := ifelse(period >= F_g, post_level, baseline)]
panel[, k_evt := period - F_g]
panel[, tau_k := 0]
# Treatment effect = 0.5 per unit of treatment change at k=0, etc.
panel[is.finite(F_g) & k_evt >= 0,
      tau_k := (post_level - baseline) * c(0.3, 0.5, 0.6, 0.5)[pmin(k_evt + 1L, 4L)]]
panel[, Y := unit_fe[unit] + time_fe[period] + tau_k + rnorm(.N, 0, 0.4)]
panel_df <- as.data.frame(panel[order(unit, period), .(unit, period, D, Y)])

cat("D distribution:\n"); print(table(panel_df$D))
cat("\nFirst-period D (d_sq):\n"); print(table(panel_df$D[panel_df$period == 1L]))

cat("\n--- reference (multivalued, no special flags) ---\n")
ref <- tryCatch(
  did_multiplegt_dyn(df = panel_df, outcome = "Y", group = "unit",
                     time = "period", treatment = "D",
                     effects = 3, placebo = 1, graph_off = TRUE),
  error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(ref)) {
  cat("ref effects: ", round(as.numeric(ref$results$Effects[, 1]), 4), "\n")
}

cat("\n--- didgpu (r backend, multivalued) ---\n")
us <- tryCatch(
  didgpu(panel_df, "Y", "unit", "period", "D",
          effects = 3L, placebo = 1L,
          bootstrap_reps = 0L, backend = "r", verbose = FALSE),
  error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(us)) {
  cat("didgpu effects: ", round(as.numeric(us$results$Effects[, "Estimate"]), 4), "\n")
}
if (!is.null(ref) && !is.null(us)) {
  cat("\nmax diff effects:",
      max(abs(as.numeric(ref$results$Effects[, 1]) -
              as.numeric(us$results$Effects[, "Estimate"]))), "\n")
}
