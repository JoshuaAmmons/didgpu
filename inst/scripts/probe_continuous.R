# Probe: does my code handle continuous treatment?
library(didgpu); library(DIDmultiplegtDYN); library(data.table)

set.seed(11)
n_units <- 80L; n_periods <- 15L
F_g <- rep(Inf, n_units)
treated <- sort(sample(seq_len(n_units), 48L))
F_g[treated] <- sample(5L:10L, 48L, replace = TRUE)
# Continuous baseline treatment intensity.
baseline_d <- runif(n_units, 0, 1)
# Switchers move to a new continuous intensity.
post_d <- baseline_d + rnorm(n_units, 0, 0.5)
unit_fe <- rnorm(n_units, 0, 1); time_fe <- rnorm(n_periods, 0, 0.3)

panel <- data.table(unit = rep(1:n_units, each = n_periods),
                     period = rep(1:n_periods, n_units))
panel[, F_g := F_g[unit]]
panel[, D := ifelse(period >= F_g, post_d[unit], baseline_d[unit])]
panel[, k_evt := period - F_g]
panel[, tau_k := 0]
panel[is.finite(F_g) & k_evt >= 0,
      tau_k := (post_d[unit] - baseline_d[unit]) *
        c(0.3, 0.5, 0.6, 0.5)[pmin(k_evt + 1L, 4L)]]
panel[, Y := unit_fe[unit] + time_fe[period] + tau_k + rnorm(.N, 0, 0.4)]
panel_df <- as.data.frame(panel[order(unit, period), .(unit, period, D, Y)])

cat("first-period D (continuous d_sq):\n"); summary(panel_df$D[panel_df$period == 1L])

cat("\n--- reference WITH continuous = 1 ---\n")
ref <- tryCatch(
  did_multiplegt_dyn(df = panel_df, outcome = "Y", group = "unit",
                     time = "period", treatment = "D",
                     effects = 3, placebo = 0, graph_off = TRUE,
                     continuous = 1),
  error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(ref)) print(as.numeric(ref$results$Effects[, 1]))

cat("\n--- didgpu (r backend, no continuous flag) ---\n")
us <- tryCatch(
  didgpu(panel_df, "Y", "unit", "period", "D",
          effects = 3L, placebo = 0L,
          bootstrap_reps = 0L, backend = "r", verbose = FALSE),
  error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(us)) print(as.numeric(us$results$Effects[, "Estimate"]))
