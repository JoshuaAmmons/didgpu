# Test if a simple FWL-style residualization of diff_y by control vars
# reproduces the reference's controls output. If yes, controls is easy.
library(didgpu); library(DIDmultiplegtDYN); library(data.table)

# Panel with one control X that correlates with both Y and treatment.
set.seed(11)
n_units <- 80L; n_periods <- 15L
F_g <- rep(Inf, n_units)
treated <- sort(sample(seq_len(n_units), 48L))
F_g[treated] <- sample(5L:10L, 48L, replace = TRUE)
unit_fe <- rnorm(n_units, 0, 1)
time_fe <- rnorm(n_periods, 0, 0.3)
X_unit  <- rnorm(n_units, 0, 1)            # time-invariant control trait
X_time  <- rnorm(n_periods, 0, 0.5)        # time-varying part

panel <- data.table(unit = rep(1:n_units, each = n_periods),
                     period = rep(1:n_periods, n_units))
panel[, F_g := F_g[unit]]
panel[, D := as.integer(period >= F_g & is.finite(F_g))]
panel[, X := X_unit[unit] + X_time[period] + rnorm(.N, 0, 0.2)]
panel[, k_evt := period - F_g]
tau <- c(0.5, 1.0, 1.2, 1.0)
panel[, tau_k := 0]
panel[is.finite(F_g) & k_evt >= 0,
      tau_k := tau[pmin(k_evt + 1L, length(tau))]]
panel[, Y := unit_fe[unit] + time_fe[period] + tau_k +
            0.7 * X +                       # X affects Y
            rnorm(.N, 0, 0.4)]
panel <- panel[order(unit, period), .(unit, period, D, Y, X)]

cat("=== ref WITHOUT controls ===\n")
ref_nc <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(panel), outcome = "Y", group = "unit",
  time = "period", treatment = "D",
  effects = 3, placebo = 1, graph_off = TRUE)))
cat("  Effects (no controls):", round(as.numeric(ref_nc$results$Effects[, 1]), 4), "\n")

cat("\n=== ref WITH controls = X ===\n")
ref_c <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(panel), outcome = "Y", group = "unit",
  time = "period", treatment = "D",
  effects = 3, placebo = 1, graph_off = TRUE,
  controls = "X")))
cat("  Effects (with controls X):", round(as.numeric(ref_c$results$Effects[, 1]), 4), "\n")

cat("\n=== didgpu WITHOUT controls (the simple FWL we have today) ===\n")
us_nc <- didgpu(as.data.frame(panel), "Y", "unit", "period", "D",
                 effects = 3L, placebo = 1L, bootstrap_reps = 0L,
                 backend = "r", verbose = FALSE)
cat("  Effects:", round(as.numeric(us_nc$results$Effects[, "Estimate"]), 4), "\n")

cat("\n=== reference-style: residualize diff_y by diff_X on controls only ===\n")
# Compute first differences per group
panel[, diff_y := Y - shift(Y), by = unit]
panel[, diff_X := X - shift(X), by = unit]
# d_sq = period-1 treatment per group
panel[, d_sq := D[1L], by = unit]
# ever-change indicator
panel[, ever_change := any(D != d_sq), by = unit]
# Regression of diff_y on diff_X among "control units" (ever_change == FALSE).
# Note: reference also residualizes diff_X by per-(time, d_sq) means first.
# Let's try the simpler version (no demean) and see how close we get.
ctrl_rows <- !panel$ever_change & !is.na(panel$diff_y) & !is.na(panel$diff_X)
m <- lm(diff_y ~ diff_X - 1, data = panel[ctrl_rows, ])
cat("  control regression slope: ", coef(m), "\n")
# Subtract from all rows (extrapolate to switchers)
panel[, diff_y_resid := diff_y - coef(m) * diff_X]
# Reconstruct Y_resid by inverse-differencing back from period 1
panel[, Y_resid := Y - cumsum(c(0, head(coef(m) * diff_X, -1))), by = unit]

us_resid2 <- didgpu(as.data.frame(panel), "Y_resid", "unit", "period", "D",
                     effects = 3L, placebo = 1L, bootstrap_reps = 0L,
                     backend = "r", verbose = FALSE)
cat("  Effects (Y - cumsum(coef * diff_X)):",
    round(as.numeric(us_resid2$results$Effects[, "Estimate"]), 4), "\n")

cat("\n--- compare: ref_c (target) vs reference-style didgpu ---\n")
cat("  max diff:",
    max(abs(as.numeric(ref_c$results$Effects[, 1]) -
            as.numeric(us_resid2$results$Effects[, "Estimate"]))), "\n")
