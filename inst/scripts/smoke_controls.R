# Test the newly-wired controls= arg in didgpu(backend = "r").
library(didgpu); library(DIDmultiplegtDYN); library(data.table)

# Panel with one control X that correlates with both Y and treatment.
set.seed(11)
n_units <- 80L; n_periods <- 15L
F_g <- rep(Inf, n_units)
treated <- sort(sample(seq_len(n_units), 48L))
F_g[treated] <- sample(5L:10L, 48L, replace = TRUE)
unit_fe <- rnorm(n_units, 0, 1)
time_fe <- rnorm(n_periods, 0, 0.3)
X_unit  <- rnorm(n_units, 0, 1)
X_time  <- rnorm(n_periods, 0, 0.5)

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
            0.7 * X +
            rnorm(.N, 0, 0.4)]
panel <- panel[order(unit, period), .(unit, period, D, Y, X)]

cat("--- ref WITHOUT controls ---\n")
ref_nc <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(panel), outcome = "Y", group = "unit",
  time = "period", treatment = "D",
  effects = 3, placebo = 1, graph_off = TRUE)))
cat("  effects:", round(as.numeric(ref_nc$results$Effects[, 1]), 6), "\n")
cat("  placebos:", round(as.numeric(ref_nc$results$Placebos[, 1]), 6), "\n")

cat("\n--- ref WITH controls = X ---\n")
ref_c <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(panel), outcome = "Y", group = "unit",
  time = "period", treatment = "D",
  effects = 3, placebo = 1, graph_off = TRUE,
  controls = "X")))
cat("  effects:", round(as.numeric(ref_c$results$Effects[, 1]), 6), "\n")
cat("  placebos:", round(as.numeric(ref_c$results$Placebos[, 1]), 6), "\n")

cat("\n--- didgpu WITH controls = 'X' ---\n")
us_c <- didgpu(as.data.frame(panel), "Y", "unit", "period", "D",
                effects = 3L, placebo = 1L, controls = "X",
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
cat("  effects:", round(as.numeric(us_c$results$Effects[, "Estimate"]), 6), "\n")
cat("  placebos:", round(as.numeric(us_c$results$Placebos[, "Estimate"]), 6), "\n")

cat("\n--- diffs ---\n")
cat("  effects max diff:",
    max(abs(as.numeric(ref_c$results$Effects[, 1]) -
            as.numeric(us_c$results$Effects[, "Estimate"]))), "\n")
cat("  placebos max diff:",
    max(abs(as.numeric(ref_c$results$Placebos[, 1]) -
            as.numeric(us_c$results$Placebos[, "Estimate"]))), "\n")
