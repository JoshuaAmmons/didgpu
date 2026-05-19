library(didgpu); library(DIDmultiplegtDYN)

# Panel with ONLY out-switchers + always-treated controls (no in-switchers, no
# never-untreated). All units start treated; some turn off, some stay on.
set.seed(11)
n_units <- 80L; n_periods <- 15L
n_treated <- 48L  # of those, all are out-switchers
treated_units <- sort(sample(seq_len(n_units), n_treated))
F_g <- rep(Inf, n_units)
F_g[treated_units] <- sample(5L:10L, n_treated, replace = TRUE)
unit_fe <- rnorm(n_units, 0, 1)
time_fe <- rnorm(n_periods, 0, 0.3)
tau_out <- c(-0.5, -0.8, -1.0, -1.1, -1.2)

panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
panel$F_g <- F_g[panel$unit]
panel$D <- ifelse(panel$period >= panel$F_g, 0L, 1L)  # baseline = 1 for all, turn off at F_g
panel$k_evt <- panel$period - panel$F_g
panel$tau_k <- 0
post <- is.finite(panel$F_g) & panel$k_evt >= 0
panel$tau_k[post] <- tau_out[pmin(panel$k_evt[post] + 1L, length(tau_out))]
panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] + panel$tau_k + rnorm(nrow(panel), 0, 0.4)
panel <- panel[order(panel$unit, panel$period), c("unit", "period", "D", "Y")]

cat("D distribution at t=1:\n"); print(table(panel$D[panel$period == 1L]))

cat("\n--- ref switchers=out (only direction; expect positive values) ---\n")
ref <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = panel, outcome = "Y", group = "unit", time = "period",
  treatment = "D", effects = 3, placebo = 0, graph_off = TRUE)))
print(as.numeric(ref$results$Effects[, 1]))

cat("\n--- didgpu r backend ---\n")
us <- didgpu(panel, "Y", "unit", "period", "D",
              effects = 3L, placebo = 0L, bootstrap_reps = 0L,
              backend = "r", verbose = FALSE)
print(as.numeric(us$results$Effects[, "Estimate"]))

cat("\nmax diff:", max(abs(as.numeric(ref$results$Effects[, 1]) -
                            as.numeric(us$results$Effects[, "Estimate"]))), "\n")
