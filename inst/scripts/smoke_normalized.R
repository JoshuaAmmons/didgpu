library(didgpu); library(DIDmultiplegtDYN)

p <- didgpu_simulate_panel(n_units = 80L, n_periods = 15L,
                            frac_treated = 0.6,
                            min_treat_period = 5L, max_treat_period = 10L,
                            tau_profile = c(0.5, 1.0, 1.2),
                            sigma = 0.4, seed = 17L)

cat("=== binary, both directions, normalized = TRUE ===\n")
ref_y <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(p), outcome = "Y", group = "unit",
  time = "period", treatment = "D",
  effects = 3, placebo = 0, graph_off = TRUE, normalized = TRUE)))
us_y <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, normalized = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
cat("ref:    "); print(as.numeric(ref_y$results$Effects[, 1]))
cat("didgpu: "); print(as.numeric(us_y$results$Effects[, "Estimate"]))
cat("max diff:",
    max(abs(as.numeric(ref_y$results$Effects[, 1]) -
            as.numeric(us_y$results$Effects[, "Estimate"]))), "\n")

cat("\n=== binary, switchers = in, normalized = TRUE ===\n")
ref_in <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(p), outcome = "Y", group = "unit",
  time = "period", treatment = "D",
  effects = 3, placebo = 0, graph_off = TRUE,
  switchers = "in", normalized = TRUE)))
us_in <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, normalized = TRUE,
                switchers = "in",
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
cat("ref:    "); print(as.numeric(ref_in$results$Effects[, 1]))
cat("didgpu: "); print(as.numeric(us_in$results$Effects[, "Estimate"]))
cat("max diff:",
    max(abs(as.numeric(ref_in$results$Effects[, 1]) -
            as.numeric(us_in$results$Effects[, "Estimate"]))), "\n")

cat("\n=== binary, placebo = 2, normalized = TRUE ===\n")
ref_pl <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(p), outcome = "Y", group = "unit",
  time = "period", treatment = "D",
  effects = 3, placebo = 2, graph_off = TRUE, normalized = TRUE)))
us_pl <- didgpu(p, "Y", "unit", "period", "D",
                effects = 3L, placebo = 2L, normalized = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
cat("ref placebos:    "); print(as.numeric(ref_pl$results$Placebos[, 1]))
cat("didgpu placebos: "); print(as.numeric(us_pl$results$Placebos[, "Estimate"]))
cat("max diff (pl):",
    max(abs(as.numeric(ref_pl$results$Placebos[, 1]) -
            as.numeric(us_pl$results$Placebos[, "Estimate"]))), "\n")

cat("\n=== continuous = 1, normalized = TRUE ===\n")
set.seed(11)
n_units <- 80L; n_periods <- 15L
F_g <- rep(Inf, n_units)
treated <- sort(sample(seq_len(n_units), 48L))
F_g[treated] <- sample(5L:10L, 48L, replace = TRUE)
baseline_d <- runif(n_units, 0, 1)
post_d <- baseline_d + rnorm(n_units, 0, 0.5)
unit_fe <- rnorm(n_units, 0, 1); time_fe <- rnorm(n_periods, 0, 0.3)
pan <- data.table::data.table(unit = rep(1:n_units, each = n_periods),
                               period = rep(1:n_periods, n_units))
pan[, F_g := F_g[unit]]
pan[, D := ifelse(period >= F_g, post_d[unit], baseline_d[unit])]
pan[, k_evt := period - F_g]
pan[, tau_k := 0]
pan[is.finite(F_g) & k_evt >= 0,
    tau_k := (post_d[unit] - baseline_d[unit]) *
      c(0.3, 0.5, 0.6, 0.5)[pmin(k_evt + 1L, 4L)]]
pan[, Y := unit_fe[unit] + time_fe[period] + tau_k + rnorm(.N, 0, 0.4)]
pan_df <- as.data.frame(pan[order(unit, period), .(unit, period, D, Y)])

ref_c <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = pan_df, outcome = "Y", group = "unit", time = "period", treatment = "D",
  effects = 3, placebo = 0, graph_off = TRUE,
  continuous = 1, normalized = TRUE)))
us_c <- didgpu(pan_df, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L,
                continuous = 1L, normalized = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
cat("ref:    "); print(as.numeric(ref_c$results$Effects[, 1]))
cat("didgpu: "); print(as.numeric(us_c$results$Effects[, "Estimate"]))
cat("max diff:",
    max(abs(as.numeric(ref_c$results$Effects[, 1]) -
            as.numeric(us_c$results$Effects[, "Estimate"]))), "\n")

cat("\n=== multivalued (D in {0,1,2,3}), normalized = TRUE ===\n")
set.seed(23)
mv <- data.table::data.table(
  unit = rep(1:80, each = 15),
  period = rep(1:15, 80))
mv[, F_g := {
  fg <- rep(Inf, 80)
  treated <- sort(sample(1:80, 48))
  fg[treated] <- sample(5:10, 48, replace = TRUE)
  fg[unit]
}]
mv[, dose := {
  doses <- sample(c(1L, 2L, 3L), 80, replace = TRUE)
  doses[unit]
}]
mv[, D := as.integer(ifelse(period >= F_g, dose, 0L))]
mv[, k_evt := period - F_g]
mv[, tau_k := 0]
mv[is.finite(F_g) & k_evt >= 0,
   tau_k := dose * c(0.3, 0.5, 0.6, 0.5)[pmin(k_evt + 1L, 4L)]]
mv_unit_fe <- rnorm(80, 0, 1)
mv_time_fe <- rnorm(15, 0, 0.3)
mv[, Y := mv_unit_fe[unit] + mv_time_fe[period] + tau_k + rnorm(.N, 0, 0.4)]
mv_df <- as.data.frame(mv[order(unit, period), .(unit, period, D, Y)])

ref_mv <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = mv_df, outcome = "Y", group = "unit", time = "period", treatment = "D",
  effects = 3, placebo = 0, graph_off = TRUE, normalized = TRUE)))
us_mv <- didgpu(mv_df, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, normalized = TRUE,
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
cat("ref:    "); print(as.numeric(ref_mv$results$Effects[, 1]))
cat("didgpu: "); print(as.numeric(us_mv$results$Effects[, "Estimate"]))
cat("max diff (mv):",
    max(abs(as.numeric(ref_mv$results$Effects[, 1]) -
            as.numeric(us_mv$results$Effects[, "Estimate"]))), "\n")
