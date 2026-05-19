library(didgpu); library(DIDmultiplegtDYN); library(data.table)

set.seed(11)
n_units <- 100L; n_periods <- 15L
F_g <- rep(Inf, n_units)
treated <- sort(sample(seq_len(n_units), 60L))
F_g[treated] <- sample(5L:10L, 60L, replace = TRUE)
unit_fe <- rnorm(n_units, 0, 1); time_fe <- rnorm(n_periods, 0, 0.3)

# Each unit belongs to one of 3 industries; industry-specific time trend.
industry <- sample.int(3L, n_units, replace = TRUE)
industry_trend <- matrix(rnorm(3L * n_periods, 0, 0.4),
                          nrow = 3L, ncol = n_periods)

panel <- data.table(unit = rep(1:n_units, each = n_periods),
                     period = rep(1:n_periods, n_units))
panel[, industry := industry[unit]]
panel[, F_g := F_g[unit]]
panel[, D := as.integer(period >= F_g & is.finite(F_g))]
panel[, k_evt := period - F_g]
tau <- c(0.5, 1.0, 1.2, 1.0)
panel[, tau_k := 0]
panel[is.finite(F_g) & k_evt >= 0, tau_k := tau[pmin(k_evt + 1L, length(tau))]]
# Y includes the industry-time interaction.
panel[, ind_time := industry_trend[cbind(industry, period)]]
panel[, Y := unit_fe[unit] + time_fe[period] + tau_k +
            ind_time + rnorm(.N, 0, 0.4)]
panel_df <- as.data.frame(panel[order(unit, period),
                                  .(unit, period, D, Y, industry)])

cat("--- ref WITHOUT trends_nonparam ---\n")
ref_nt <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = panel_df, outcome = "Y", group = "unit", time = "period",
  treatment = "D", effects = 3, placebo = 0, graph_off = TRUE)))
print(as.numeric(ref_nt$results$Effects[, 1]))

cat("\n--- ref WITH trends_nonparam = 'industry' ---\n")
ref_t <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = panel_df, outcome = "Y", group = "unit", time = "period",
  treatment = "D", effects = 3, placebo = 0, graph_off = TRUE,
  trends_nonparam = "industry")))
print(as.numeric(ref_t$results$Effects[, 1]))

cat("\n--- didgpu WITH trends_nonparam = 'industry' ---\n")
us <- didgpu(panel_df, "Y", "unit", "period", "D",
              effects = 3L, placebo = 0L, trends_nonparam = "industry",
              bootstrap_reps = 0L, backend = "r", verbose = FALSE)
print(as.numeric(us$results$Effects[, "Estimate"]))

cat("\nmax diff:",
    max(abs(as.numeric(ref_t$results$Effects[, 1]) -
            as.numeric(us$results$Effects[, "Estimate"]))), "\n")
