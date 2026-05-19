# Test that cluster argument (distinct from group) works.
library(didgpu); library(DIDmultiplegtDYN); library(data.table)

set.seed(11)
n_units <- 80L; n_periods <- 15L
# Each unit is in one of 20 clusters (so cluster nests group).
clust <- sample.int(20L, n_units, replace = TRUE)
F_g <- rep(Inf, n_units)
treated <- sort(sample(seq_len(n_units), 48L))
F_g[treated] <- sample(5L:10L, 48L, replace = TRUE)
unit_fe <- rnorm(n_units, 0, 1)
time_fe <- rnorm(n_periods, 0, 0.3)

panel <- expand.grid(unit = 1:n_units, period = 1:n_periods)
panel$cluster <- clust[panel$unit]
panel$F_g <- F_g[panel$unit]
panel$D <- as.integer(panel$period >= panel$F_g & is.finite(panel$F_g))
panel$k_evt <- panel$period - panel$F_g
panel$tau_k <- 0
post <- is.finite(panel$F_g) & panel$k_evt >= 0
panel$tau_k[post] <- c(0.5, 1.0, 1.2, 1.0)[pmin(panel$k_evt[post] + 1L, 4L)]
panel$Y <- unit_fe[panel$unit] + time_fe[panel$period] +
           panel$tau_k + rnorm(nrow(panel), 0, 0.4)
panel <- panel[order(panel$unit, panel$period), ]

cat("--- point estimate, cluster=NULL ---\n")
us_no <- didgpu(panel, "Y", "unit", "period", "D",
                 effects = 3L, placebo = 0L, bootstrap_reps = 0L,
                 backend = "r", verbose = FALSE)
cat("  ", round(as.numeric(us_no$results$Effects[, "Estimate"]), 4), "\n")

cat("\n--- point estimate, cluster='cluster' (should be same; clusters don't affect the point estimate) ---\n")
us_c <- didgpu(panel, "Y", "unit", "period", "D",
                effects = 3L, placebo = 0L, cluster = "cluster",
                bootstrap_reps = 0L, backend = "r", verbose = FALSE)
cat("  ", round(as.numeric(us_c$results$Effects[, "Estimate"]), 4), "\n")
cat("  diff vs no-cluster: ", max(abs(as.numeric(us_no$results$Effects[, "Estimate"]) -
                                       as.numeric(us_c$results$Effects[, "Estimate"]))), "\n")

cat("\n--- bootstrap, cluster=NULL (cluster by group) ---\n")
us_b_g <- didgpu(panel, "Y", "unit", "period", "D",
                  effects = 3L, placebo = 0L, bootstrap_reps = 30L, seed = 1L,
                  backend = "r", verbose = FALSE)
cat("  Effects:\n"); print(us_b_g$results$Effects)

cat("\n--- bootstrap, cluster='cluster' ---\n")
us_b_c <- didgpu(panel, "Y", "unit", "period", "D",
                  effects = 3L, placebo = 0L, cluster = "cluster",
                  bootstrap_reps = 30L, seed = 1L,
                  backend = "r", verbose = FALSE)
cat("  Effects:\n"); print(us_b_c$results$Effects)

cat("\n  SEs differ between group-cluster and cluster-cluster bootstrap?\n")
cat("    group-cluster:   ", round(us_b_g$results$Effects[, "SE"], 4), "\n")
cat("    custom-cluster:  ", round(us_b_c$results$Effects[, "SE"], 4), "\n")
