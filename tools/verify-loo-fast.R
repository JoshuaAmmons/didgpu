# Verify the fast cohort-LOO path (re-aggregation, no refit) is
# bit-identical to the refit-based LOO for never-treated controls, and
# measure the speedup. Also confirm notyet controls fall back to refit.
suppressPackageStartupMessages({ library(didgpu) })

make_panel <- function(n_units, n_periods, seed = 17L) {
  p <- didgpu_simulate_panel(n_units = n_units, n_periods = n_periods,
                              tau_profile = c(0.5, 1.0), seed = seed)
  p$D <- as.integer(p$D >= 0.5); p
}

# Reference: explicit refit-based cohort LOO (what the generic loop does).
refit_loo <- function(p, args) {
  d <- data.table::as.data.table(p)
  d[, F_g := if (any(D == 1L)) min(period[D == 1L]) else Inf, by = unit]
  cohorts <- sort(unique(d$F_g[is.finite(d$F_g)]))
  full <- didgpu_cs(p, "Y","unit","period","D", est_method = args$est_method,
                    control_group = args$control_group,
                    aggregation = args$aggregation, covariates = args$covariates,
                    bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  full_est <- full$aggregation$estimate[1]
  ests <- sapply(cohorts, function(g) {
    keep_units <- unique(d$unit[d$F_g != g | !is.finite(d$F_g)])
    pm <- p[p$unit %in% keep_units, , drop = FALSE]
    fb <- didgpu_cs(pm, "Y","unit","period","D", est_method = args$est_method,
                    control_group = args$control_group,
                    aggregation = args$aggregation, covariates = args$covariates,
                    bootstrap_reps = 0L, backend = "r", verbose = FALSE)
    fb$aggregation$estimate[1]
  })
  data.frame(cohort = cohorts, estimate = ests)
}

cat("=== Correctness: fast (re-aggregate) vs refit, never controls ===\n")
for (agg in c("overall", "event", "group")) {
  p <- make_panel(80L, 12L)
  fit <- didgpu_cs(p, "Y","unit","period","D", est_method = "OR",
                   control_group = "never", aggregation = agg,
                   bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  fast <- didgpu_loo(fit, by = "cohort", df = p, verbose = FALSE)
  ref  <- refit_loo(p, fit$args)
  # align by leave_out / cohort
  fast_ord <- fast[order(as.numeric(fast$leave_out)), ]
  ref_ord  <- ref[order(ref$cohort), ]
  maxdiff <- max(abs(fast_ord$estimate - ref_ord$estimate), na.rm = TRUE)
  cat(sprintf("  aggregation=%-8s  max|fast - refit| = %.3e  method=%s\n",
              agg, maxdiff, attr(fast, "method")))
}

cat("\n=== notyet controls must fall back to refit (method != reaggregate) ===\n")
p <- make_panel(80L, 12L)
fit_ny <- didgpu_cs(p, "Y","unit","period","D", est_method = "OR",
                    control_group = "notyet", aggregation = "overall",
                    bootstrap_reps = 0L, backend = "r", verbose = FALSE)
loo_ny <- didgpu_loo(fit_ny, by = "cohort", df = p, verbose = FALSE)
cat(sprintf("  notyet method = '%s' (should be the refit path)\n",
            attr(loo_ny, "method") %||% "refit (generic)"))

cat("\n=== Speedup: fast vs refit ===\n")
for (nu in c(80L, 200L)) {
  p <- make_panel(nu, 12L)
  fit <- didgpu_cs(p, "Y","unit","period","D", est_method = "OR",
                   control_group = "never", aggregation = "overall",
                   bootstrap_reps = 0L, backend = "r", verbose = FALSE)
  t0 <- Sys.time(); invisible(didgpu_loo(fit, by="cohort", df=p, verbose=FALSE))
  t_fast <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t0 <- Sys.time(); invisible(refit_loo(p, fit$args))
  t_refit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  n_units=%3d  fast=%.4fs  refit=%.3fs  speedup=%.0fx\n",
              nu, t_fast, t_refit, t_refit / t_fast))
}
