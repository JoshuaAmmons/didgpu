library(didgpu); library(DIDmultiplegtDYN)
set.seed(11)
p <- didgpu_simulate_panel_bidir(n_units = 100L, n_periods = 18L,
                                  frac_treated = 0.6, frac_in = 0.5,
                                  min_treat_period = 6L, max_treat_period = 12L,
                                  seed = 11L)
p$cluster <- ((p$unit - 1L) %/% 4L) + 1L

# What does the cluster-resampled panel look like for one bootstrap iter?
args <- list(group = "unit", time = "period", cluster = "cluster")
df_boot <- didgpu:::.cluster_resample(p, args, iter_seed = 1L)
cat("resampled rows:", nrow(df_boot), "\n")
cat("unique units:", length(unique(df_boot$unit)), "\n")
cat("unique clusters:", length(unique(df_boot$cluster)), "\n")
cat("unit range: [", min(df_boot$unit), ",", max(df_boot$unit), "]\n")
cat("cluster range: [", min(df_boot$cluster), ",", max(df_boot$cluster), "]\n")
cat("any dup (unit, period)?", anyDuplicated(df_boot[, c("unit", "period")]) > 0L, "\n")

# Run the reference on the resampled panel directly with cluster arg.
cat("\n--- ref on resampled panel ---\n")
res <- tryCatch({
  fit <- did_multiplegt_dyn(
    df = as.data.frame(df_boot), outcome = "Y", group = "unit",
    time = "period", treatment = "D", cluster = "cluster",
    effects = 4, placebo = 2, graph_off = TRUE)
  cat("ref ran OK\n"); print(as.numeric(fit$results$Effects[, 1]))
}, error = function(e) cat("FAILED:", conditionMessage(e), "\n"))
