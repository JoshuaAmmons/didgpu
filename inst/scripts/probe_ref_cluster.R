library(DIDmultiplegtDYN); library(didgpu)
set.seed(11)
p <- didgpu_simulate_panel_bidir(n_units = 100L, n_periods = 18L,
                                  frac_treated = 0.6, frac_in = 0.5,
                                  min_treat_period = 6L, max_treat_period = 12L,
                                  seed = 11L)
p$cluster <- ((p$unit - 1L) %/% 4L) + 1L
cat("--- ref on ORIGINAL panel with cluster arg ---\n")
fit <- tryCatch(
  did_multiplegt_dyn(df = as.data.frame(p), outcome = "Y", group = "unit",
                     time = "period", treatment = "D", cluster = "cluster",
                     effects = 4, placebo = 2, graph_off = TRUE),
  error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(fit)) {
  cat("ref OK with cluster\n"); print(as.numeric(fit$results$Effects[, 1]))
}

cat("\n--- ref on ORIGINAL panel WITHOUT cluster (effects=4 only) ---\n")
fit2 <- tryCatch(
  did_multiplegt_dyn(df = as.data.frame(p), outcome = "Y", group = "unit",
                     time = "period", treatment = "D",
                     effects = 4, placebo = 2, graph_off = TRUE),
  error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(fit2)) {
  cat("ref OK without cluster\n"); print(as.numeric(fit2$results$Effects[, 1]))
}

cat("\n--- ref on ORIGINAL panel with cluster + placebo=0 ---\n")
fit3 <- tryCatch(
  did_multiplegt_dyn(df = as.data.frame(p), outcome = "Y", group = "unit",
                     time = "period", treatment = "D", cluster = "cluster",
                     effects = 4, placebo = 0, graph_off = TRUE),
  error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(fit3)) {
  cat("ref OK with cluster + placebo=0\n"); print(as.numeric(fit3$results$Effects[, 1]))
}
