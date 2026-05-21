library(didgpu); library(data.table)
p <- as.data.table(didgpu_simulate_panel_bidir(
  n_units = 80L, n_periods = 15L, frac_treated = 0.6, frac_in = 0.5,
  min_treat_period = 5L, max_treat_period = 10L, seed = 11L
))
truth <- attr(didgpu_simulate_panel_bidir(
  n_units = 80L, n_periods = 15L, frac_treated = 0.6, frac_in = 0.5,
  min_treat_period = 5L, max_treat_period = 10L, seed = 11L
), "truth")
p[, F_g := truth$F_g[as.character(unit)]]
p[, direction := truth$direction[as.character(unit)]]
p[, k_evt := period - F_g]
p[!is.finite(F_g), k_evt := NA]
cat("Mean Y by (direction, event_time) for switchers:\n")
print(p[!is.na(direction) & k_evt >= -3 & k_evt <= 3,
        list(mean_Y = mean(Y), n = .N), by = list(direction, k_evt)][order(direction, k_evt)])

cat("\nFirst-period D by direction:\n")
print(p[period == 1L, list(unit_fe_mean = mean(Y)), by = list(direction, D)])
