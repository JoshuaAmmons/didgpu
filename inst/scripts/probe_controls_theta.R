# Instrument the reference to dump its coefs_sq_l_XX, then compare to
# .prefit_controls()'s theta.
library(didgpu); library(DIDmultiplegtDYN); library(data.table)

set.seed(11)
n_units <- 80L; n_periods <- 15L
F_g <- rep(Inf, n_units); treated <- sort(sample(seq_len(n_units), 48L))
F_g[treated] <- sample(5L:10L, 48L, replace = TRUE)
unit_fe <- rnorm(n_units, 0, 1); time_fe <- rnorm(n_periods, 0, 0.3)
X_unit <- rnorm(n_units, 0, 1); X_time <- rnorm(n_periods, 0, 0.5)
panel <- data.table(unit = rep(1:n_units, each = n_periods),
                     period = rep(1:n_periods, n_units))
panel[, F_g := F_g[unit]]; panel[, D := as.integer(period >= F_g & is.finite(F_g))]
panel[, X := X_unit[unit] + X_time[period] + rnorm(.N, 0, 0.2)]
panel[, k_evt := period - F_g]
tau <- c(0.5, 1.0, 1.2, 1.0); panel[, tau_k := 0]
panel[is.finite(F_g) & k_evt >= 0, tau_k := tau[pmin(k_evt + 1L, length(tau))]]
panel[, Y := unit_fe[unit] + time_fe[period] + tau_k + 0.7 * X + rnorm(.N, 0, 0.4)]
panel_df <- as.data.frame(panel[order(unit, period), .(unit, period, D, Y, X)])

# Trace approach: insert assignment to a global env from inside the
# reference's call. The reference assigns coefs_sq_l_XX into its main
# function's environment, then passes it via controls_globals to the
# core. We trace did_multiplegt_main and capture coefs_sq_1_XX (level l=1).
trace_env <- new.env()
trace(DIDmultiplegtDYN:::did_multiplegt_main, exit = quote({
  for (nm in ls(environment())) {
    if (grepl("^coefs_sq_", nm) || grepl("^useful_res_", nm) ||
        grepl("^inv_Denom_", nm)) {
      assign(nm, get(nm), envir = trace_env)
    }
  }
  # Also capture the regression inputs for l=1: Y_vec and X_mat
  if (exists("Y_vec", inherits = FALSE)) {
    assign("Y_vec", Y_vec, envir = trace_env)
  }
  if (exists("X_mat", inherits = FALSE)) {
    assign("X_mat", X_mat, envir = trace_env)
  }
  if (exists("data_XX", inherits = FALSE)) {
    assign("data_XX_rows", nrow(data_XX), envir = trace_env)
    if ("diff_y_wXX" %in% names(data_XX)) {
      assign("ref_Y_diff_y_w", data_XX$diff_y_wXX, envir = trace_env)
    }
    keep_cols <- intersect(c("group_XX", "time_XX", "diff_y_XX", "diff_X1_XX",
                              "avg_diff_X1_XX", "resid_X1_time_FE_XX", "N_gt_XX"),
                            names(data_XX))
    if (length(keep_cols)) {
      assign("ref_data_subset", data_XX[, keep_cols, with = FALSE],
             envir = trace_env)
    }
  }
}), print = FALSE)

ref <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = panel_df, outcome = "Y", group = "unit", time = "period",
  treatment = "D", effects = 3, placebo = 0, graph_off = TRUE,
  controls = "X")))
untrace(DIDmultiplegtDYN:::did_multiplegt_main)

cat("--- captured from reference's main() ---\n")
for (nm in ls(trace_env)) {
  cat(nm, ":\n")
  print(get(nm, envir = trace_env))
  cat("\n")
}

# Now compute MY theta via .prefit_controls.
prepped <- didgpu:::.prep_panel(panel_df, "Y", "unit", "period", "D", controls = "X")
prefit  <- didgpu:::.prefit_controls(prepped, controls = "X")
cat("\n--- didgpu .prefit_controls()$theta ---\n")
print(prefit$theta)

# Compare the regression inputs row-counts
cat("\n--- reference data_XX rows used in regression: ",
    trace_env$data_XX_rows, "---\n")
# How many rows do I use?
T_max_XX <- max(prepped$time_XX)
prepped[, ever_change_d_XX := as.integer(F_g_XX <= T_max_XX)]
prepped[, fd_X_all_non_missing_XX := as.integer(!is.na(get("diff_X_1_XX")))]
my_mask <- prepped$ever_change_d_XX == 0L &
           prepped$d_sq_XX == 0L &
           !is.na(prepped$diff_y_XX) &
           prepped$fd_X_all_non_missing_XX == 1L &
           prepped$N_gt_XX > 0
cat("--- didgpu rows used in regression for l=0: ", sum(my_mask), "---\n")

if (!is.null(trace_env$ref_data_subset)) {
  cat("\n--- first 10 rows of reference's regression data ---\n")
  print(head(trace_env$ref_data_subset, 10))
}
my_reg <- prepped[my_mask, list(group_XX, time_XX, diff_y_XX, diff_X1 = get("diff_X_1_XX"))]
cat("\n--- first 10 rows of MY regression data ---\n")
print(head(my_reg, 10))
