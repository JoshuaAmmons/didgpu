# Instrument .core_one_event_time to dump intermediates.
library(didgpu); library(data.table)

# Tiny panel: 3 out-switchers, balanced, no noise, known DGP.
# baseline = 1, F_g = 3, 4, 5 respectively. T_max = 7.
# tau_out = -0.5, -0.8 (k=0 and k=1 effects)
panel <- data.table(
  unit = rep(1:3, each = 7),
  period = rep(1:7, 3),
  D = c(rep(1, 2), rep(0, 5),   # unit 1: F_g=3
        rep(1, 3), rep(0, 4),   # unit 2: F_g=4
        rep(1, 4), rep(0, 3)),  # unit 3: F_g=5
  Y = c(
    # unit 1: F_g=3, pre Y=0, post Y=-0.5, -0.8 (and beyond -1.0)
    0, 0, -0.5, -0.8, -1.0, -1.0, -1.0,
    # unit 2: F_g=4
    0, 0, 0, -0.5, -0.8, -1.0, -1.0,
    # unit 3: F_g=5
    0, 0, 0, 0, -0.5, -0.8, -1.0
  )
)
cat("Panel:\n"); print(panel)

prep <- didgpu:::.prep_panel(panel, "Y", "unit", "period", "D")
cat("\nprep columns:", paste(names(prep), collapse = ", "), "\n")
# Print one row per group to verify per-group quantities.
g_summary <- unique(prep[, list(group_XX, d_sq_XX, F_g_XX, S_g_XX, T_g_XX, L_g_XX, L_g_placebo_XX)])
cat("\nper-group summary:\n"); print(g_summary)

cat("\n=== reference on this panel ===\n")
library(DIDmultiplegtDYN)
ref <- suppressMessages(suppressWarnings(did_multiplegt_dyn(
  df = as.data.frame(panel), outcome = "Y", group = "unit",
  time = "period", treatment = "D",
  effects = 2, placebo = 0, graph_off = TRUE)))
cat("ref Effects:\n"); print(ref$results$Effects)

cat("\n=== k=1, direction=0 (out) ===\n")
res <- didgpu:::.core_one_event_time(prep, k = 1L, direction = 0L)
cat("att =", res$att, "\n")
cat("N_inc =", res$N_inc, "\n")
cat("U_g (per group):", res$U_g, "\n")
cat("sum(U_g)/G =", sum(res$U_g) / length(res$U_g), "\n")

# Manual expected calculation:
# Each unit is its own switcher at k=1 row: t = F_g + 0 = F_g
# Actually wait, dist mask uses t == F_g + k - 1 = F_g for k=1. So switcher
# row is AT F_g (the switch period).
# Wait that's the switch period itself, where D first becomes 0. The diff_y
# at that row is Y(F_g) - Y(F_g - 1) = 0 - 0 = -0.5 - 0 = -0.5 (for unit 1
# with F_g=3, Y(2)=0, Y(3)=-0.5). diff_y = -0.5.
# For k=1 at switcher rows: diff_y = -0.5 for unit 1, -0.5 for unit 2, -0.5
# for unit 3 (all pre-switch Y=0, post-switch first period Y=-0.5).
#
# G = 3 unique groups; N_inc = 3 switcher cells (each unit contributes 1).
# (G/N_inc) = 1.
#
# Switcher row contribution: 1 * 1 * 1 * (1 - 0) * (-0.5) = -0.5
# Sum over 3 switcher rows = -1.5
#
# Control rows: for each unit at pre-switch times, never_change=1.
# unit 1 (F_g=3) is pre-switch at t=1,2. At t=2, diff_y_1 = Y(2)-Y(1) = 0-0=0.
# unit 2 (F_g=4) is pre-switch at t=1,2,3. At t=3, diff_y_1 = Y(3)-Y(2)=0-0=0.
# unit 3 (F_g=5) is pre-switch at t=1,2,3,4. At t=4, diff_y_1 = Y(4)-Y(3)=0-0=0.
#
# Control contributions are all 0 (diff_y = 0 for unobserved-switch controls).
# So sum = -1.5. DID = -1.5 / 3 = -0.5. ✓
#
# Reference (sign-flipped) would report 0.5.
