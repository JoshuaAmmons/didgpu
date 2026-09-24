# Analytic standard errors for the dCDH estimator.
#
# DIDmultiplegtDYN reports SEs from the estimator's asymptotic linear
# representation in a single pass. didgpu used to get them by bootstrap,
# which was both a different number (the tests said so outright: "SEs come
# from different estimators (bootstrap vs. analytic) and are not compared
# here") and ~100x the work, since the default 100 reps means 101 full
# fits against the reference's one.
#
# The variance kernel is the SAME kernel as the point estimate, with the
# outcome difference replaced by its within-cell residual and a small-cell
# degrees-of-freedom inflation:
#
#     U^var_gt = (G / N_inc)
#                * 1{k+1 <= t <= T_g} * N_gt
#                * (dist_k - (N_t_switch / N_t_control) * never_change_k)
#                * DOF_gt * (diff_y_k - E_hat_gt)
#
# against the estimate's
#
#     U_gt     = (G / N_inc)
#                * 1{k+1 <= t <= T_g} * N_gt
#                * (dist_k - (N_t_switch / N_t_control) * never_change_k)
#                * diff_y_k
#
# E_hat_gt is the estimated mean of diff_y in the cell the row belongs to
# -- switcher cells and not-yet-switched cells are pooled separately --
# and DOF_gt is sqrt(n / (n - 1)) for that cell. Both fall back to a
# combined switcher + non-switcher cell when their own cell is a
# singleton, and to zero / NA when neither is estimable.
#
# Summing U^var by group gives one number per group. The SE is then
#
#     se_k = sqrt( sum_g (U^var_g)^2 / G^2 )
#
# (or, clustered, the same after summing U^var within cluster first).
# Keeping the per-group vectors also gives the full covariance matrix by
# polarisation, which is what the joint nullity tests need:
#
#     cov(k, j) = ( sum_g (U_k + U_j)^2 / G^2 - se_k^2 - se_j^2 ) / 2
#
# Reference: did_multiplegt_dyn_core.R:296-535 (cells, E_hat, DOF and the
# variance kernel) and did_multiplegt_main.R:1005-1060 (pooling and se),
# :1148-1190 (the covariance identity).
#
# Only the default `less_conservative_se = FALSE` branch is implemented,
# which is the only one didgpu exposes.


# The effects and placebo kernels are the same computation over
# differently-named columns; this keeps one implementation of each.
#' @keywords internal
#' @noRd
.se_cols <- function(which = c("effect", "placebo")) {
  which <- match.arg(which)
  if (identical(which, "effect")) {
    list(diff_y = "diff_y_k_XX", dist = "dist_k_XX",
         never = "never_change_k_XX", n_switch = "N_t_switch",
         ratio = "ratio_XX")
  } else {
    list(diff_y = "diff_y_pl_k_XX", dist = "dist_k_pl_XX",
         never = "never_change_k_pl_XX", n_switch = "N_t_switch_pl",
         ratio = "ratio_pl_XX")
  }
}


# Cell means and degrees of freedom for one (event-time, direction).
# Adds scratch columns to `d` by reference; the caller drops them.
#' @keywords internal
#' @noRd
.se_build_cells <- function(d, cohort_cols, cluster_col = NULL,
                             cols = .se_cols()) {
  clustered <- !is.null(cluster_col) && nzchar(cluster_col) &&
                 cluster_col %in% names(d)

  # Rows eligible to inform a cell mean.
  #   ns: a not-yet-switched row in a cohort that has switcher mass now.
  #   s : a switcher row at its own event-time-k cell.
  dy <- cols$diff_y; dk <- cols$dist; nc <- cols$never; nts <- cols$n_switch
  d[, dof_ns_XX := as.integer(!is.na(N_gt_XX) & N_gt_XX != 0 &
                                !is.na(get(dy)) &
                                !is.na(get(nc)) & get(nc) == 1L &
                                !is.na(get(nts)) & get(nts) > 0)]
  d[, dof_s_XX  := as.integer(!is.na(N_gt_XX) & N_gt_XX != 0 &
                                !is.na(get(dk)) & get(dk) == 1L)]
  d[, dof_ns_s_XX := as.integer(dof_s_XX == 1L | dof_ns_XX == 1L)]
  d[, diff_y_N_XX := N_gt_XX * get(dy)]

  # Switcher cells are grouped by the treatment PATH (baseline dose and
  # the dose switched into), not just by period -- reference
  # did_multiplegt_dyn_core.R:326.
  s_cols <- c("d_sq_XX", "F_g_XX", "d_fg_XX",
              setdiff(cohort_cols, c("time_XX", "d_sq_XX")))

  .cell <- function(flag, cnt, tot, mean_c, dofc, by_cols) {
    d[get(flag) == 1L, (cnt) := sum(N_gt_XX, na.rm = TRUE), by = by_cols]
    d[get(flag) == 1L, (tot) := sum(diff_y_N_XX, na.rm = TRUE), by = by_cols]
    d[, (mean_c) := get(tot) / get(cnt)]
    if (clustered) {
      d[, cl_tmp_XX := ifelse(get(flag) == 1L, get(cluster_col), NA)]
      d[!is.na(cl_tmp_XX), (dofc) := data.table::uniqueN(cl_tmp_XX),
        by = by_cols]
      d[, "cl_tmp_XX" := NULL]
    } else {
      d[get(flag) == 1L, (dofc) := sum(get(flag), na.rm = TRUE),
        by = by_cols]
    }
  }

  for (nm in c("cnt_ns_XX", "tot_ns_XX", "dofc_ns_XX",
               "cnt_s_XX", "tot_s_XX", "dofc_s_XX",
               "cnt_nss_XX", "tot_nss_XX", "dofc_nss_XX")) {
    d[, (nm) := NA_real_]
  }
  .cell("dof_ns_XX",   "cnt_ns_XX",  "tot_ns_XX",  "mean_ns_XX",  "dofc_ns_XX",  cohort_cols)
  .cell("dof_s_XX",    "cnt_s_XX",   "tot_s_XX",   "mean_s_XX",   "dofc_s_XX",   s_cols)
  .cell("dof_ns_s_XX", "cnt_nss_XX", "tot_nss_XX", "mean_nss_XX", "dofc_nss_XX", cohort_cols)
  invisible(d)
}


# E_hat_gt: the cell mean that each row's diff_y is measured against.
# The assignment ORDER matters and mirrors the reference statement for
# statement (did_multiplegt_dyn_core.R:448-480), including the two
# unconditional NA sweeps, which are not no-ops.
#' @keywords internal
#' @noRd
.se_build_ehat <- function(d, k) {
  k <- as.integer(k)
  d[, E_hat_XX := NA_real_]
  d[time_XX < F_g_XX | (F_g_XX - 1L + k) == time_XX, E_hat_XX := 0]
  d[time_XX < F_g_XX & (is.na(dofc_ns_XX) | dofc_ns_XX >= 2),
    E_hat_XX := mean_ns_XX]
  d[is.na(mean_ns_XX), E_hat_XX := NA_real_]
  d[(F_g_XX - 1L + k) == time_XX & (is.na(dofc_s_XX) | dofc_s_XX >= 2),
    E_hat_XX := mean_s_XX]
  d[(!is.na(dofc_nss_XX) & dofc_nss_XX >= 2) &
      ((((F_g_XX - 1L + k) == time_XX) & !is.na(dofc_s_XX) & dofc_s_XX == 1) |
       ((time_XX < F_g_XX) & !is.na(dofc_ns_XX) & dofc_ns_XX == 1)),
    E_hat_XX := mean_nss_XX]
  d[is.na(mean_nss_XX), E_hat_XX := NA_real_]
  invisible(d)
}


# DOF_gt: sqrt(n / (n - 1)) for the cell backing this row, with the same
# singleton fallback as E_hat. Reference: core.R:421-447.
#' @keywords internal
#' @noRd
.se_build_dof <- function(d, k) {
  k <- as.integer(k)
  d[, DOF_XX := NA_real_]
  d[time_XX < F_g_XX | (F_g_XX - 1L + k) == time_XX, DOF_XX := 1]
  d[(F_g_XX - 1L + k) == time_XX & dofc_s_XX > 1,
    DOF_XX := sqrt(dofc_s_XX / (dofc_s_XX - 1))]
  d[time_XX < F_g_XX & dofc_ns_XX > 1,
    DOF_XX := sqrt(dofc_ns_XX / (dofc_ns_XX - 1))]
  d[dofc_nss_XX >= 2 &
      ((((F_g_XX - 1L + k) == time_XX) & dofc_s_XX == 1) |
       ((time_XX < F_g_XX) & dofc_ns_XX == 1)),
    DOF_XX := sqrt(dofc_nss_XX / (dofc_nss_XX - 1))]
  d[is.na(dofc_s_XX) & is.na(dofc_ns_XX) & is.na(dofc_nss_XX),
    DOF_XX := NA_real_]
  invisible(d)
}


# The per-group influence-function contribution for one (event-time,
# direction). Returns one value per group, in group_XX order -- the same
# order `.core_one_event_time` returns U_g in.
#' @keywords internal
#' @noRd
.se_u_g_var <- function(d, k, G, N_inc, cohort_cols, cluster_col = NULL,
                         cols = .se_cols()) {
  if (!is.finite(N_inc) || N_inc == 0) return(numeric(G))
  .se_build_cells(d, cohort_cols, cluster_col, cols)
  .se_build_ehat(d, k)
  .se_build_dof(d, k)
  k <- as.integer(k)
  dy <- cols$diff_y; dk <- cols$dist; nc <- cols$never; rt <- cols$ratio
  d[, kern_var_XX := (G / N_inc) *
       as.integer(time_XX >= (k + 1L) & time_XX <= T_g_XX) *
       N_gt_XX *
       (get(dk) - get(rt) * get(nc)) *
       DOF_XX * (get(dy) - E_hat_XX)]
  d[is.na(kern_var_XX), kern_var_XX := 0]
  out <- d[, list(v = sum(kern_var_XX)), by = group_XX]$v
  d[, c("dof_ns_XX", "dof_s_XX", "dof_ns_s_XX", "diff_y_N_XX",
        "cnt_ns_XX", "tot_ns_XX", "mean_ns_XX", "dofc_ns_XX",
        "cnt_s_XX", "tot_s_XX", "mean_s_XX", "dofc_s_XX",
        "cnt_nss_XX", "tot_nss_XX", "mean_nss_XX", "dofc_nss_XX",
        "E_hat_XX", "DOF_XX", "kern_var_XX") := NULL]
  out
}


# se from a per-group (or per-cluster) influence vector.
#   unclustered: sqrt( sum_g U_g^2 / G^2 )
#   clustered  : sum U within cluster first, then the same.
# Reference: did_multiplegt_main.R:1032-1047.
#' @keywords internal
#' @noRd
.se_from_u <- function(u, G, cluster_of_group = NULL) {
  if (is.null(u) || !any(is.finite(u))) return(NA_real_)
  v <- .se_collapse_u(u, cluster_of_group)
  sqrt(sum(v^2, na.rm = TRUE) / (G^2))
}

#' @keywords internal
#' @noRd
.se_collapse_u <- function(u, cluster_of_group = NULL) {
  if (is.null(cluster_of_group)) return(u)
  as.numeric(rowsum(u, cluster_of_group, reorder = FALSE))
}


# Covariance between two influence vectors, by polarisation. The
# reference never forms the cross-product directly; it squares the SUM
# and subtracts, so we do the same to stay bit-comparable.
# Reference: did_multiplegt_main.R:1170-1181.
#' @keywords internal
#' @noRd
.se_cov_from_u <- function(ui, uj, se_i, se_j, G, cluster_of_group = NULL) {
  if (!is.finite(se_i) || !is.finite(se_j)) return(NA_real_)
  v <- .se_collapse_u(ui + uj, cluster_of_group)
  var_sum <- sum(v^2, na.rm = TRUE) / (G^2)
  (var_sum - se_i^2 - se_j^2) / 2
}


# Joint nullity test over a block of estimates, from the analytic
# covariance. Mirrors did_multiplegt_main.R:1163-1192: a generalised
# inverse, then a chi-square on the full rank of the block.
#' @keywords internal
#' @noRd
.se_joint_test <- function(est, U, se, G, cluster_of_group = NULL) {
  n <- length(est)
  if (n < 2L || is.null(U)) return(NA_real_)
  if (!is.numeric(G) || length(G) != 1L || !is.finite(G) || G <= 0) {
    return(NA_real_)
  }
  if (!is.matrix(U) || ncol(U) != n) return(NA_real_)
  if (!all(is.finite(est)) || !all(is.finite(se))) return(NA_real_)
  V <- matrix(0, n, n)
  diag(V) <- se^2
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      cv <- .se_cov_from_u(U[, i], U[, j], se[i], se[j], G,
                            cluster_of_group)
      if (!is.finite(cv)) return(NA_real_)
      V[i, j] <- cv
      V[j, i] <- cv
    }
  }
  Vinv <- tryCatch(MASS::ginv(V), error = function(e) NULL)
  if (is.null(Vinv)) return(NA_real_)
  chi2 <- as.numeric(t(est) %*% Vinv %*% est)
  if (!is.finite(chi2)) return(NA_real_)
  stats::pchisq(chi2, df = n, lower.tail = FALSE)
}
