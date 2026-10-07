# ============================================================================
# DIDmultiplegtDYN 2.4.0's `reset` and `avg_time_periods` options.
#
# Ports of did_multiplegt_reset() and did_multiplegt_avg_cumul()
# (MIT-licensed, Copyright (c) 2024 Diego Ciccia, Felix Knau, Melitine
# Malezieux, Doulo Sow, Clement de Chaisemartin). The two C++ routines
# behind avg_time_periods are in src/avg_cumul.cpp.
# ============================================================================


# reset = k: keep the groups with the most non-missing treatment
# observations, then start a new sub-group whenever a group's treatment
# has stayed unchanged for k periods after a change. Without a cluster
# variable the original group becomes the cluster (old_group_XX).
#' @keywords internal
#' @noRd
.dcdh_reset <- function(df, group, time, treatment, cluster, reset) {
  dt <- data.table::as.data.table(df)
  dt[, count_d_non_miss_XX := sum(!is.na(get(treatment))), by = group]
  max_count <- max(dt$count_d_non_miss_XX, na.rm = TRUE)
  dt <- dt[count_d_non_miss_XX >= max_count]
  dt[, count_d_non_miss_XX := NULL]
  data.table::setorderv(dt, c(group, time))
  dt[, changed_XX := as.integer(get(treatment) != data.table::shift(get(treatment))), by = group]
  dt[is.na(changed_XX), changed_XX := 0L]
  dt[, spell_XX := cumsum(changed_XX), by = group]
  dt[, periods_since_change_XX := seq_len(.N) - 1L, by = c(group, "spell_XX")]
  dt[spell_XX == 0L, periods_since_change_XX := 0L]
  dt[, reset_dummy_XX := as.integer(periods_since_change_XX == reset)]
  dt[, reset_count_XX := cumsum(reset_dummy_XX), by = group]
  dt[, new_group_XX := .GRP, by = c(group, "reset_count_XX")]
  if (is.null(cluster)) {
    cluster_out <- "old_group_XX"
    dt[, old_group_XX := get(group)]
  } else {
    cluster_out <- cluster
  }
  dt[[group]] <- dt$new_group_XX
  dt[, c("changed_XX", "spell_XX", "periods_since_change_XX",
         "reset_dummy_XX", "reset_count_XX", "new_group_XX") := NULL]
  list(df = as.data.frame(dt), cluster = cluster_out)
}


# avg_time_periods on the prepared panel (did_multiplegt_avg_cumul).
#' @keywords internal
#' @noRd
.dcdh_avg_cumul <- function(prepped, l_eff, same_switchers = FALSE,
                            switchers = "", continuous = NULL) {
  d <- data.table::copy(prepped[, intersect(c("group_XX", "time_XX", "treatment_XX",
                                               "outcome_XX", "d_sq_XX", "d_sq_XX_orig",
                                               "F_g_XX", "T_g_XX", "S_g_XX", "N_gt_XX"),
                                             names(prepped)), with = FALSE])
  data.table::setorderv(d, c("group_XX", "time_XX"))
  d1 <- if (!is.null(continuous) && continuous > 0 && "d_sq_XX_orig" %in% names(d))
          d$d_sq_XX_orig else d$d_sq_XX
  T_max <- max(d$time_XX)
  Fv <- as.numeric(d$F_g_XX)
  EV <- as.integer(!is.na(Fv) & Fv < (T_max + 1L))
  NGT <- as.integer(!is.na(d$outcome_XX))
  CLS <- as.integer(factor(d$d_sq_XX))
  Tg_ph <- didgpu_compute_Tg_cpp(as.integer(d$group_XX), as.integer(d$time_XX),
                                 as.numeric(d$outcome_XX), Fv,
                                 as.numeric(d$T_g_XX), EV, CLS, NGT)
  Mg_ph <- pmin(as.numeric(l_eff), Tg_ph - Fv + 1)
  didgpu_avg_cumul_cpp(as.integer(d$group_XX), as.integer(d$time_XX),
                       as.numeric(d$treatment_XX), as.numeric(d$outcome_XX),
                       as.numeric(d1), Fv, Tg_ph, EV, CLS, Mg_ph, NGT,
                       as.numeric(d$S_g_XX), as.numeric(d$N_gt_XX),
                       as.integer(l_eff), if (isTRUE(same_switchers)) 1L else 0L,
                       switchers %||% "")
}
