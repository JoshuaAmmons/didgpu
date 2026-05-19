# ============================================================================
# Pure-R reimplementation of the binary non-absorbing core (Phase 1 / "r" backend).
#
# STATUS: partial. As of this commit:
#   [x] Panel preparation: balancing, F_g, S_g, L_g, never-switcher tagging
#   [x] Single-direction switchers ("in" only)
#   [x] effects = 1, placebo = 0, no controls
#   [-] effects > 1: stub (TODO)
#   [-] placebos: stub (TODO)
#   [-] switcher-out direction: stub (TODO)
#   [-] controls / continuous / normalized / trends_lin: NOT planned for "r" backend.
#       Use backend = "reference" if you need those.
#
# Algorithm reference: vignettes/reference_internals.md sections 2-3.
# The agent-derived spec there is the source of truth; this file is the
# implementation. Any divergence from the reference output for the
# supported subset is a bug in this file.
# ============================================================================


#' Prepare a raw panel for the core estimator
#'
#' Reproduces the binary, non-absorbing-relevant subset of the
#' preparation done by `did_multiplegt_main` (sections 2 of
#' reference_internals.md). Output columns mirror the reference's
#' `_XX` suffix convention so the algorithmic code below can be
#' compared to the reference line-by-line.
#'
#' @param df A panel data.frame.
#' @param outcome,group,time,treatment Column names.
#'
#' @return A data.table with the prep columns (`group_XX`, `time_XX`,
#'   `outcome_XX`, `treatment_XX`, `d_sq_XX`, `F_g_XX`, `S_g_XX`,
#'   `T_g_XX`, `L_g_XX`, `never_change_d_XX`, `N_gt_XX`, `diff_y_XX`).
#'   The panel is balanced (one row per (group, time) over the full
#'   time range); originally missing cells are present with NA
#'   `outcome_XX`/`treatment_XX` and `N_gt_XX = 0`.
#'
#' @keywords internal
#' @noRd
.prep_panel <- function(df, outcome, group, time, treatment,
                         controls = NULL, weight = NULL,
                         trends_nonparam = NULL,
                         dont_drop_larger_lower = FALSE,
                         continuous = NULL,
                         trends_lin = FALSE) {
  d <- data.table::as.data.table(df)
  if (!is.null(controls)) {
    missing_ctrl <- setdiff(controls, names(d))
    if (length(missing_ctrl)) {
      stop("controls not found in df: ", paste(missing_ctrl, collapse = ", "))
    }
  }
  if (!is.null(weight) && !weight %in% names(d)) {
    stop("weight column not found in df: ", weight)
  }
  if (!is.null(trends_nonparam) && !trends_nonparam %in% names(d)) {
    stop("trends_nonparam column not found in df: ", trends_nonparam)
  }
  # Stash the weight column under a stable name so the rest of prep can
  # find it after we rename outcome/group/time/treatment.
  if (!is.null(weight)) {
    d[, "weight_XX_input" := as.numeric(get(weight))]
  } else {
    d[, "weight_XX_input" := 1.0]
  }
  # Stash trends_nonparam under a stable name (keeps original too).
  if (!is.null(trends_nonparam)) {
    d[, "trends_np_XX" := get(trends_nonparam)]
  }
  data.table::setnames(d, c(outcome, group, time, treatment),
                       c("outcome_XX", "group_XX", "time_XX", "treatment_XX"),
                       skip_absent = FALSE)
  # Recode group and time to consecutive integers so the U-statistic
  # arithmetic doesn't depend on raw labels.
  d[, group_XX := as.integer(factor(group_XX, levels = sort(unique(group_XX))))]
  d[, time_XX  := as.integer(factor(time_XX,  levels = sort(unique(time_XX))))]
  data.table::setkeyv(d, c("group_XX", "time_XX"))

  # Balance the panel: one row per (group, time) over the full
  # observed (group x time) grid.
  full_grid <- data.table::CJ(group_XX = sort(unique(d$group_XX)),
                              time_XX  = sort(unique(d$time_XX)))
  d <- merge(full_grid, d, by = c("group_XX", "time_XX"), all.x = TRUE)
  data.table::setkeyv(d, c("group_XX", "time_XX"))

  # Effective weight: 0 if outcome or treatment missing (excluded from
  # sums); otherwise the user's weight column (default 1).
  d[, N_gt_XX := ifelse(!is.na(outcome_XX) & !is.na(treatment_XX),
                         weight_XX_input, 0)]
  d[is.na(N_gt_XX), N_gt_XX := 0]
  d[, "weight_XX_input" := NULL]

  # dont_drop_larger_lower: by default, drop rows from groups that have
  # had BOTH a strict increase AND a strict decrease in treatment
  # (non-monotone path). Reference: main.R:185-191. The drop is partial
  # — only the rows AFTER the second-direction switch are removed, not
  # the entire group. Skip this whole block when dont_drop_larger_lower
  # is TRUE.
  if (!isTRUE(dont_drop_larger_lower)) {
    t_min0 <- min(d$time_XX)
    # Per-group baseline (first-period) treatment, broadcast in place.
    d[, d_sq_tmp := treatment_XX[time_XX == t_min0][1L], by = group_XX]
    d[, diff_from_sq_tmp := treatment_XX - d_sq_tmp]
    data.table::setorder(d, group_XX, time_XX)
    d[, ever_strict_increase_tmp := as.integer(pmin(1L,
        cumsum(diff_from_sq_tmp > 0 & !is.na(treatment_XX)))),
      by = group_XX]
    d[, ever_strict_decrease_tmp := as.integer(pmin(1L,
        cumsum(diff_from_sq_tmp < 0 & !is.na(treatment_XX)))),
      by = group_XX]
    d <- d[!(ever_strict_increase_tmp == 1L & ever_strict_decrease_tmp == 1L)]
    d[, c("d_sq_tmp", "diff_from_sq_tmp",
          "ever_strict_increase_tmp", "ever_strict_decrease_tmp") := NULL]
  }

  # Baseline (period-1) treatment per group. NA if the group is unobserved
  # at the first period of the balanced grid. Computed in place with
  # grouped assignment (faster than merge: no second allocation + bmerge).
  t_min <- min(d$time_XX)
  d[, d_sq_XX := treatment_XX[time_XX == t_min][1L], by = group_XX]

  # First-switch period F_g_XX: the smallest t with treatment_XX != d_sq_XX
  # (and both observed). For never-switchers, F_g_XX = T_max + 1.
  T_max <- max(d$time_XX)
  d[, F_g_XX := {
      ok <- !is.na(treatment_XX) & !is.na(d_sq_XX) & treatment_XX != d_sq_XX
      if (any(ok)) as.integer(min(time_XX[ok])) else NA_integer_
    },
    by = group_XX]
  d[is.na(F_g_XX), F_g_XX := as.integer(T_max + 1L)]

  # T_g_XX: last usable time per group (= last time with N_gt > 0).
  d[, T_g_XX := {
      ok <- N_gt_XX > 0
      if (any(ok)) max(time_XX[ok]) else NA_integer_
    },
    by = group_XX]

  # Switcher direction S_g_XX based on average post-switch treatment.
  # 1 = switcher-in (avg post > baseline). 0 = switcher-out (avg post < baseline).
  # NA = never-switcher or unobserved post-switch.
  d[, avg_post := {
      ok <- time_XX >= F_g_XX & N_gt_XX > 0
      if (any(ok)) mean(treatment_XX[ok]) else NA_real_
    },
    by = group_XX]
  d[, S_g_XX := ifelse(is.na(avg_post) | !is.finite(F_g_XX), NA_integer_,
                ifelse(avg_post > d_sq_XX, 1L,
                ifelse(avg_post < d_sq_XX, 0L, NA_integer_)))]

  # Continuous treatment. Reference: main.R:202-208 + 306-318. Three
  # things happen when `continuous` is set:
  #
  #  (1) Save originals, collapse d_sq cohorts (d_sq_XX = 0 for all).
  #      F_g, S_g, L_g were already computed on the ORIGINAL d_sq, which
  #      is correct (they mean "first period D differed from baseline").
  #
  #  (2) BINARIZE treatment_XX to a +1/-1 signed indicator:
  #      treatment_XX := (F_g <= time) * S_g_het, where S_g_het = +1 for
  #      in-switchers and -1 for out-switchers. This is what flows into
  #      the kernel's diff_y / dist calculations. The original treatment
  #      values are kept as treatment_XX_orig.
  #
  #  (3) Add (time-onset >= j) * baseline^k interaction controls for
  #      j in 2..T_max and k in 1..degree_pol. These absorb time-varying
  #      effects of the continuous baseline.
  if (!is.null(continuous)) {
    degree_pol <- as.integer(continuous)
    stopifnot(degree_pol >= 1L)
    d[, d_sq_XX_orig := d_sq_XX]
    for (pol in seq_len(degree_pol)) {
      d[, paste0("d_sq_", pol, "_XX") := d_sq_XX_orig^pol]
    }
    d[, d_sq_XX := 0]
    # Binarize treatment to a +1/-1 indicator (S_g_het = +1 in, -1 out).
    d[, S_g_het_XX := ifelse(is.na(S_g_XX), NA_integer_,
                              ifelse(S_g_XX == 0L, -1L, S_g_XX))]
    d[, treatment_XX_orig := treatment_XX]
    d[, treatment_XX := ifelse(is.na(S_g_het_XX), NA_real_,
                                as.numeric((F_g_XX <= time_XX) * S_g_het_XX))]
    # Add time x baseline^k interaction controls.
    T_max_full <- max(d$time_XX)
    auto_ctrl_names <- character(0)
    for (j in 2L:T_max_full) {
      for (k in seq_len(degree_pol)) {
        col <- sprintf("time_fe_XX_%d_bt%d_XX", j, k)
        d[, (col) := as.numeric(time_XX >= j) *
            get(paste0("d_sq_", k, "_XX"))]
        auto_ctrl_names <- c(auto_ctrl_names, col)
      }
    }
    data.table::setattr(d, "didgpu_auto_controls", auto_ctrl_names)
  }

  # L_g_XX = horizon of post-switch periods available for this group.
  d[, L_g_XX := pmax(0L, as.integer(T_g_XX - F_g_XX + 1L))]
  d[is.na(L_g_XX), L_g_XX := 0L]

  # L_g_placebo_XX = max placebo horizon supported for this group.
  # Reference did_multiplegt_main.R:324: min(L_g, F_g - 2) when F_g >= 3.
  # This is the largest k such that BOTH t = F_g + k - 1 fits in the
  # observed window AND t - 2k = F_g - k - 1 >= 1 (placebo diff exists).
  d[, L_g_placebo_XX := ifelse(F_g_XX >= 3L,
                                pmin(L_g_XX, as.integer(F_g_XX - 2L)),
                                NA_integer_)]

  # never_change_d_XX: at this (g, t), the group has not yet switched
  # (i.e., is still a valid control for any switcher reaching event-time
  # k such that t < F_g_other for that other group).
  d[, never_change_d_XX := as.integer(time_XX < F_g_XX & N_gt_XX > 0)]

  # trends_lin = TRUE: first-difference outcome (and any user controls)
  # per group, drop the t = t_min rows (FD is NA there), and drop
  # F_g == 2 cohorts (no pre-period FD exists for them). Reference:
  # did_multiplegt_main.R:289-298. After this transformation, the
  # estimator runs on Y_fd and the final cumulative-recovery step lives
  # in the orchestrator (.compute_effects_trends_lin).
  if (isTRUE(trends_lin)) {
    t_min_fd <- min(d$time_XX)
    # Drop F_g == t_min + 1 (i.e. switchers at the first usable period
    # after FD); they would have no pre-treatment FD observation.
    d <- d[F_g_XX != (t_min_fd + 1L)]
    data.table::setorder(d, group_XX, time_XX)
    # First-difference outcome and any user-supplied controls per group.
    fd_cols <- c("outcome_XX", controls)
    for (v in fd_cols) {
      if (!is.null(v) && v %in% names(d)) {
        d[, (v) := get(v) - data.table::shift(get(v), 1L, type = "lag"),
          by = group_XX]
      }
    }
    # Drop the t = t_min rows (FD is NA there).
    d <- d[time_XX != t_min_fd]
    # Re-derive N_gt_XX on the post-FD data — rows that are now NA in
    # outcome_XX or treatment_XX should drop out of all sums.
    d[, N_gt_XX := ifelse(!is.na(outcome_XX) & !is.na(treatment_XX),
                           N_gt_XX, 0)]
    # Tag prepped panel so .compute_effects knows to use trends_lin logic.
    data.table::setattr(d, "didgpu_trends_lin", TRUE)
  }

  # First difference of outcome (used for k=1; longer differences built
  # per-event-time inside the core).
  d[, diff_y_XX := outcome_XX - data.table::shift(outcome_XX, 1L,
                                                   type = "lag"),
    by = group_XX]

  d
}


#' Compute the per-event-time U-statistic for a single switcher direction
#'
#' Implements the core formula in section 3 of reference_internals.md
#' for event-time `k` only, single direction, no controls, no placebos.
#' This is the kernel that would later be moved to C++ (cpu backend)
#' and CUDA (cuda backend).
#'
#' @param d The prepared panel (output of .prep_panel).
#' @param k Integer event-time (1, 2, ...).
#' @param direction 1 for switcher-in, 0 for switcher-out.
#' @return A list with `att` (the event-time-k ATT estimate),
#'   `N_inc` (number of contributing switcher cells), and
#'   `U_g` (per-group contribution to the numerator).
#'
#' @keywords internal
#' @noRd
.core_one_event_time <- function(d_in, k, direction = 1L, prefit = NULL,
                                  only_never_switchers = FALSE,
                                  normalized = FALSE,
                                  skip_prep = FALSE) {
  k <- as.integer(k)
  stopifnot(k >= 1L)
  # We operate directly on `d_in` (data.table by reference). Per-call
  # scratch columns (diff_y_k_XX, dist_k_XX, never_change_k_XX,
  # N_t_control, N_t_switch, ratio_XX, kernel_XX, contrib_XX) are
  # overwritten on every call with the SAME column names, so successive
  # calls (different k or direction) don't see each other's residue.
  # Copying defensively was costing ~12% of total time at scale.
  #
  # The chunk above the `if (skip_prep)` gate is direction-independent —
  # only diff_y_k_XX, controls adjustment, never_change_k_XX, and
  # N_t_control. .compute_effects sets `skip_prep = TRUE` on the second
  # of the in/out direction calls at the same `k` to avoid redoing
  # this work, which saves ~20% wall time at large panels.
  d <- d_in
  G <- length(unique(d$group_XX))
  T_max_XX <- max(d$time_XX)

  # Cohort grouping columns: (time, d_sq) plus trends_nonparam if present.
  cohort_cols <- c("time_XX", "d_sq_XX")
  if ("trends_np_XX" %in% names(d)) cohort_cols <- c(cohort_cols, "trends_np_XX")

  if (!isTRUE(skip_prep)) {
    # Long difference at horizon k.
    if (k == 1L) {
      d[, diff_y_k_XX := diff_y_XX]
    } else {
      d[, diff_y_k_XX := outcome_XX -
          data.table::shift(outcome_XX, k, type = "lag"),
        by = group_XX]
    }

    # Apply controls adjustment: for each baseline level l, subtract
    # sum_c theta_d[l, c] * (X_c_t - lag(X_c_t, k)) from diff_y_k for
    # rows with d_sq == l. Reference: did_multiplegt_dyn_core.R:267-269.
    if (!is.null(prefit) && length(prefit$controls) > 0L) {
      for (l in prefit$dsq_levels) {
        key <- as.character(l)
        theta_l <- prefit$theta[[key]]
        if (is.null(theta_l) || all(theta_l == 0)) next
        mask <- d$d_sq_XX == l
        adj <- rep(0, nrow(d))
        for (c in seq_along(prefit$controls)) {
          cname <- prefit$controls[c]
          # Long-difference the c-th control at horizon k.
          d[, ".__diff_X_k_XX" := get(cname) -
              data.table::shift(get(cname), k, type = "lag"),
            by = group_XX]
          dx <- d$.__diff_X_k_XX
          adj <- adj + theta_l[c] * ifelse(is.na(dx), 0, dx)
        }
        d[mask, diff_y_k_XX := diff_y_k_XX - adj[mask]]
        d[, ".__diff_X_k_XX" := NULL]
      }
    }

    # Controls first (needed to gate the switcher mask).
    d[, never_change_k_XX := as.integer(time_XX < F_g_XX &
                                         N_gt_XX > 0 &
                                         !is.na(diff_y_k_XX))]
    d[is.na(never_change_k_XX), never_change_k_XX := 0L]

    # only_never_switchers: restrict controls to TRUE never-switchers
    # (drop pre-switch rows of units that eventually switch). Reference:
    # did_multiplegt_dyn_core.R:132-134.
    if (isTRUE(only_never_switchers)) {
      d[F_g_XX < T_max_XX + 1L, never_change_k_XX := 0L]
    }

    # Per-(time, baseline-cohort) control mass. Reference:
    # did_multiplegt_dyn_core.R:137. Grouped by (time, d_sq, trends_nonparam).
    # Using in-place grouped assignment instead of merge — same result, ~5x
    # faster (no second data.table allocation + bmerge).
    d[, N_t_control := sum(N_gt_XX * never_change_k_XX),
      by = cohort_cols]
  }

  # Switcher-cell indicator: at time t == F_g + k - 1, group is a
  # switcher in this direction, the horizon is reachable, a valid
  # diff_y exists at this row, AND there are controls in the same
  # (time, d_sq) cohort. The N_t_control > 0 condition mirrors the
  # reference at did_multiplegt_dyn_core.R:197 and is essential when
  # a switcher's cohort has no concurrent controls.
  d[, dist_k_XX := as.integer(
        time_XX == (F_g_XX + k - 1L) &
        k <= L_g_XX &
        !is.na(S_g_XX) & S_g_XX == direction &
        !is.na(diff_y_k_XX) &
        N_gt_XX > 0 &
        !is.na(N_t_control) & N_t_control > 0
      )]
  d[is.na(dist_k_XX), dist_k_XX := 0L]
  # If same_switchers is in force, also gate on still_switcher_XX == 1
  # (the unit must qualify at EVERY event-time q in 1..effects).
  if ("still_switcher_XX" %in% names(d)) {
    d[still_switcher_XX != 1L, dist_k_XX := 0L]
  }

  # Per-(time, cohort) switcher mass (used in the ratio N_t_g / N_t_control).
  d[, N_t_switch := sum(N_gt_XX * dist_k_XX),
    by = cohort_cols]

  # Total incidence: total switcher contribution across all valid times.
  N_inc <- sum(d$N_gt_XX * d$dist_k_XX, na.rm = TRUE)

  if (N_inc == 0) {
    # No switchers reach event-time k; ATT is undefined.
    return(list(att = NA_real_, N_inc = 0L, U_g = numeric(G)))
  }

  # The U-statistic kernel. From section 3 of reference_internals.md:
  #   U_gt_temp = dummy * (G/N_inc) * 1[t in (k+1)..T_g] * N_gt
  #               * (dist - (N_t_switch / N_t_control) * never_change_k)
  #               * diff_y_k
  # For our simplified case (no `same_switchers` filter, no controls),
  # `dummy` is just `S_g_XX == direction` (i.e., the group is a
  # switcher in this direction at all). For controls and switchers
  # alike the indicator is on; the per-row contribution is non-zero
  # only when at least one of `dist_k_XX` or `never_change_k_XX` is 1.
  d[, ratio_XX := ifelse(N_t_control > 0,
                          N_t_switch / N_t_control, 0)]
  d[, kernel_XX := (G / N_inc) *
                   as.integer(time_XX >= (k + 1L) & time_XX <= T_g_XX) *
                   N_gt_XX *
                   (dist_k_XX - ratio_XX * never_change_k_XX) *
                   diff_y_k_XX]
  d[is.na(kernel_XX), kernel_XX := 0]

  U_g <- d[, list(U_g = sum(kernel_XX)), by = group_XX]

  # Per the reference (did_multiplegt_main.R:892-893):
  #   DID_i = sum_g U_Gg_i / G
  # i.e., average per-group U-statistic across ALL groups (including
  # controls and never-switchers; their U_g entries are 0).
  att <- sum(U_g$U_g) / G

  # N_eff = count of contributing observations (rows). Mirrors the
  # reference's count_i_core_XX summation. A row contributes if it
  # is a switcher cell OR an active control cell (never_change in a
  # cohort that has at least one switcher at this time and within the
  # post-window). For binary panels with no weights, N_gt is always
  # 0 or 1.
  d[, contrib_XX := as.integer(
        as.integer(time_XX >= (k + 1L) & time_XX <= T_g_XX) *
        N_gt_XX *
        ((dist_k_XX == 1L) |
         (never_change_k_XX == 1L &
          !is.na(N_t_switch) & N_t_switch > 0))
      )]
  d[is.na(contrib_XX), contrib_XX := 0L]
  N_eff <- sum(d$contrib_XX, na.rm = TRUE)

  # delta_norm: average cumulative treatment-change magnitude for switchers
  # at event-time k in this direction. Reference: did_multiplegt_dyn_core.R
  # lines 772-796. Used by main.R:1086-1090 to build delta_D_i_global, which
  # divides DID_i_XX (line 1124-1126) when normalized=TRUE.
  #
  # Per-row contribution: at switcher cells (dist == 1):
  #   (N_gt / N_inc) * sign * sum_treat_until_k_g
  # where sign = +1 for in (S_g=1), -1 for out (S_g=0) — i.e. sign = 2*S_g-1.
  # sum_treat_until_k_g = sum over t in [F_g, F_g+k-1] of (D_t - D_baseline).
  # For continuous treatment, use ORIGINAL D values (treatment_XX has been
  # binarized in prep) — reference at lines 773-784.
  delta_norm <- if (isTRUE(normalized)) {
    has_orig <- all(c("treatment_XX_orig", "d_sq_XX_orig") %in% names(d))
    t_col <- if (has_orig) "treatment_XX_orig" else "treatment_XX"
    b_col <- if (has_orig) "d_sq_XX_orig"      else "d_sq_XX"
    d[, sum_temp_XX := ifelse(
          time_XX >= F_g_XX &
          time_XX <= (F_g_XX - 1L + k) &
          !is.na(S_g_XX) & S_g_XX == direction,
          get(t_col) - get(b_col), NA_real_)]
    d[, sum_treat_until_XX := sum(sum_temp_XX, na.rm = TRUE),
      by = group_XX]
    sign_dir <- if (direction == 1L) 1 else -1
    d[, delta_cum_XX := ifelse(
          dist_k_XX == 1L,
          (N_gt_XX / N_inc) * sign_dir * sum_treat_until_XX,
          NA_real_)]
    val <- sum(d$delta_cum_XX, na.rm = TRUE)
    d[, c("sum_temp_XX", "sum_treat_until_XX", "delta_cum_XX") := NULL]
    val
  } else NA_real_

  list(att = att,
       N_inc = as.integer(N_inc),
       N_eff = as.integer(N_eff),
       U_g = U_g$U_g,
       delta_norm = delta_norm)
}


# -------- same_switchers pre-pass --------

#' Compute the per-group still_switcher_XX indicator
#'
#' Mirrors did_multiplegt_dyn_core.R:139-199 (the same_switchers branch).
#' For each event-time q in 1..effects, checks that the switcher at row
#' (t = F_g + q - 1) has positive control mass in its (time, d_sq)
#' cohort AND a valid long-difference diff_y_q. A unit qualifies as a
#' still_switcher iff:
#'   (a) F_g + effects - 1 <= T_g (enough post-treatment data), AND
#'   (b) ALL `effects` event-time checks pass.
#'
#' Writes `still_switcher_XX` (integer 0/1) onto `prepped` in place.
#'
#' Compute the per-group still_switcher_pl_XX indicator (placebo-side)
#'
#' Mirrors did_multiplegt_dyn_core.R:177-215 (the same_switchers_pl
#' branch). For each placebo horizon q in 1..placebo, checks that the
#' switcher at row (t = F_g - 1 - q) has positive control mass in its
#' (time, d_sq) cohort AND a valid pre-period diff_y (computed as
#' outcome - lead(outcome, q), which equals Y_{t} - Y_{t+q}, so at
#' t = F_g - 1 - q this is Y_{F_g - 1 - q} - Y_{F_g - 1}, the negative
#' of the standard placebo diff).
#'
#' Writes `still_switcher_pl_XX` (integer 0/1) onto `prepped` in place.
#'
#' @keywords internal
#' @noRd
.compute_still_switcher_pl <- function(prepped, placebo,
                                         only_never_switchers = FALSE) {
  d <- prepped
  T_max_XX <- max(d$time_XX)
  cohort_cols <- c("time_XX", "d_sq_XX")
  if ("trends_np_XX" %in% names(d)) cohort_cols <- c(cohort_cols, "trends_np_XX")

  d[, N_g_control_check_pl_XX := 0L]
  for (q in seq_len(placebo)) {
    d[, diff_y_last_pl_XX := outcome_XX -
        data.table::shift(outcome_XX, n = q, type = "lead"),
      by = group_XX]
    d[, never_change_d_last_pl_XX := as.integer(
          !is.na(diff_y_last_pl_XX) & F_g_XX > time_XX)]
    if (isTRUE(only_never_switchers)) {
      d[F_g_XX > time_XX & F_g_XX < (T_max_XX + 1L) & !is.na(diff_y_last_pl_XX),
        never_change_d_last_pl_XX := 0L]
    }
    d[is.na(never_change_d_last_pl_XX), never_change_d_last_pl_XX := 0L]
    d[, N_gt_control_last_pl_XX := sum(never_change_d_last_pl_XX * N_gt_XX,
                                         na.rm = TRUE),
      by = cohort_cols]
    # Per-group: the value at the row where t == F_g - 1 - q (NA elsewhere).
    d[, N_g_control_last_m_pl_XX := mean(
          ifelse(time_XX == (F_g_XX - 1L - q),
                  N_gt_control_last_pl_XX, NA_real_),
          na.rm = TRUE),
      by = group_XX]
    d[, diff_y_relev_pl_XX := mean(
          ifelse(time_XX == (F_g_XX - 1L - q),
                  diff_y_last_pl_XX, NA_real_),
          na.rm = TRUE),
      by = group_XX]
    d[, N_g_control_check_pl_XX := N_g_control_check_pl_XX + as.integer(
          !is.na(N_g_control_last_m_pl_XX) & N_g_control_last_m_pl_XX > 0 &
          !is.na(diff_y_relev_pl_XX))]
  }
  d[, still_switcher_pl_XX := as.integer(
        N_g_control_check_pl_XX == placebo)]
  d[, c("N_g_control_check_pl_XX", "diff_y_last_pl_XX",
        "never_change_d_last_pl_XX", "N_gt_control_last_pl_XX",
        "N_g_control_last_m_pl_XX", "diff_y_relev_pl_XX") := NULL]
  invisible(d)
}


#' @keywords internal
#' @noRd
.compute_still_switcher <- function(prepped, effects, only_never_switchers = FALSE) {
  d <- prepped
  T_max_XX <- max(d$time_XX)
  cohort_cols <- c("time_XX", "d_sq_XX")
  if ("trends_np_XX" %in% names(d)) cohort_cols <- c(cohort_cols, "trends_np_XX")

  d[, N_g_control_check_XX := 0L]
  for (q in seq_len(effects)) {
    d[, diff_y_last_XX := outcome_XX -
        data.table::shift(outcome_XX, n = q, type = "lag"),
      by = group_XX]
    d[, never_change_d_last_XX := as.integer(
          !is.na(diff_y_last_XX) & F_g_XX > time_XX)]
    if (isTRUE(only_never_switchers)) {
      d[F_g_XX > time_XX & F_g_XX < (T_max_XX + 1L) & !is.na(diff_y_last_XX),
        never_change_d_last_XX := 0L]
    }
    d[is.na(never_change_d_last_XX), never_change_d_last_XX := 0L]
    d[, N_gt_control_last_XX := sum(never_change_d_last_XX * N_gt_XX,
                                     na.rm = TRUE),
      by = cohort_cols]
    # Per-group: the value at the row where t == F_g + q - 1 (NA elsewhere).
    d[, N_g_control_last_m_XX := mean(
          ifelse(time_XX == (F_g_XX + q - 1L), N_gt_control_last_XX, NA_real_),
          na.rm = TRUE),
      by = group_XX]
    d[, diff_y_relev_XX := mean(
          ifelse(time_XX == (F_g_XX + q - 1L), diff_y_last_XX, NA_real_),
          na.rm = TRUE),
      by = group_XX]
    # Increment the per-group check counter.
    d[, N_g_control_check_XX := N_g_control_check_XX + as.integer(
          !is.na(N_g_control_last_m_XX) & N_g_control_last_m_XX > 0 &
          !is.na(diff_y_relev_XX))]
  }
  d[, still_switcher_XX := as.integer(
        (F_g_XX + effects - 1L) <= T_g_XX &
        N_g_control_check_XX == effects)]
  # Clean up scratch.
  d[, c("N_g_control_check_XX", "diff_y_last_XX",
        "never_change_d_last_XX", "N_gt_control_last_XX",
        "N_g_control_last_m_XX", "diff_y_relev_XX") := NULL]
  invisible(d)
}


# -------- placebos (mirrors effects with diff_y_pl_k = Y_{t-2k} - Y_{t-k}) --------

#' Compute the per-event-time placebo for one switcher direction
#'
#' Placebo horizon `k` uses the same switcher mask as effect horizon
#' `k` (units at `t = F_g + k - 1`) but replaces `diff_y_k` with the
#' "pre-treatment mirror" `diff_y_pl_k = Y_{t-2k} - Y_{t-k}`. For a
#' switcher this is `Y_{F_g - k - 1} - Y_{F_g - 1}`: the difference
#' between two pre-treatment outcomes spanning the same `k`-period
#' gap that the effect difference spans on the post-treatment side.
#'
#' Reference: did_multiplegt_dyn_core.R:576-891.
#'
#' @keywords internal
#' @noRd
.core_one_placebo <- function(d_in, k, direction = 1L, prefit = NULL,
                               only_never_switchers = FALSE,
                               normalized = FALSE,
                               skip_prep = FALSE) {
  k <- as.integer(k)
  stopifnot(k >= 1L)
  # Same in-place mutation pattern as .core_one_event_time — see notes there.
  d <- d_in
  G <- length(unique(d$group_XX))
  T_max_XX <- max(d$time_XX)

  cohort_cols <- c("time_XX", "d_sq_XX")
  if ("trends_np_XX" %in% names(d)) cohort_cols <- c(cohort_cols, "trends_np_XX")

  if (!isTRUE(skip_prep)) {
    # Placebo difference: Y_{t-2k} - Y_{t-k}. Reference line 624.
    d[, diff_y_pl_k_XX := data.table::shift(outcome_XX, 2L * k, type = "lag") -
                          data.table::shift(outcome_XX, k, type = "lag"),
      by = group_XX]

    # Apply controls adjustment to the placebo difference, analogous to
    # the effects path. Reference: did_multiplegt_dyn_core.R:681-683.
    if (!is.null(prefit) && length(prefit$controls) > 0L) {
      for (l in prefit$dsq_levels) {
        key <- as.character(l)
        theta_l <- prefit$theta[[key]]
        if (is.null(theta_l) || all(theta_l == 0)) next
        mask <- d$d_sq_XX == l
        adj <- rep(0, nrow(d))
        for (c in seq_along(prefit$controls)) {
          cname <- prefit$controls[c]
          d[, ".__diff_X_pl_k_XX" :=
              data.table::shift(get(cname), 2L * k, type = "lag") -
              data.table::shift(get(cname), k, type = "lag"),
            by = group_XX]
          dx <- d$.__diff_X_pl_k_XX
          adj <- adj + theta_l[c] * ifelse(is.na(dx), 0, dx)
        }
        d[mask, diff_y_pl_k_XX := diff_y_pl_k_XX - adj[mask]]
        d[, ".__diff_X_pl_k_XX" := NULL]
      }
    }

    # Controls first (needed to gate the switcher mask).
    d[, never_change_k_pl_XX := as.integer(time_XX < F_g_XX &
                                            N_gt_XX > 0 &
                                            !is.na(diff_y_pl_k_XX))]
    d[is.na(never_change_k_pl_XX), never_change_k_pl_XX := 0L]
    if (isTRUE(only_never_switchers)) {
      d[F_g_XX < T_max_XX + 1L, never_change_k_pl_XX := 0L]
    }

    # Per-cohort control mass (in-place; see core_one_event_time).
    d[, N_t_control_pl := sum(N_gt_XX * never_change_k_pl_XX),
      by = cohort_cols]
  }

  # Switcher-cell indicator: same time/direction/horizon mask as the
  # effect at k, plus a valid placebo diff and concurrent controls in
  # the same (time, d_sq) cohort.
  d[, dist_k_pl_XX := as.integer(
        time_XX == (F_g_XX + k - 1L) &
        k <= L_g_XX &
        !is.na(S_g_XX) & S_g_XX == direction &
        !is.na(diff_y_pl_k_XX) &
        N_gt_XX > 0 &
        !is.na(N_t_control_pl) & N_t_control_pl > 0
      )]
  d[is.na(dist_k_pl_XX), dist_k_pl_XX := 0L]
  # same_switchers_pl gate: only switchers that qualify at EVERY placebo
  # horizon q in 1..placebo contribute. Reference: did_multiplegt_dyn_core.R:177-215.
  if ("still_switcher_pl_XX" %in% names(d)) {
    d[still_switcher_pl_XX != 1L, dist_k_pl_XX := 0L]
  }

  d[, N_t_switch_pl := sum(N_gt_XX * dist_k_pl_XX),
    by = cohort_cols]

  N_inc <- sum(d$N_gt_XX * d$dist_k_pl_XX, na.rm = TRUE)
  if (N_inc == 0) {
    return(list(att = NA_real_, N_inc = 0L, U_g = numeric(G)))
  }

  d[, ratio_pl_XX := ifelse(N_t_control_pl > 0,
                             N_t_switch_pl / N_t_control_pl, 0)]
  d[, kernel_pl_XX := (G / N_inc) *
                      as.integer(time_XX >= (k + 1L) & time_XX <= T_g_XX) *
                      N_gt_XX *
                      (dist_k_pl_XX - ratio_pl_XX * never_change_k_pl_XX) *
                      diff_y_pl_k_XX]
  d[is.na(kernel_pl_XX), kernel_pl_XX := 0]

  U_g <- d[, list(U_g = sum(kernel_pl_XX)), by = group_XX]
  att <- sum(U_g$U_g) / G

  d[, contrib_pl_XX := as.integer(
        as.integer(time_XX >= (k + 1L) & time_XX <= T_g_XX) *
        N_gt_XX *
        ((dist_k_pl_XX == 1L) |
         (never_change_k_pl_XX == 1L &
          !is.na(N_t_switch_pl) & N_t_switch_pl > 0))
      )]
  d[is.na(contrib_pl_XX), contrib_pl_XX := 0L]
  N_eff <- sum(d$contrib_pl_XX, na.rm = TRUE)

  # delta_norm for placebos. Same formula as effects (reference: did_multiplegt_dyn_core.R
  # lines 1234-1260; main.R:1238-1242 does the corresponding division).
  delta_norm <- if (isTRUE(normalized)) {
    has_orig <- all(c("treatment_XX_orig", "d_sq_XX_orig") %in% names(d))
    t_col <- if (has_orig) "treatment_XX_orig" else "treatment_XX"
    b_col <- if (has_orig) "d_sq_XX_orig"      else "d_sq_XX"
    d[, sum_temp_pl_XX := ifelse(
          time_XX >= F_g_XX &
          time_XX <= (F_g_XX - 1L + k) &
          !is.na(S_g_XX) & S_g_XX == direction,
          get(t_col) - get(b_col), NA_real_)]
    d[, sum_treat_until_pl_XX := sum(sum_temp_pl_XX, na.rm = TRUE),
      by = group_XX]
    sign_dir <- if (direction == 1L) 1 else -1
    d[, delta_cum_pl_XX := ifelse(
          dist_k_pl_XX == 1L,
          (N_gt_XX / N_inc) * sign_dir * sum_treat_until_pl_XX,
          NA_real_)]
    val <- sum(d$delta_cum_pl_XX, na.rm = TRUE)
    d[, c("sum_temp_pl_XX", "sum_treat_until_pl_XX", "delta_cum_pl_XX") := NULL]
    val
  } else NA_real_

  list(att = att,
       N_inc = as.integer(N_inc),
       N_eff = as.integer(N_eff),
       U_g = U_g$U_g,
       delta_norm = delta_norm)
}


#' Compute per-event-time placebos for k = 1..placebo, both directions
#' @keywords internal
#' @noRd
.compute_placebos <- function(prepped, placebo, switchers = "", prefit = NULL,
                               only_never_switchers = FALSE,
                               normalized = FALSE,
                               same_switchers_pl = FALSE) {
  if (placebo == 0L) return(list(placebos = numeric(0), n_inc = integer(0),
                                  n_eff = integer(0),
                                  delta_D = numeric(0)))
  # same_switchers_pl: only switchers who have valid pre-period diff_y at
  # every placebo horizon q in 1..placebo contribute. Reference:
  # did_multiplegt_dyn_core.R:177-215.
  if (isTRUE(same_switchers_pl)) {
    .compute_still_switcher_pl(prepped, placebo, only_never_switchers)
  } else if ("still_switcher_pl_XX" %in% names(prepped)) {
    prepped[, still_switcher_pl_XX := NULL]
  }
  out <- numeric(placebo)
  n_inc <- integer(placebo)
  n_eff <- integer(placebo)
  delta_D <- numeric(placebo)
  both_dirs <- (switchers == "")
  for (k in seq_len(placebo)) {
    res_in  <- if (switchers != "out") .core_one_placebo(prepped, k = k, direction = 1L, prefit = prefit, only_never_switchers = only_never_switchers, normalized = normalized)
               else list(att = NA_real_, N_inc = 0L, N_eff = 0L, delta_norm = NA_real_)
    res_out <- if (switchers != "in")  .core_one_placebo(prepped, k = k, direction = 0L, prefit = prefit, only_never_switchers = only_never_switchers, normalized = normalized, skip_prep = both_dirs)
               else list(att = NA_real_, N_inc = 0L, N_eff = 0L, delta_norm = NA_real_)
    n_in  <- res_in$N_inc
    n_out <- res_out$N_inc
    if (n_in + n_out == 0L) {
      out[k]   <- NA_real_
      n_inc[k] <- 0L
      n_eff[k] <- 0L
      delta_D[k] <- NA_real_
      next
    }
    # Sign-flip convention for the out direction (see .compute_effects).
    att_in       <- if (n_in  > 0L)  res_in$att  else 0
    att_out_pool <- if (n_out > 0L) -res_out$att else 0
    w_in <- n_in / (n_in + n_out)
    out[k]   <- w_in * att_in + (1 - w_in) * att_out_pool
    n_inc[k] <- n_in + n_out
    n_eff[k] <- (res_in$N_eff %||% 0L) + (res_out$N_eff %||% 0L)
    if (isTRUE(normalized)) {
      dn_in  <- if (n_in  > 0L) res_in$delta_norm  else NA_real_
      dn_out <- if (n_out > 0L) res_out$delta_norm else NA_real_
      delta_k <- if (n_in == 0L) dn_out
                 else if (n_out == 0L) dn_in
                 else w_in * dn_in + (1 - w_in) * dn_out
      delta_D[k] <- delta_k
      if (!is.na(delta_k) && delta_k != 0) {
        out[k] <- out[k] / delta_k
      } else {
        out[k] <- NA_real_
      }
    }
  }
  list(placebos = out, n_inc = n_inc, n_eff = n_eff, delta_D = delta_D)
}


# -------- multi-event-time orchestration --------

#' Compute per-event-time effects for k = 1..effects, both directions
#'
#' Returns a numeric vector of length `effects` where element i is the
#' per-event-time DID_i = sum_g U_Gg_i_global / G,
#' with U_Gg_i_global = w_in * U_Gg_i_in + (1 - w_in) * U_Gg_i_out
#' and w_in the Neyman weight (relative incidence count for direction "in").
#'
#' @keywords internal
#' @noRd
.compute_effects <- function(prepped, effects, switchers = "", prefit = NULL,
                              only_never_switchers = FALSE,
                              same_switchers = FALSE,
                              normalized = FALSE) {
  out <- numeric(effects)
  n_inc <- integer(effects)
  n_eff <- integer(effects)
  delta_D <- numeric(effects)
  # same_switchers gates dist on a per-group "still_switcher" indicator
  # that requires the switcher to qualify at EVERY event-time q in
  # 1..effects (reference: did_multiplegt_dyn_core.R:139-199).
  # Compute it once before the per-k loop and pass via prepped column.
  if (isTRUE(same_switchers)) {
    .compute_still_switcher(prepped, effects, only_never_switchers)
  } else if ("still_switcher_XX" %in% names(prepped)) {
    prepped[, still_switcher_XX := NULL]
  }
  out_raw <- numeric(effects)
  # When both directions are requested at the same k, the prep work
  # (diff_y_k, never_change_k, N_t_control) is direction-independent.
  # We compute it during the first call and pass skip_prep=TRUE on the
  # second to avoid duplicating it. Saves ~15-20% wall time on
  # both-directions runs.
  both_dirs <- (switchers == "")
  for (k in seq_len(effects)) {
    res_in  <- if (switchers != "out") .core_one_event_time(prepped, k = k, direction = 1L, prefit = prefit, only_never_switchers = only_never_switchers, normalized = normalized)
               else list(att = NA_real_, N_inc = 0L, N_eff = 0L, delta_norm = NA_real_)
    res_out <- if (switchers != "in")  .core_one_event_time(prepped, k = k, direction = 0L, prefit = prefit, only_never_switchers = only_never_switchers, normalized = normalized, skip_prep = both_dirs)
               else list(att = NA_real_, N_inc = 0L, N_eff = 0L, delta_norm = NA_real_)
    # Neyman pooling across directions, exactly as the reference does at
    # did_multiplegt_main.R:859-861. SIGN CONVENTION: at line 804 of
    # main, the reference negates the out-direction U_Gg
    # (`U_Gg_i_minus_XX = -U_Gg_i_XX`) so that a "treatment removal"
    # whose Y drops by 0.5 is reported as a +0.5 treatment effect.
    # Hence att_out_for_pooling = -att_out.
    n_in  <- res_in$N_inc
    n_out <- res_out$N_inc
    if (n_in + n_out == 0L) {
      out[k]      <- NA_real_
      out_raw[k]  <- NA_real_
      n_inc[k]    <- 0L
      n_eff[k]    <- 0L
      delta_D[k]  <- NA_real_
      next
    }
    att_in       <- if (n_in  > 0L)  res_in$att      else 0
    att_out_pool <- if (n_out > 0L) -res_out$att     else 0
    w_in <- n_in / (n_in + n_out)
    out[k]      <- w_in * att_in + (1 - w_in) * att_out_pool
    out_raw[k]  <- out[k]
    n_inc[k]    <- n_in + n_out
    n_eff[k]    <- (res_in$N_eff %||% 0L) + (res_out$N_eff %||% 0L)
    # Pool delta_norm across directions using the same Neyman weights, then
    # divide DID_k by it (reference: did_multiplegt_main.R:1086-1090 +
    # 1124-1126). Per-direction delta_norm is already a positive magnitude
    # of treatment change (the sign flip is baked into .core_one_event_time
    # via sign_dir), so we pool directly without flipping. Note: the
    # reference does NOT normalize Av_tot_eff (ATE), only the per-event-time
    # DIDs. We therefore keep an unnormalized copy in `out_raw` and use that
    # for the ATE computation upstream.
    if (isTRUE(normalized)) {
      dn_in  <- if (n_in  > 0L) res_in$delta_norm  else NA_real_
      dn_out <- if (n_out > 0L) res_out$delta_norm else NA_real_
      delta_k <- if (n_in == 0L) dn_out
                 else if (n_out == 0L) dn_in
                 else w_in * dn_in + (1 - w_in) * dn_out
      delta_D[k] <- delta_k
      if (!is.na(delta_k) && delta_k != 0) {
        out[k] <- out[k] / delta_k
      } else {
        out[k] <- NA_real_
      }
    }
  }
  list(effects = out, effects_raw = out_raw,
       n_inc = n_inc, n_eff = n_eff, delta_D = delta_D)
}


# -------- trends_lin variant: cumulative FD effects --------

#' Compute trends_lin (linear cohort trends) effects
#'
#' When `trends_lin = TRUE`, the panel has already been first-differenced
#' in `.prep_panel`. The published Effect_k is the CUMULATIVE SUM of
#' per-event-time DIDs from j = 1 to k, computed with `same_switchers = TRUE`
#' (gated on whether the switcher qualifies at every j in 1..k).
#'
#' Reference: did_multiplegt_main.R:831-857 calls
#' `did_multiplegt_dyn_core(..., effects = i, same_switchers = TRUE)`
#' separately for each i in 1..effects. We mirror that by recomputing
#' still_switcher_XX per outer-k and summing U_g across the inner j.
#'
#' @keywords internal
#' @noRd
.compute_effects_trends_lin <- function(prepped, effects, switchers = "",
                                          prefit = NULL,
                                          only_never_switchers = FALSE,
                                          normalized = FALSE) {
  out     <- numeric(effects)
  out_raw <- numeric(effects)
  n_inc   <- integer(effects)
  n_eff   <- integer(effects)
  delta_D <- numeric(effects)
  G <- length(unique(prepped$group_XX))

  for (k in seq_len(effects)) {
    # Recompute still_switcher_XX for this outer-k (effects = k).
    if ("still_switcher_XX" %in% names(prepped)) {
      prepped[, still_switcher_XX := NULL]
    }
    .compute_still_switcher(prepped, effects = k,
                             only_never_switchers = only_never_switchers)

    # Per-direction cumulative U_g across j = 1..k.
    res_in  <- if (switchers != "out") .accumulate_u_g(prepped, k_max = k, direction = 1L, prefit = prefit, only_never_switchers = only_never_switchers, normalized = normalized)
               else list(att = NA_real_, N_inc = 0L, N_eff = 0L, delta_norm = NA_real_)
    res_out <- if (switchers != "in")  .accumulate_u_g(prepped, k_max = k, direction = 0L, prefit = prefit, only_never_switchers = only_never_switchers, normalized = normalized)
               else list(att = NA_real_, N_inc = 0L, N_eff = 0L, delta_norm = NA_real_)

    n_in  <- res_in$N_inc
    n_out <- res_out$N_inc
    if (n_in + n_out == 0L) {
      out[k] <- NA_real_; out_raw[k] <- NA_real_
      n_inc[k] <- 0L; n_eff[k] <- 0L
      delta_D[k] <- NA_real_
      next
    }
    att_in       <- if (n_in  > 0L)  res_in$att      else 0
    att_out_pool <- if (n_out > 0L) -res_out$att     else 0
    w_in <- n_in / (n_in + n_out)
    out[k]      <- w_in * att_in + (1 - w_in) * att_out_pool
    out_raw[k]  <- out[k]
    n_inc[k]    <- n_in + n_out
    n_eff[k]    <- (res_in$N_eff %||% 0L) + (res_out$N_eff %||% 0L)
    if (isTRUE(normalized)) {
      dn_in  <- if (n_in  > 0L) res_in$delta_norm  else NA_real_
      dn_out <- if (n_out > 0L) res_out$delta_norm else NA_real_
      delta_k <- if (n_in == 0L) dn_out
                 else if (n_out == 0L) dn_in
                 else w_in * dn_in + (1 - w_in) * dn_out
      delta_D[k] <- delta_k
      if (!is.na(delta_k) && delta_k != 0) {
        out[k] <- out[k] / delta_k
      } else {
        out[k] <- NA_real_
      }
    }
  }
  list(effects = out, effects_raw = out_raw,
       n_inc = n_inc, n_eff = n_eff, delta_D = delta_D)
}


# Placebo equivalent of .compute_effects_trends_lin. Reference:
# did_multiplegt_main.R:882-906 calls
# did_multiplegt_dyn_core(..., effects = i, placebo = i,
#                          same_switchers = TRUE, same_switchers_pl = TRUE)
# separately for each i in 1..placebo. We mirror by recomputing
# still_switcher_XX per outer-k and summing U_g across the inner j.
#' @keywords internal
#' @noRd
.compute_placebos_trends_lin <- function(prepped, placebo, switchers = "",
                                           prefit = NULL,
                                           only_never_switchers = FALSE,
                                           normalized = FALSE) {
  if (placebo == 0L) return(list(placebos = numeric(0), n_inc = integer(0),
                                  n_eff = integer(0),
                                  delta_D = numeric(0)))
  out     <- numeric(placebo)
  n_inc   <- integer(placebo)
  n_eff   <- integer(placebo)
  delta_D <- numeric(placebo)
  G <- length(unique(prepped$group_XX))

  for (k in seq_len(placebo)) {
    if ("still_switcher_XX" %in% names(prepped)) {
      prepped[, still_switcher_XX := NULL]
    }
    .compute_still_switcher(prepped, effects = k,
                             only_never_switchers = only_never_switchers)

    res_in  <- if (switchers != "out") .accumulate_u_g_placebo(prepped, k_max = k, direction = 1L, prefit = prefit, only_never_switchers = only_never_switchers, normalized = normalized)
               else list(att = NA_real_, N_inc = 0L, N_eff = 0L, delta_norm = NA_real_)
    res_out <- if (switchers != "in")  .accumulate_u_g_placebo(prepped, k_max = k, direction = 0L, prefit = prefit, only_never_switchers = only_never_switchers, normalized = normalized)
               else list(att = NA_real_, N_inc = 0L, N_eff = 0L, delta_norm = NA_real_)

    n_in  <- res_in$N_inc
    n_out <- res_out$N_inc
    if (n_in + n_out == 0L) {
      out[k] <- NA_real_; n_inc[k] <- 0L; n_eff[k] <- 0L
      delta_D[k] <- NA_real_
      next
    }
    att_in       <- if (n_in  > 0L)  res_in$att      else 0
    att_out_pool <- if (n_out > 0L) -res_out$att     else 0
    w_in <- n_in / (n_in + n_out)
    out[k]   <- w_in * att_in + (1 - w_in) * att_out_pool
    n_inc[k] <- n_in + n_out
    n_eff[k] <- (res_in$N_eff %||% 0L) + (res_out$N_eff %||% 0L)
    if (isTRUE(normalized)) {
      dn_in  <- if (n_in  > 0L) res_in$delta_norm  else NA_real_
      dn_out <- if (n_out > 0L) res_out$delta_norm else NA_real_
      delta_k <- if (n_in == 0L) dn_out
                 else if (n_out == 0L) dn_in
                 else w_in * dn_in + (1 - w_in) * dn_out
      delta_D[k] <- delta_k
      if (!is.na(delta_k) && delta_k != 0) {
        out[k] <- out[k] / delta_k
      } else {
        out[k] <- NA_real_
      }
    }
  }
  list(placebos = out, n_inc = n_inc, n_eff = n_eff, delta_D = delta_D)
}


# Compute sum_{j=1..k_max} U_g[j] (placebo variant) for one direction.
#' @keywords internal
#' @noRd
.accumulate_u_g_placebo <- function(prepped, k_max, direction = 1L,
                                      prefit = NULL,
                                      only_never_switchers = FALSE,
                                      normalized = FALSE) {
  G <- length(unique(prepped$group_XX))
  cum_U <- numeric(G)
  any_valid <- FALSE
  last_N_inc <- 0L
  last_N_eff <- 0L
  last_delta_norm <- NA_real_
  for (j in seq_len(k_max)) {
    r <- .core_one_placebo(prepped, k = j, direction = direction,
                             prefit = prefit,
                             only_never_switchers = only_never_switchers,
                             normalized = normalized)
    if (r$N_inc > 0L) {
      any_valid <- TRUE
      cum_U <- cum_U + r$U_g
    }
    if (j == k_max) {
      last_N_inc <- r$N_inc %||% 0L
      last_N_eff <- r$N_eff %||% 0L
      if (isTRUE(normalized)) last_delta_norm <- r$delta_norm
    }
  }
  if (!any_valid) {
    return(list(att = NA_real_, N_inc = 0L, N_eff = 0L,
                 delta_norm = NA_real_))
  }
  list(att = sum(cum_U) / G,
       N_inc = as.integer(last_N_inc),
       N_eff = as.integer(last_N_eff),
       delta_norm = last_delta_norm)
}


# Compute sum_{j=1..k_max} U_g[j] for one direction, then att = sum(U_g) / G.
# Reused .core_one_event_time for the per-j scoring (it already returns U_g).
#' @keywords internal
#' @noRd
.accumulate_u_g <- function(prepped, k_max, direction = 1L, prefit = NULL,
                              only_never_switchers = FALSE,
                              normalized = FALSE) {
  G <- length(unique(prepped$group_XX))
  cum_U <- numeric(G)
  any_valid <- FALSE
  # Reported counts: the reference publishes N_k = count of switcher
  # cells at event-time k (with same_switchers gated on effects=k), so
  # we use the LAST j's count rather than the sum across j. Same logic
  # for delta_norm under normalized=TRUE (reference: main.R:870-874
  # uses delta_norm_i_XX = delta_norm at j=i, not the sum across j).
  last_N_inc <- 0L
  last_N_eff <- 0L
  last_delta_norm <- NA_real_
  for (j in seq_len(k_max)) {
    r <- .core_one_event_time(prepped, k = j, direction = direction,
                                prefit = prefit,
                                only_never_switchers = only_never_switchers,
                                normalized = normalized)
    if (r$N_inc > 0L) {
      any_valid <- TRUE
      cum_U <- cum_U + r$U_g
    }
    if (j == k_max) {
      last_N_inc <- r$N_inc %||% 0L
      last_N_eff <- r$N_eff %||% 0L
      if (isTRUE(normalized)) last_delta_norm <- r$delta_norm
    }
  }
  if (!any_valid) {
    return(list(att = NA_real_, N_inc = 0L, N_eff = 0L,
                 delta_norm = NA_real_))
  }
  list(att = sum(cum_U) / G,
       N_inc = as.integer(last_N_inc),
       N_eff = as.integer(last_N_eff),
       delta_norm = last_delta_norm)
}


# -------- backend hook: wire the r-port into the backend dispatcher --------

#' Compute the clamped horizons the reference would use
#'
#' Matches did_multiplegt_main.R:526-557 (switchers == "" branch with
#' no trends_lin). Returns the actually-feasible (`l_eff`, `l_pl`)
#' subject to:
#'   l_eff <= effects
#'   l_pl  <= placebo
#'   l_pl  <= effects        (definitional: placebos at k > effects are not
#'                            comparable to any computed effect)
#'   l_eff <= max(L_u, L_a)  (data-availability for the dist mask)
#'   l_pl  <= max(L_placebo_u, L_placebo_a)
#'
#' @keywords internal
#' @noRd
.clamp_horizons <- function(prepped, effects, placebo, switchers = "") {
  # One row per group for the group-level quantities.
  g <- unique(prepped[, list(group_XX, S_g_XX, L_g_XX, L_g_placebo_XX)])

  L_u <- suppressWarnings(max(g$L_g_XX[!is.na(g$S_g_XX) & g$S_g_XX == 1L], na.rm = TRUE))
  L_a <- suppressWarnings(max(g$L_g_XX[!is.na(g$S_g_XX) & g$S_g_XX == 0L], na.rm = TRUE))
  if (!is.finite(L_u)) L_u <- 0L
  if (!is.finite(L_a)) L_a <- 0L

  # Reference main.R:526-557: cap by max horizon of the requested
  # switcher direction(s).
  L_for_eff <- switch(switchers,
                      "in"  = L_u,
                      "out" = L_a,
                      max(L_u, L_a, 0L))
  l_eff <- min(as.integer(effects), L_for_eff)

  if (placebo == 0L) {
    l_pl <- 0L
  } else {
    Lp_u <- suppressWarnings(max(g$L_g_placebo_XX[!is.na(g$S_g_XX) & g$S_g_XX == 1L], na.rm = TRUE))
    Lp_a <- suppressWarnings(max(g$L_g_placebo_XX[!is.na(g$S_g_XX) & g$S_g_XX == 0L], na.rm = TRUE))
    if (!is.finite(Lp_u)) Lp_u <- 0L
    if (!is.finite(Lp_a)) Lp_a <- 0L
    L_for_pl <- switch(switchers,
                       "in"  = Lp_u,
                       "out" = Lp_a,
                       max(Lp_u, Lp_a, 0L))
    l_pl <- min(as.integer(placebo), L_for_pl, as.integer(effects))
  }

  list(l_eff = as.integer(l_eff), l_pl = as.integer(l_pl))
}


.backend_r_impl <- function() {
  function(df, args, iter_seed) {
    df_use <- if (iter_seed == 0L) df else .cluster_resample(df, args, iter_seed)

    t0 <- Sys.time()
    tl <- isTRUE(args$trends_lin)
    prepped <- .prep_panel(df_use, args$outcome, args$group, args$time, args$treatment,
                           controls = args$controls, weight = args$weight,
                           trends_nonparam = args$trends_nonparam,
                           dont_drop_larger_lower = isTRUE(args$dont_drop_larger_lower),
                           continuous = args$continuous,
                           trends_lin = tl)
    # Auto-include polynomial features (from continuous=) in the FWL
    # control set. Reference: main.R:202-208 + the controls-residualization
    # block — the polynomial baseline becomes additional regressors.
    auto_ctrl <- attr(prepped, "didgpu_auto_controls")
    eff_controls <- unique(c(args$controls, auto_ctrl))
    prefit <- if (length(eff_controls) > 0L) {
      .prefit_controls(prepped, eff_controls)
    } else NULL
    # Auto-clamp to the feasible horizons, matching the reference.
    sw   <- args$switchers %||% ""
    ons  <- isTRUE(args$only_never_switchers)
    # trends_lin forces same_switchers = TRUE (reference: main.R:842).
    ss   <- isTRUE(args$same_switchers) || tl
    sspl <- isTRUE(args$same_switchers_pl)
    nrm  <- isTRUE(args$normalized)
    h <- .clamp_horizons(prepped, args$effects, args$placebo, switchers = sw)
    if (tl) {
      ce <- .compute_effects_trends_lin(prepped, h$l_eff, switchers = sw,
                                          prefit = prefit,
                                          only_never_switchers = ons,
                                          normalized = nrm)
      # Placebos under trends_lin: same accumulator pattern as effects,
      # summing per-event-time placebo U_g across j = 1..k. Reference:
      # main.R:882-906.
      cp <- .compute_placebos_trends_lin(prepped, h$l_pl, switchers = sw,
                                          prefit = prefit,
                                          only_never_switchers = ons,
                                          normalized = nrm)
    } else {
      ce <- .compute_effects(prepped, h$l_eff, switchers = sw, prefit = prefit,
                              only_never_switchers = ons,
                              same_switchers = ss,
                              normalized = nrm)
      cp <- .compute_placebos(prepped, h$l_pl, switchers = sw, prefit = prefit,
                               only_never_switchers = ons,
                               normalized = nrm,
                               same_switchers_pl = sspl)
    }
    wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    # ATE = N_inc-weighted average of per-event-time DIDs.
    # Derivation (binary, delta_D = 1):
    #   ATE = sum_g (sum_k w_k * U_Gg_k_g) / G
    #       = sum_k w_k * (sum_g U_Gg_k_g) / G
    #       = sum_k w_k * DID_k
    # where w_k = N_inc_k / sum_j N_inc_j.
    # For the both-directions case, the per-event-time DID_k_published
    # already includes the cross-direction Neyman pooling, and the
    # cumulative N_inc weight matches (see notes in .compute_effects).
    #
    # IMPORTANT: when normalized=TRUE the reference does NOT normalize the
    # ATE (Av_tot_eff stays the same; reference main.R:1147-1163 builds it
    # from per-row U_Gg contributions that are summed across event-times
    # without ever being divided by delta_D). We mirror that by computing
    # the ATE from `effects_raw` (pre-normalization) rather than `effects`.
    ate_src <- if (isTRUE(nrm)) ce$effects_raw else ce$effects
    ate <- if (tl) {
      # Reference suppresses ATE under trends_lin (the linear-trends
      # identification doesn't pin down a single cumulative average).
      NA_real_
    } else if (h$l_eff == 0L || sum(ce$n_inc) == 0L) {
      NA_real_
    } else {
      valid <- !is.na(ate_src) & ce$n_inc > 0L
      if (!any(valid)) NA_real_ else
        sum(ate_src[valid] * ce$n_inc[valid]) / sum(ce$n_inc[valid])
    }

    # predict_het: post-fit heterogeneity regression. Only run on the
    # point-estimate iter (iter_seed == 0); bootstrap iters skip it to
    # save time (we don't bootstrap the het coefficients separately
    # because the reference's reported SEs are analytical HC1 anyway).
    het_block <- if (!is.null(args$predict_het) && iter_seed == 0L) {
      het_vars    <- unlist(args$predict_het[[1L]])
      het_effects <- as.integer(unlist(args$predict_het[[2L]]))
      tryCatch(
        .compute_predict_het(prepped, het_vars, het_effects,
                              l_eff = h$l_eff,
                              trends_nonparam_col = args$trends_nonparam,
                              ci_level = args$ci_level %||% 95),
        error = function(e) {
          warning("predict_het failed: ", conditionMessage(e))
          NULL
        })
    } else NULL

    list(
      effects        = ce$effects,
      ate            = ate,
      placebos       = cp$placebos,
      n_effects      = h$l_eff,
      n_placebos     = h$l_pl,
      n_inc_effects  = ce$n_inc,
      n_inc_placebos = cp$n_inc,
      n_eff_effects  = ce$n_eff,
      n_eff_placebos = cp$n_eff,
      predict_het    = het_block,
      iter_seed      = as.integer(iter_seed),
      wall_seconds   = wall,
      backend        = "r"
    )
  }
}
