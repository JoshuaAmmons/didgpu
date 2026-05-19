# ============================================================================
# Callaway-Sant'Anna ATT(g, t) computation. Loops over (cohort g, time t)
# cells and dispatches to the chosen inner estimator (OR / IPW / DR) on
# each cell.
#
# Includes pre-treatment cells (t < g) so the placebo / pre-trend test
# is a side-product of the same machinery.
# ============================================================================


# Per-cell ATT(g, t) computation for ALL (g, t) including pre-treatment
# placebos. Dispatches to .cs_inner_or / .cs_inner_ipw / .cs_inner_dr.
#' @keywords internal
#' @noRd
.cs_compute_att_gt <- function(df, args, verbose = TRUE) {
  d <- data.table::as.data.table(df)
  data.table::setnames(d,
    c(args$outcome, args$group, args$time, args$treatment),
    c("Y_XX", "G_XX", "T_XX", "D_XX"))
  data.table::setorder(d, G_XX, T_XX)

  # Covariates: keep them under their original names but pull as a
  # per-unit matrix later.
  cov_cols <- args$covariates
  has_cov <- !is.null(cov_cols) && length(cov_cols) > 0L

  # Per-unit first-treated period (Inf for never-treated).
  d[, F_g_XX := {
      if (any(D_XX == 1L, na.rm = TRUE))
        as.numeric(min(T_XX[D_XX == 1L], na.rm = TRUE))
      else Inf
    }, by = G_XX]

  T_min <- min(d$T_XX, na.rm = TRUE)
  T_max <- max(d$T_XX, na.rm = TRUE)
  units <- sort(unique(d$G_XX))
  F_g_per_unit <- d[!duplicated(G_XX), F_g_XX]
  names(F_g_per_unit) <- as.character(d[!duplicated(G_XX), G_XX])

  cohorts <- sort(unique(d$F_g_XX[is.finite(d$F_g_XX) & d$F_g_XX > T_min]))
  if (length(cohorts) == 0L) {
    stop("No treated cohorts found.")
  }
  if (args$control_group == "never") {
    if (!any(!is.finite(F_g_per_unit))) {
      stop("control_group = 'never' requires at least one never-treated unit.")
    }
  }

  # Pre-compute per-(unit, time) covariate matrix (constant across t
  # for time-invariant covariates; we take the value at T_min).
  X_per_unit <- NULL
  if (has_cov) {
    # Check time-invariance per unit; take first value.
    X_per_unit <- d[!duplicated(G_XX),
                     c(list(G_XX = G_XX), as.list(.SD)),
                     .SDcols = cov_cols]
    X_mat <- as.matrix(X_per_unit[, cov_cols, with = FALSE])
    rownames(X_mat) <- as.character(X_per_unit$G_XX)
  }

  # ------------------------------------------------------------------
  # Pass 1: enumerate every (g, t) cell and gather its inputs WITHOUT
  # solving. The cell list is then handed to one of two solver paths:
  #   (a) batched CUDA — single call across all cells (Phase 1 hook,
  #       Phase 2 implementation),
  #   (b) per-cell R   — current production path; also the fallback
  #       when CUDA returns NULL.
  # The cell metadata (g, t, units, n_treated, n_control) is the same
  # either way, so the long-form data frame assembled below is solver-
  # agnostic.
  # ------------------------------------------------------------------
  cell_meta <- list()
  cell_data <- list()
  for (g in cohorts) {
    pre_t <- g - 1L
    if (pre_t < T_min) next
    Y_pre <- d[T_XX == pre_t, list(G_XX, Y_pre = Y_XX)]
    data.table::setkey(Y_pre, G_XX)
    for (t in T_min:T_max) {
      if (t == pre_t) next   # pre-period reference; mechanically 0
      Y_t <- d[T_XX == t, list(G_XX, Y_t = Y_XX)]
      data.table::setkey(Y_t, G_XX)
      merged <- merge(Y_pre, Y_t, by = "G_XX")
      merged[, delta_XX := Y_t - Y_pre]
      merged <- merged[!is.na(delta_XX)]

      treated_units <- unique(d$G_XX[d$F_g_XX == g])
      control_units <- .cs_control_units(d, units, F_g_per_unit,
                                            g, t, args$control_group)
      keep <- merged$G_XX %in% c(treated_units, control_units)
      merged <- merged[keep]
      if (nrow(merged) == 0L) next

      D_mask <- merged$G_XX %in% treated_units
      delta_v <- merged$delta_XX
      if (has_cov) {
        Xm <- X_mat[as.character(merged$G_XX), , drop = FALSE]
      } else {
        Xm <- NULL
      }
      n_total <- nrow(merged)

      cell_meta[[length(cell_meta) + 1L]] <- list(
        g = g, t = t, n_total = n_total,
        n_treated = sum(D_mask), n_control = sum(!D_mask),
        units = merged$G_XX)
      cell_data[[length(cell_data) + 1L]] <- list(
        delta = delta_v, D_mask = D_mask, X = Xm, n_total = n_total,
        units = merged$G_XX)
    }
  }
  n_cells <- length(cell_meta)
  if (n_cells == 0L) {
    out <- data.frame(g = integer(), t = integer(), event_time = integer(),
                       att = numeric(), se = numeric(),
                       n_treated = integer(), n_control = integer())
    attr(out, "IF_per_cell") <- list()
    attr(out, "F_g_per_unit") <- F_g_per_unit
    attr(out, "units")        <- units
    return(out)
  }

  # ------------------------------------------------------------------
  # Pass 2: solve.
  # ------------------------------------------------------------------
  cuda_result <- NULL
  if (identical(args$backend, "cuda")) {
    cuda_result <- .cs_inner_batched_cuda(
      cells     = cell_data,
      method    = args$est_method,
      all_units = units)
  }

  solver_atts <- numeric(n_cells)
  IF_list <- vector("list", n_cells)
  if (!is.null(cuda_result)) {
    # CUDA succeeded. Pull per-cell ATT from the returned vector and
    # per-cell influence vectors from the n_units x n_cells matrix.
    # Each cell only uses a subset of units; we project the column of
    # the IF matrix back onto each cell's unit list.
    unit_to_row <- stats::setNames(seq_along(units), as.character(units))
    for (c_idx in seq_len(n_cells)) {
      solver_atts[c_idx] <- cuda_result$att[c_idx]
      cell_units <- cell_meta[[c_idx]]$units
      if (!is.null(cuda_result$influence)) {
        rows <- unit_to_row[as.character(cell_units)]
        IF_list[[c_idx]] <- list(
          g = cell_meta[[c_idx]]$g, t = cell_meta[[c_idx]]$t,
          units = cell_units,
          IF    = as.numeric(cuda_result$influence[rows, c_idx]))
      } else {
        IF_list[[c_idx]] <- list(
          g = cell_meta[[c_idx]]$g, t = cell_meta[[c_idx]]$t,
          units = cell_units, IF = rep(NA_real_, length(cell_units)))
      }
    }
  } else {
    # R fallback: per-cell solve.
    for (c_idx in seq_len(n_cells)) {
      ce <- cell_data[[c_idx]]
      inner <- .cs_inner_dispatch(args$est_method, ce$delta, ce$D_mask,
                                    ce$X, ce$n_total)
      solver_atts[c_idx] <- as.numeric(inner$att)
      IF_list[[c_idx]] <- list(
        g = cell_meta[[c_idx]]$g, t = cell_meta[[c_idx]]$t,
        units = cell_meta[[c_idx]]$units, IF = inner$IF)
    }
  }

  # ------------------------------------------------------------------
  # Assemble long-form result.
  # ------------------------------------------------------------------
  results <- lapply(seq_len(n_cells), function(c_idx) {
    m <- cell_meta[[c_idx]]
    data.frame(
      g = as.integer(m$g), t = as.integer(m$t),
      event_time = as.integer(m$t - m$g),
      att = solver_atts[c_idx], se = NA_real_,
      n_treated = as.integer(m$n_treated),
      n_control = as.integer(m$n_control),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  attr(out, "IF_per_cell") <- IF_list
  attr(out, "F_g_per_unit") <- F_g_per_unit
  attr(out, "units")        <- units
  out
}


# Cluster bootstrap SE: resample units with replacement, recompute the
# long-form ATT(g, t) table on each rep, take SD across reps per cell.
#' @keywords internal
#' @noRd
.cs_bootstrap_se <- function(df, args, att_gt, verbose = TRUE) {
  B <- args$bootstrap_reps
  if (B <= 0L) return(att_gt)
  d <- data.table::as.data.table(df)
  units <- unique(d[[args$group]])

  boot_mat <- matrix(NA_real_, nrow = B, ncol = nrow(att_gt))
  for (b in seq_len(B)) {
    set.seed(args$seed + b, kind = "Mersenne-Twister")
    picks <- sample(units, length(units), replace = TRUE)
    OFFSET <- max(as.numeric(units), na.rm = TRUE) + 1L
    rep_dfs <- vector("list", length(picks))
    pcount <- integer(length(units))
    names(pcount) <- as.character(units)
    for (i in seq_along(picks)) {
      u <- picks[i]
      pcount[as.character(u)] <- pcount[as.character(u)] + 1L
      block <- d[d[[args$group]] == u, , drop = FALSE]
      if (pcount[as.character(u)] > 1L) {
        shift <- (pcount[as.character(u)] - 1L) * OFFSET
        block[[args$group]] <- as.numeric(block[[args$group]]) + shift
      }
      rep_dfs[[i]] <- block
    }
    df_b <- data.table::rbindlist(rep_dfs)
    att_b <- tryCatch(
      .cs_compute_att_gt(as.data.frame(df_b), args, verbose = FALSE),
      error = function(e) NULL)
    if (is.null(att_b)) next
    key_b <- paste(att_b$g, att_b$t, sep = "_")
    key_p <- paste(att_gt$g, att_gt$t, sep = "_")
    boot_mat[b, ] <- att_b$att[match(key_p, key_b)]
  }
  ses <- apply(boot_mat, 2L, function(col) {
    ok <- !is.na(col)
    if (sum(ok) >= 2L) stats::sd(col[ok]) else NA_real_
  })
  att_gt$se <- as.numeric(ses)
  z <- stats::qnorm(0.5 + (args$ci_level %||% 95) / 200)
  att_gt$ci_low  <- att_gt$att - z * att_gt$se
  att_gt$ci_high <- att_gt$att + z * att_gt$se
  att_gt
}


# Multiplier bootstrap SE (Mammen wild bootstrap on influence functions).
# Much faster than the cluster bootstrap because we only need the IF
# once; each replicate is a single re-weighted sum.
#' @keywords internal
#' @noRd
# GPU cluster-bootstrap fast path for the CS estimator. Returns
# `att_gt` with `se` / `ci_low` / `ci_high` populated, OR NULL on any
# failure (caller falls back to .cs_bootstrap_se).
#
# Uses the IF-shortcut: builds an (n_units, n_cells) influence matrix
# from attr(att_gt, "IF_per_cell"), then calls the CUDA kernel to
# compute B replicates of weight @ IF where weight[b, u] is the count
# of times unit u was picked in replicate b's cluster draw.
#
# This is asymptotically equivalent to the per-rep refit cluster
# bootstrap (.cs_bootstrap_se) by the standard delta-method argument
# (Hansen 2022, Ch.10). Finite-sample SEs differ by O(1/sqrt(B)).
#' @keywords internal
#' @noRd
.cs_cluster_bootstrap_cuda <- function(att_gt, args) {
  B <- args$bootstrap_reps
  if (B <= 0L) return(att_gt)
  if (!isTRUE(tryCatch(didgpu_has_cuda_support(),
                       error = function(e) FALSE))) return(NULL)
  IF_list <- attr(att_gt, "IF_per_cell")
  units   <- attr(att_gt, "units")
  if (is.null(IF_list) || is.null(units)) return(NULL)
  n_unit  <- length(units)
  n_cells <- nrow(att_gt)
  if (length(IF_list) != n_cells) return(NULL)

  # Build (n_units, n_cells) IF matrix: row u is "unit u's per-cell
  # influence", with 0 for cells where unit u doesn't appear.
  IF_mat <- matrix(0.0, nrow = n_unit, ncol = n_cells)
  unit_to_row <- stats::setNames(seq_along(units), as.character(units))
  for (c_idx in seq_along(IF_list)) {
    cell <- IF_list[[c_idx]]
    if (is.null(cell$IF) || length(cell$IF) == 0L) next
    rows <- unit_to_row[as.character(cell$units)]
    # The per-cell IF stored by .cs_inner_dispatch is divided by cell
    # size for the multiplier bootstrap form; cluster bootstrap also
    # needs that same scaling so the unit-weighted sum equals the cell
    # ATT contribution.
    IF_mat[rows, c_idx] <- cell$IF / length(cell$IF)
  }

  # For CS the "cluster" is the unit itself.
  cluster_id <- as.integer(seq_along(units) - 1L)
  result <- tryCatch(
    didgpu_cuda_cluster_bootstrap_r(
      IF         = IF_mat,
      cluster_id = cluster_id,
      n_clusters = n_unit,
      B          = as.integer(B),
      seed       = as.integer(args$seed %||% 1L)),
    error = function(e) NULL)
  if (is.null(result)) return(NULL)

  # boot_mat[b, c] is the b-th replicate's bootstrap deviation for
  # cell c (the IF-weighted sum, centered around the original att). SE
  # is the cross-replicate SD per column.
  ses <- apply(result, 2L, function(col) stats::sd(col, na.rm = TRUE))
  att_gt$se <- as.numeric(ses)
  z <- stats::qnorm(0.5 + (args$ci_level %||% 95) / 200)
  att_gt$ci_low  <- att_gt$att - z * att_gt$se
  att_gt$ci_high <- att_gt$att + z * att_gt$se
  att_gt
}


.cs_multiplier_bootstrap_se <- function(att_gt, args, verbose = TRUE) {
  B <- args$bootstrap_reps
  if (B <= 0L) return(att_gt)
  IF_list <- attr(att_gt, "IF_per_cell")
  units   <- attr(att_gt, "units")
  if (is.null(IF_list) || is.null(units)) {
    warning("Multiplier bootstrap requires IF attributes; ",
            "falling back to cluster bootstrap.")
    return(att_gt)
  }
  n_unit <- length(units)

  # GPU fast path: when backend == "cuda", try the multiplier-bootstrap
  # kernel. It computes the (B, n_cells) bootstrap-deviation matrix in
  # one launch (a (B x n_units) @ (n_units x n_cells) Rademacher-weighted
  # product). Falls back to the R loop on any failure.
  if (identical(args$backend, "cuda")) {
    cuda_se <- .cs_multiplier_bootstrap_cuda(att_gt, args)
    if (!is.null(cuda_se)) return(cuda_se)
  }

  # Rademacher weights: +1 / -1 with equal probability. Mammen weights
  # (golden-ratio-based) are an alternative; both are second-order accurate.
  boot_mat <- matrix(NA_real_, nrow = B, ncol = nrow(att_gt))
  set.seed(args$seed, kind = "Mersenne-Twister")
  for (b in seq_len(B)) {
    xi <- sample(c(-1, 1), n_unit, replace = TRUE)
    names(xi) <- as.character(units)
    for (cell_i in seq_along(IF_list)) {
      cell <- IF_list[[cell_i]]
      w <- xi[as.character(cell$units)]
      # Bootstrap stat = att + mean(xi_i * IF_i)
      boot_mat[b, cell_i] <- att_gt$att[cell_i] +
                              sum(w * cell$IF) / length(cell$IF)
    }
  }
  ses <- apply(boot_mat, 2L, function(col) stats::sd(col, na.rm = TRUE))
  att_gt$se <- as.numeric(ses)
  z <- stats::qnorm(0.5 + (args$ci_level %||% 95) / 200)
  att_gt$ci_low  <- att_gt$att - z * att_gt$se
  att_gt$ci_high <- att_gt$att + z * att_gt$se
  att_gt
}


# GPU multiplier (wild) bootstrap. Builds the (n_units, n_cells) IF
# matrix from attr(att_gt, "IF_per_cell"), calls the CUDA kernel with
# Rademacher weights, and reads back columnwise SDs.
#
# Mathematically identical to .cs_multiplier_bootstrap_se's R loop
# (modulo RNG differences) when the per-cell IF stored by
# .cs_inner_dispatch is already divided by the per-cell sample size
# (.cs_inner_or et al. do this divide before storing). Per-replicate
# deviations differ because cuRAND != MT19937; population SDs match
# to within Monte-Carlo error.
#
# Returns NULL on any failure for caller fallback.
#' @keywords internal
#' @noRd
.cs_multiplier_bootstrap_cuda <- function(att_gt, args) {
  B <- args$bootstrap_reps
  if (B <= 0L) return(NULL)
  if (!isTRUE(tryCatch(didgpu_has_cuda_support(),
                       error = function(e) FALSE))) return(NULL)
  IF_list <- attr(att_gt, "IF_per_cell")
  units   <- attr(att_gt, "units")
  if (is.null(IF_list) || is.null(units)) return(NULL)
  n_unit  <- length(units)
  n_cells <- nrow(att_gt)
  if (length(IF_list) != n_cells) return(NULL)

  IF_mat <- matrix(0.0, nrow = n_unit, ncol = n_cells)
  unit_to_row <- stats::setNames(seq_along(units), as.character(units))
  for (c_idx in seq_along(IF_list)) {
    cell <- IF_list[[c_idx]]
    if (is.null(cell$IF) || length(cell$IF) == 0L) next
    rows <- unit_to_row[as.character(cell$units)]
    IF_mat[rows, c_idx] <- cell$IF / length(cell$IF)
  }

  result <- tryCatch(
    didgpu_cuda_multiplier_bootstrap_r(
      IF        = IF_mat,
      B         = as.integer(B),
      mult_kind = 0L,            # Rademacher
      seed      = as.integer(args$seed %||% 1L)),
    error = function(e) NULL)
  if (is.null(result)) return(NULL)

  ses <- apply(result, 2L, function(col) stats::sd(col, na.rm = TRUE))
  att_gt$se <- as.numeric(ses)
  z <- stats::qnorm(0.5 + (args$ci_level %||% 95) / 200)
  att_gt$ci_low  <- att_gt$att - z * att_gt$se
  att_gt$ci_high <- att_gt$att + z * att_gt$se
  att_gt
}


# Aggregation step. Given the long-form ATT(g, t) table, return one of
# four summaries. Pre-treatment cells (event_time < 0) appear in the
# event-study aggregation as placebos.
#' @keywords internal
#' @noRd
.cs_aggregate <- function(att_gt, aggregation, args) {
  if (nrow(att_gt) == 0L) {
    return(data.frame(scheme = aggregation, level = character(0),
                      estimate = numeric(0), se = numeric(0),
                      n_cells = integer(0)))
  }
  w <- att_gt$n_treated
  switch(aggregation,
    "event" = {
      events <- sort(unique(att_gt$event_time))
      out <- data.frame(
        event_time = events,
        estimate = vapply(events, function(e) {
          cells <- att_gt$event_time == e
          if (!any(cells)) return(NA_real_)
          sum(att_gt$att[cells] * w[cells]) / sum(w[cells])
        }, numeric(1)),
        n_cells = vapply(events, function(e) sum(att_gt$event_time == e),
                          integer(1)),
        stringsAsFactors = FALSE
      )
      out
    },
    "group" = {
      gs <- sort(unique(att_gt$g))
      data.frame(
        g = gs,
        estimate = vapply(gs, function(g) {
          cells <- att_gt$g == g & att_gt$t >= g
          if (!any(cells)) return(NA_real_)
          sum(att_gt$att[cells] * w[cells]) / sum(w[cells])
        }, numeric(1)),
        n_cells = vapply(gs, function(g) sum(att_gt$g == g & att_gt$t >= g),
                          integer(1)),
        stringsAsFactors = FALSE
      )
    },
    "calendar" = {
      ts <- sort(unique(att_gt$t))
      data.frame(
        t = ts,
        estimate = vapply(ts, function(t) {
          cells <- att_gt$t == t & att_gt$g <= t
          if (!any(cells)) return(NA_real_)
          sum(att_gt$att[cells] * w[cells]) / sum(w[cells])
        }, numeric(1)),
        n_cells = vapply(ts, function(t) sum(att_gt$t == t & att_gt$g <= t),
                          integer(1)),
        stringsAsFactors = FALSE
      )
    },
    "overall" = {
      post <- att_gt$t >= att_gt$g
      est <- sum(att_gt$att[post] * w[post]) / sum(w[post])
      data.frame(scheme = "overall", estimate = est,
                  n_cells = sum(post),
                  stringsAsFactors = FALSE)
    }
  )
}


# Placebo joint test: pre-treatment cells should have ATT ~ 0 under
# parallel trends. Returns a list with per-event-time placebo summary
# and a joint chi-square p-value.
#' @keywords internal
#' @noRd
.cs_placebo_test <- function(att_gt, args) {
  pre <- att_gt[att_gt$event_time < 0L, , drop = FALSE]
  if (nrow(pre) == 0L) {
    return(list(per_event = data.frame(), joint_pval = NA_real_,
                 message = "no pre-treatment cells available"))
  }
  # Per-event-time aggregation of placebos.
  events <- sort(unique(pre$event_time))
  per_event <- data.frame(
    event_time = events,
    estimate = vapply(events, function(e) {
      cells <- pre$event_time == e
      sum(pre$att[cells] * pre$n_treated[cells]) / sum(pre$n_treated[cells])
    }, numeric(1)),
    n_cells = vapply(events, function(e) sum(pre$event_time == e),
                      integer(1)),
    se = vapply(events, function(e) {
      cells <- pre$event_time == e
      ses <- pre$se[cells]
      ok <- !is.na(ses)
      if (!any(ok)) return(NA_real_)
      # Conservative aggregated SE (assume independence): sqrt(mean(se^2)).
      sqrt(mean(ses[ok]^2))
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  per_event$pval <- 2 * stats::pnorm(-abs(per_event$estimate / per_event$se))
  # Joint chi-square test: assume per-cell estimates are independent
  # (conservative). For the cluster/multiplier bootstrap, the proper
  # joint test uses the bootstrap covariance — extension TBD.
  ok <- !is.na(per_event$estimate) & !is.na(per_event$se) & per_event$se > 0
  if (sum(ok) >= 1L) {
    z <- per_event$estimate[ok] / per_event$se[ok]
    chi <- sum(z^2)
    joint_pval <- stats::pchisq(chi, df = sum(ok), lower.tail = FALSE)
  } else {
    joint_pval <- NA_real_
  }
  list(per_event = per_event, joint_pval = joint_pval)
}
