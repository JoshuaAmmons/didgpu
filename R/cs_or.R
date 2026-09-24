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
  #
  # With `first_treat` (did's `gname`: the period a unit is first treated,
  # 0 or NA if never) the cohort is read from the data, exactly as did
  # does. Without it the cohort is inferred as the first period in which
  # the unit is OBSERVED with treatment == 1 -- which is the same thing on
  # a balanced panel, but not when the adoption-year row itself is
  # missing: the unit then lands in a later cohort. That is harmless when
  # units with gaps are dropped (the default, balancing, mode) and wrong
  # under allow_unbalanced_panel = TRUE, where it is warned about below.
  ft_col <- args$first_treat
  if (!is.null(ft_col)) {
    ft_raw <- as.numeric(df[[ft_col]])
    ft_by_row <- ft_raw[match(paste(d$G_XX, d$T_XX),
                              paste(df[[args$group]], df[[args$time]]))]
    d[, ft_row_XX := ft_by_row]
    d[, F_g_XX := {
        v <- ft_row_XX[!is.na(ft_row_XX)]
        if (!length(v) || v[1L] == 0) Inf else v[1L]
      }, by = G_XX]
    d[, ft_row_XX := NULL]
  } else {
    d[, F_g_XX := {
        if (any(D_XX == 1L, na.rm = TRUE))
          as.numeric(min(T_XX[D_XX == 1L], na.rm = TRUE))
        else Inf
      }, by = G_XX]
  }

  # ---- did's preprocessing (pre_process_did) -------------------------
  # Rows with a missing outcome, treatment or covariate are dropped first.
  keep_cols <- c("Y_XX", "T_XX", "D_XX", if (has_cov) cov_cols)
  d <- d[stats::complete.cases(d[, keep_cols, with = FALSE])]
  T_min <- min(d$T_XX)
  # Units already treated in the first period have no pre-period: did
  # drops them from the data entirely, so they are neither cohorts nor
  # controls, and do not count in n.
  first_treated <- unique(d$G_XX[is.finite(d$F_g_XX) & d$F_g_XX <= T_min])
  if (length(first_treated)) {
    if (isTRUE(verbose)) {
      message(sprintf("[didgpu] Dropped %d units that were already treated in the first period.",
                      length(first_treated)))
    }
    d <- d[!G_XX %in% first_treated]
  }
  tlist <- sort(unique(d$T_XX))
  nT <- length(tlist)
  T_min <- tlist[1L]; T_max <- tlist[nT]
  # A cohort first treated after the last observed period is never
  # treated within the sample (did: asif_never_treated).
  d[is.finite(F_g_XX) & F_g_XX > T_max, F_g_XX := Inf]
  # Balanced? Every unit observed in every period.
  per_unit <- d[, list(n = .N), by = G_XX]
  balanced <- all(per_unit$n == nT)
  allow_unb <- isTRUE(args$allow_unbalanced_panel)
  rc_mode <- !balanced && allow_unb
  if (rc_mode && is.null(ft_col)) {
    # A unit whose adoption period is unobserved gets a later inferred
    # cohort than its true one. Detect exactly that case.
    # Ambiguous exactly when the period just before the first observed
    # treated period is itself unobserved for that unit.
    late <- d[is.finite(F_g_XX), list(miss = {
                 k <- match(F_g_XX[1L], tlist)
                 !is.na(k) && k > 1L && !(tlist[k - 1L] %in% T_XX) }),
              by = G_XX]
    if (any(late$miss, na.rm = TRUE)) {
      warning("allow_unbalanced_panel = TRUE with cohorts inferred from `treatment`: ",
              sum(late$miss, na.rm = TRUE), " treated units have a gap right before ",
              "their first observed treated period, so their true cohort may be earlier. ",
              "did reads the cohort from its gname column; pass first_treat = <that column> ",
              "to match it.", call. = FALSE)
    }
  }
  if (!balanced && !allow_unb) {
    # did's default: convert to a balanced panel by dropping every unit
    # not observed in all periods (BMisc::makeBalancedPanel). didgpu used
    # to keep, cell by cell, whatever units happened to be observed in
    # both periods -- which matched neither of did's modes.
    keep_u <- per_unit$G_XX[per_unit$n == nT]
    n_drop <- nrow(per_unit) - length(keep_u)
    if (isTRUE(verbose)) {
      message(sprintf("[didgpu] Dropped %d units while converting to a balanced panel ", n_drop),
              "(as did::att_gt does by default); pass allow_unbalanced_panel = TRUE ",
              "for did's repeated-cross-section estimators instead.")
    }
    d <- d[G_XX %in% keep_u]
    if (nrow(d) == 0L) stop("All units dropped converting to a balanced panel.")
  }

  units <- sort(unique(d$G_XX))
  n_units <- length(units)
  F_g_per_unit <- d[!duplicated(G_XX), F_g_XX]
  names(F_g_per_unit) <- as.character(d[!duplicated(G_XX), G_XX])
  F_g_per_unit <- F_g_per_unit[as.character(units)]

  cohorts <- sort(unique(d$F_g_XX[is.finite(d$F_g_XX) & d$F_g_XX > T_min]))
  if (length(cohorts) == 0L) {
    stop("No treated cohorts found.")
  }
  if (args$control_group == "never") {
    if (!any(!is.finite(F_g_per_unit))) {
      stop("control_group = 'never' requires at least one never-treated unit.")
    }
  }
  # Cohort sizes in UNITS: did weights every ATT(g,t) by its cohort's
  # share of units (pg in compute.aggte), not by the cell's row count.
  cohort_size <- table(F_g_per_unit[is.finite(F_g_per_unit)])

  X_mat <- NULL
  if (has_cov) {
    X_per_unit <- d[!duplicated(G_XX),
                     c(list(G_XX = G_XX), as.list(.SD)),
                     .SDcols = cov_cols]
    X_mat <- as.matrix(X_per_unit[, cov_cols, with = FALSE])
    rownames(X_mat) <- as.character(X_per_unit$G_XX)
  }

  # ------------------------------------------------------------------
  # Pass 1: enumerate the cells exactly as did::att_gt does
  # (compute.att_gt), for both base periods:
  #   varying   (did's default): a pre-treatment cell compares period t
  #             with the period before it; cells run over every period
  #             but the first.
  #   universal: every cell compares with the last period before g, and
  #             the cell AT that period is reported as exactly 0.
  # Post-treatment cells use the last period before g either way. For the
  # not-yet-treated control group the threshold is the LATER of the two
  # periods compared -- under a universal base that is the base period
  # for an early pre-treatment cell, which didgpu previously got wrong.
  # ------------------------------------------------------------------
  bp <- args$base_period %||% "varying"
  tfac <- if (identical(bp, "universal")) 0L else 1L
  cell_meta <- list()
  cell_data <- list()
  for (g in cohorts) {
    pre_g <- utils::tail(which(tlist < g), 1L)
    if (!length(pre_g)) next
    treated_units <- units[is.finite(F_g_per_unit) & F_g_per_unit == g]
    for (ti in seq_len(nT - tfac)) {
      cur <- tlist[ti + tfac]
      pret <- if (tfac == 0L) pre_g else ti
      if (g <= cur) pret <- pre_g
      base <- tlist[pret]
      if (tfac == 0L && base == cur) {
        cell_meta[[length(cell_meta) + 1L]] <- list(
          g = g, t = cur, n_total = length(treated_units),
          n_treated = length(treated_units), n_control = NA_integer_,
          n_cohort = as.integer(cohort_size[as.character(g)]),
          units = treated_units, zero = TRUE)
        cell_data[[length(cell_data) + 1L]] <- list(zero = TRUE,
          units = treated_units)
        next
      }
      thr <- tlist[max(ti, pret) + tfac]
      control_units <- if (args$control_group == "never") {
        units[!is.finite(F_g_per_unit)]
      } else {
        units[!is.finite(F_g_per_unit) |
              (F_g_per_unit > thr & F_g_per_unit != g)]
      }
      if (!rc_mode) {
        Y_base <- d[T_XX == base, list(G_XX, Y_pre = Y_XX)]
        Y_t    <- d[T_XX == cur,  list(G_XX, Y_t = Y_XX)]
        merged <- merge(Y_base, Y_t, by = "G_XX")
        merged[, delta_XX := Y_t - Y_pre]
        merged <- merged[!is.na(delta_XX)]
        merged <- merged[G_XX %in% c(treated_units, control_units)]
        if (nrow(merged) == 0L) next
        D_mask <- merged$G_XX %in% treated_units
        Xm <- if (has_cov) X_mat[as.character(merged$G_XX), , drop = FALSE] else NULL
        cell_meta[[length(cell_meta) + 1L]] <- list(
          g = g, t = cur, n_total = nrow(merged),
          n_treated = sum(D_mask), n_control = sum(!D_mask),
          n_cohort = as.integer(cohort_size[as.character(g)]),
          units = merged$G_XX)
        cell_data[[length(cell_data) + 1L]] <- list(
          delta = merged$delta_XX, D_mask = D_mask, X = Xm,
          n_total = nrow(merged), units = merged$G_XX)
      } else {
        rows <- d[T_XX %in% c(base, cur) &
                    G_XX %in% c(treated_units, control_units)]
        if (nrow(rows) == 0L) next
        Dr <- as.numeric(rows$G_XX %in% treated_units)
        Xr <- if (has_cov) as.matrix(rows[, cov_cols, with = FALSE]) else NULL
        cell_meta[[length(cell_meta) + 1L]] <- list(
          g = g, t = cur, n_total = nrow(rows),
          n_treated = sum(Dr), n_control = sum(1 - Dr),
          n_cohort = as.integer(cohort_size[as.character(g)]),
          units = unique(rows$G_XX))
        cell_data[[length(cell_data) + 1L]] <- list(
          rc = TRUE, y = rows$Y_XX, post = as.numeric(rows$T_XX == cur),
          D = Dr, C = 1 - Dr, X = Xr, row_unit = rows$G_XX)
      }
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
  special <- any(vapply(cell_data, function(z) isTRUE(z$rc) || isTRUE(z$zero),
                        logical(1)))
  if (!special && identical(args$backend, "cuda") &&
      isTRUE(getOption("didgpu.cs_cuda_inner", FALSE))) {
    # OFF BY DEFAULT, deliberately. The batched CUDA inner kernel still
    # computes the OLD influence functions: treated-arm only, no
    # nuisance-estimation terms, wrong normalisers. Its ATT agrees with
    # the CPU path to ~4e-16, but its influence functions do not -- and
    # they now feed BOTH the multiplier bootstrap and the aggregation
    # SEs. Measured against did::att_gt() after the CPU rewrite, per-cell
    # max |dSE|:
    #     backend "r"    1.11e-16
    #     backend "cpu"  1.11e-16
    #     backend "cuda" 1.96e-01   <- the stale kernel
    #
    # Until src/ implements the DRDID influence functions on device (a
    # per-cell propensity Hessian inverse plus the OLS estimation-effect
    # term), CS routes through the validated CPU path so that every
    # backend returns identical numbers. Correctness before speed --
    # and the GPU is not winning for CS at realistic panel sizes anyway.
    # Re-enable for kernel development with
    # options(didgpu.cs_cuda_inner = TRUE).
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
      if (isTRUE(ce$zero)) {
        # The universal base period's own cell: 0 by construction, with a
        # zero influence function (did reports it the same way).
        solver_atts[c_idx] <- 0
        IF_list[[c_idx]] <- list(
          g = cell_meta[[c_idx]]$g, t = cell_meta[[c_idx]]$t,
          units = ce$units, IF = rep(0, length(ce$units)))
        next
      }
      if (isTRUE(ce$rc)) {
        rc <- .cs_rc_cell(args$est_method, ce$y, ce$post, ce$D, ce$C,
                          ce$X, ce$row_unit, n_units)
        solver_atts[c_idx] <- as.numeric(rc$att)
        IF_list[[c_idx]] <- list(
          g = cell_meta[[c_idx]]$g, t = cell_meta[[c_idx]]$t,
          units = rc$units, IF = rc$IF)
        next
      }
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
      n_cohort  = as.integer(m$n_cohort %||% m$n_treated),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  attr(out, "IF_per_cell") <- IF_list
  attr(out, "F_g_per_unit") <- F_g_per_unit
  attr(out, "units")        <- units
  attr(out, "panel_mode")   <- if (rc_mode) "repeated_cross_section" else "panel"
  out
}


# Cluster bootstrap SE: resample units with replacement, recompute the
# long-form ATT(g, t) table on each rep, take SD across reps per cell.
#' @keywords internal
#' @noRd
.cs_bootstrap_se <- function(df, args, att_gt, verbose = TRUE) {
  B <- args$bootstrap_reps
  if (B <= 0L) return(att_gt)
  # Keep the panel a plain data.frame and resolve every column lookup
  # BEFORE indexing into it. `[.data.table` evaluates its `i` expression
  # with the table's COLUMNS in scope, so on a panel carrying a column
  # literally named `d` the local `d` was shadowed by the treatment
  # vector and `d[[args$group]]` became treatment[["g"]], failing with
  # "subscript out of bounds". Every panel whose treatment column is
  # named `d` crashed here -- and only here, since bootstrap_reps = 0
  # skips this function, which is why point estimates were unaffected.
  panel   <- as.data.frame(df, stringsAsFactors = FALSE)
  grp_vec <- panel[[args$group]]
  units   <- unique(grp_vec)
  # Row positions per unit, computed once. This also removes a full-panel
  # scan per (replicate x pick), which was O(B * n_units * nrow).
  rows_by_unit <- split(seq_len(nrow(panel)),
                        factor(grp_vec, levels = units))

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
      key <- as.character(u)
      pcount[key] <- pcount[key] + 1L
      block <- panel[rows_by_unit[[key]], , drop = FALSE]
      if (pcount[key] > 1L) {
        shift <- (pcount[key] - 1L) * OFFSET
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
# ---------------------------------------------------------------------------
# Aggregation of ATT(g, t) into event-time / group / calendar / overall
# summaries, WITH standard errors derived from the influence functions.
#
# This mirrors did::aggte(). The aggregate's influence function is the same
# weighted combination of the cells' influence functions, PLUS a correction
# for the fact that the aggregation weights are themselves estimated
# (did:::wif). Omitting that correction understates the SE badly at long
# event times, where few cohorts contribute: measured against
# did::aggte(type = "dynamic"), fixed-weight SEs ran from 0.98x of the
# correct value at event 0 down to 0.33x at event 9.
#
# Scaling convention, taken from did:::compute.att_gt, which rescales each
# cell's influence onto the FULL sample before aggregating:
#     psi_full(i, c) = (n / n_c) * psi_c(i)   for i in cell c, else 0
# and did:::getSE, which reports se = sqrt(mean(psi^2) / n). For a single
# cell this collapses to sqrt(sum(psi_c^2)) / n_c, exactly the per-cell SE
# that DRDID reports -- so cell and aggregate SEs sit on one footing.
#' @keywords internal
#' @noRd
.cs_agg_se <- function(cell_idx, att_cells, pg_cells, g_cells,
                        IF_list, units, F_g_per_unit) {
  K <- length(cell_idx)
  if (K == 0L) return(NA_real_)
  if (is.null(IF_list) || is.null(units)) return(NA_real_)
  n <- length(units)
  if (n == 0L) return(NA_real_)
  key <- as.character(units)
  pos <- stats::setNames(seq_len(n), key)

  # Cell influence functions, rescaled onto the full sample.
  infl <- matrix(0, nrow = n, ncol = K)
  for (k in seq_len(K)) {
    cl <- IF_list[[cell_idx[k]]]
    if (is.null(cl) || is.null(cl$IF) || !length(cl$IF)) return(NA_real_)
    idx <- pos[as.character(cl$units)]
    if (anyNA(idx)) return(NA_real_)
    nc <- length(cl$IF)
    infl[idx, k] <- (n / nc) * cl$IF
  }

  pg <- as.numeric(pg_cells)
  spg <- sum(pg)
  if (!is.finite(spg) || spg <= 0) return(NA_real_)
  a <- pg / spg

  psi <- as.numeric(infl %*% a)

  # Weight-estimation correction (did:::wif). Vanishes when every cell in
  # the level shares one cohort (group aggregation), because then the
  # numerator and denominator effects cancel exactly.
  if (!is.null(F_g_per_unit) && !anyNA(g_cells)) {
    Gu <- F_g_per_unit[key]
    Gind <- matrix(0, nrow = n, ncol = K)
    for (k in seq_len(K)) {
      Gind[, k] <- as.numeric(!is.na(Gu) & Gu == g_cells[k])
    }
    dev <- sweep(Gind, 2L, pg, "-")
    if1 <- dev / spg
    if2 <- rowSums(dev) %*% t(pg / (spg^2))
    wif <- if1 - if2
    psi <- psi + as.numeric(wif %*% as.matrix(att_cells))
  }

  sqrt(mean(psi^2) / n)
}


# SE for one aggregation level, given the cells it covers.
#' @keywords internal
#' @noRd
.cs_level_se <- function(cells, att_gt, IF_list, units, F_g_per_unit, n_units) {
  idx <- which(cells)
  if (!length(idx)) return(NA_real_)
  .cs_agg_se(cell_idx = idx,
             att_cells = att_gt$att[idx],
             pg_cells  = (if ("n_cohort" %in% names(att_gt)) att_gt$n_cohort
                          else att_gt$n_treated)[idx] / n_units,
             g_cells   = att_gt$g[idx],
             IF_list   = IF_list, units = units,
             F_g_per_unit = F_g_per_unit)
}


.cs_aggregate <- function(att_gt, aggregation, args) {
  if (nrow(att_gt) == 0L) {
    return(data.frame(scheme = aggregation, level = character(0),
                      estimate = numeric(0), se = numeric(0),
                      n_cells = integer(0)))
  }
  # did weights each cell by its cohort's share of units (pg). On a
  # balanced panel that equals the cell's treated count; under the
  # repeated-cross-section path the cell count is in ROWS, so the cohort
  # size is used explicitly.
  w <- if ("n_cohort" %in% names(att_gt)) att_gt$n_cohort else att_gt$n_treated
  # Cells that could not be estimated -- typically a cohort with no
  # treated unit left in period t on an unbalanced panel -- carry
  # att = NA with weight n_treated = 0. In R, NA * 0 is still NA, so a
  # single such cell turned every aggregate containing it into NA: on the
  # trade subsamples of a real application the overall ATT was NA on all
  # twelve. did::aggte either stops on NA cells or, with na.rm = TRUE,
  # removes them before weighting (compute.aggte, lines 52-76); this is
  # the latter. `ok` masks every aggregation below, event times and
  # calendar periods left with no cell are dropped, and a cohort with no
  # estimable post-treatment cell is dropped from the group aggregation,
  # all as did does.
  ok <- !is.na(att_gt$att)
  n_na <- sum(!ok)
  if (n_na > 0L && !isFALSE(args$verbose)) {
    message(sprintf("[didgpu] %d of %d ATT(g,t) cells could not be estimated ",
                    n_na, nrow(att_gt)),
            "(no treated or no control units observed); they are dropped ",
            "from the aggregation, as did::aggte(na.rm = TRUE) does.")
  }
  IF_list <- attr(att_gt, "IF_per_cell")
  units   <- attr(att_gt, "units")
  F_g_per_unit <- attr(att_gt, "F_g_per_unit")
  n_units <- if (is.null(units)) NA_integer_ else length(units)
  z <- stats::qnorm(0.5 + (args$ci_level %||% 95) / 200)
  se_of <- function(cells) {
    if (is.null(IF_list) || is.null(units)) return(NA_real_)
    .cs_level_se(cells, att_gt, IF_list, units, F_g_per_unit, n_units)
  }
  # did::aggte reports an aggregated SE that is numerically zero as NA
  # (compute.aggte: `se[se <= sqrt(.Machine$double.eps) * 10] <- NA`) --
  # e.g. the universal base period's own event time, whose influence
  # function is identically 0.
  na_tiny <- function(se) {
    se[!is.na(se) & se <= sqrt(.Machine$double.eps) * 10] <- NA_real_
    se
  }
  finish <- function(out) {
    out$se <- na_tiny(out$se)
    out$ci_low  <- out$estimate - z * out$se
    out$ci_high <- out$estimate + z * out$se
    out
  }

  .drop_empty <- function(out) {
    if (is.data.frame(out) && "n_cells" %in% names(out) && nrow(out) > 1L) {
      out <- out[out$n_cells > 0L, , drop = FALSE]
      rownames(out) <- NULL
    }
    out
  }

  out <- switch(aggregation,
    "event" = {
      events <- sort(unique(att_gt$event_time))
      finish(data.frame(
        event_time = events,
        estimate = vapply(events, function(e) {
          cells <- att_gt$event_time == e & ok
          if (!any(cells)) return(NA_real_)
          sum(att_gt$att[cells] * w[cells]) / sum(w[cells])
        }, numeric(1)),
        se = vapply(events, function(e) se_of(att_gt$event_time == e & ok),
                     numeric(1)),
        n_cells = vapply(events, function(e) sum(att_gt$event_time == e & ok),
                          integer(1)),
        stringsAsFactors = FALSE
      ))
    },
    "group" = {
      gs <- sort(unique(att_gt$g))
      finish(data.frame(
        g = gs,
        estimate = vapply(gs, function(g) {
          cells <- att_gt$g == g & att_gt$t >= g & ok
          if (!any(cells)) return(NA_real_)
          sum(att_gt$att[cells] * w[cells]) / sum(w[cells])
        }, numeric(1)),
        se = vapply(gs, function(g) se_of(att_gt$g == g & att_gt$t >= g & ok),
                     numeric(1)),
        n_cells = vapply(gs, function(g) sum(att_gt$g == g & att_gt$t >= g & ok),
                          integer(1)),
        stringsAsFactors = FALSE
      ))
    },
    "calendar" = {
      ts <- sort(unique(att_gt$t))
      finish(data.frame(
        t = ts,
        estimate = vapply(ts, function(t) {
          cells <- att_gt$t == t & att_gt$g <= t & ok
          if (!any(cells)) return(NA_real_)
          sum(att_gt$att[cells] * w[cells]) / sum(w[cells])
        }, numeric(1)),
        se = vapply(ts, function(t) se_of(att_gt$t == t & att_gt$g <= t & ok),
                     numeric(1)),
        n_cells = vapply(ts, function(t) sum(att_gt$t == t & att_gt$g <= t & ok),
                          integer(1)),
        stringsAsFactors = FALSE
      ))
    },
    "overall" = {
      post <- att_gt$t >= att_gt$g & ok
      est <- if (any(post)) sum(att_gt$att[post] * w[post]) / sum(w[post])
             else NA_real_
      se  <- na_tiny(se_of(post))
      data.frame(scheme = "overall", estimate = est, se = se,
                  ci_low = est - z * se, ci_high = est + z * se,
                  n_cells = sum(post), stringsAsFactors = FALSE)
    }
  )
  .drop_empty(out)
}



# Placebo joint test: pre-treatment cells should have ATT ~ 0 under
# parallel trends. Returns a list with per-event-time placebo summary
# and a joint chi-square p-value.
#' @keywords internal
#' @noRd
.cs_placebo_test <- function(att_gt, args) {
  # Unestimable (NA) cells are excluded here too, as in .cs_aggregate.
  pre <- att_gt[att_gt$event_time < 0L & !is.na(att_gt$att), , drop = FALSE]
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
