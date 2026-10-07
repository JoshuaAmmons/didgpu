# ============================================================================
# Controls support — WORK IN PROGRESS, not wired into the backend.
#
# This file holds the pre-fit logic that computes the per-baseline-level
# regression coefficients (theta_d) from the FWL machinery the reference
# uses. To complete:
#
# 1. (Done in this file) `.prefit_controls(prepped, controls)` returns
#    theta_d[l, k] from the OLS of sqrt(N_gt)*diff_y on
#    sqrt(N_gt)*(diff_X - avg_diff_X[time, d_sq]) restricted to
#    never-switcher rows of baseline level l.
#
# 2. (Not done) `.core_one_event_time` must accept the prefit and, for
#    each event-time k, build `diff_X_k = X - shift(X, k)` per control,
#    then adjust `diff_y_k := diff_y_k - sum_c theta_d[l, c] * diff_X_c_k`
#    for rows with d_sq == l. Reference: core.R:267-269. The current
#    `.apply_controls_adjustment` only handles k = 1 (it modifies
#    diff_y_XX in place, which only flows through to the k=1 case).
#
# 3. (Not done) Add `controls` arg to `didgpu()`, plumb through
#    `args$controls` to `.backend_r_impl`, call `.prefit_controls` after
#    `.prep_panel`, pass the prefit to the per-k kernel calls.
#
# 4. (Not done) Variance correction `part2_switch` (core.R:504-533).
#    Affects SE only; bootstrap SEs work without it.
#
# Algorithm (verified against did_multiplegt_dyn_core.R:225-273 and
# did_multiplegt_main.R:348-491):
#
#   For each control variable c in 1..K:
#     diff_X_c    = X_c - lag(X_c)         by group
#     diff_X_c_N  = sqrt(N_gt) * diff_X_c
#
#   For each (time, d_sq):
#     avg_diff_X_c[time, d_sq] = mean(diff_X_c restricted to
#                                       ever_change_d == 0 & non-missing)
#     resid_X_c = sqrt(N_gt) * (diff_X_c - avg_diff_X_c[time, d_sq])
#
#   For each baseline level l:
#     Restrict to ever_change_d == 0 & d_sq == l & non-missing rows.
#     X_mat = [resid_X_1, ..., resid_X_K]
#     y     = sqrt(N_gt) * diff_y
#     theta_d[l, ] = ginv(X_mat' X_mat) X_mat' y
#
#   Inside the per-event-time core, BEFORE building diff_y_k:
#     For each row with d_sq == l and each control c, compute
#       diff_X_c_k = X_c - shift(X_c, k) by group
#     and subtract sum_c theta_d[l, c] * diff_X_c_k from diff_y_k.
#
# Variance correction (part2_switch) is deferred — bootstrap SEs work fine
# without it.
# ============================================================================


#' Compute pre-fit controls coefficients per baseline level
#'
#' @param prepped Output of `.prep_panel()`.
#' @param controls Character vector of column names in prepped.
#' @return A list with `theta` (named list keyed by baseline level,
#'   each a numeric vector of length K) and `diff_X` (a list of
#'   the per-row first-differenced control columns, length K, each
#'   length n_rows; the kernel uses these to adjust diff_y).
#' @keywords internal
#' @noRd
.prefit_controls <- function(prepped, controls) {
  d <- prepped
  K <- length(controls)
  if (K == 0L) return(list(theta = list(), diff_X = list(), controls = controls))

  # ever_change_d_XX is PER-(group, time): 1 if this row is at or after
  # the group's first switch, 0 before. So a switcher unit has
  # ever_change=0 for its pre-switch rows AND those rows count as
  # "controls" for the regression below. Reference: main.R:193-196.
  d[, ever_change_d_XX := as.integer(time_XX >= F_g_XX)]
  d[is.na(F_g_XX), ever_change_d_XX := 0L]

  # First-difference each control by group.
  diff_X <- vector("list", K)
  for (k in seq_len(K)) {
    cname <- controls[k]
    if (!cname %in% names(d))
      stop("control variable not in prepped panel: ", cname)
    d[, paste0("diff_X_", k, "_XX") := get(cname) -
        data.table::shift(get(cname), 1L, type = "lag"),
      by = group_XX]
    diff_X[[k]] <- d[[paste0("diff_X_", k, "_XX")]]
  }

  # `fd_X_all_non_missing_XX`: 1 iff every control's diff_X is non-NA.
  # The reference uses this as part of the row mask for the control
  # regression and the per-(time, d_sq) means (main.R:350-355, 367-380).
  d[, fd_X_all_non_missing_XX := 1L]
  for (k in seq_len(K)) {
    d[is.na(get(paste0("diff_X_", k, "_XX"))), fd_X_all_non_missing_XX := 0L]
  }

  # Per (time, d_sq[, trends_nonparam]) mean of diff_X_k among control
  # units. Mask matches the reference's mask: ever_change == 0 AND
  # non-NA diff_y AND all control diffs present.
  cohort_cols <- c("time_XX", "d_sq_XX")
  if ("trends_np_XX" %in% names(d)) cohort_cols <- c(cohort_cols, "trends_np_XX")
  resid_X <- vector("list", K)
  mask_ctrl <- d$ever_change_d_XX == 0L &
               !is.na(d$diff_y_XX) &
               d$fd_X_all_non_missing_XX == 1L &
               d$N_gt_XX > 0
  for (k in seq_len(K)) {
    diff_col   <- paste0("diff_X_", k, "_XX")
    avg_col    <- paste0("avg_diff_X_", k, "_XX")
    resid_col  <- paste0("resid_X_", k, "_XX")

    d[, "._num_XX" := ifelse(mask_ctrl, N_gt_XX * get(diff_col), 0)]
    d[, "._den_XX" := ifelse(mask_ctrl, N_gt_XX, 0)]
    d[, "._sumnum_XX" := sum(._num_XX, na.rm = TRUE),
      by = cohort_cols]
    d[, "._sumden_XX" := sum(._den_XX, na.rm = TRUE),
      by = cohort_cols]
    d[, (avg_col) := ifelse(._sumden_XX > 0, ._sumnum_XX / ._sumden_XX, NA_real_)]
    d[, (resid_col) := sqrt(N_gt_XX) * (get(diff_col) - get(avg_col))]
    d[is.na(get(resid_col)), (resid_col) := 0]
    d[, c("._num_XX", "._den_XX", "._sumnum_XX", "._sumden_XX") := NULL]

    resid_X[[k]] <- d[[resid_col]]
  }

  # Per baseline level: OLS of diff_y_w on resid_X among controls.
  diff_y_w <- sqrt(d$N_gt_XX) * d$diff_y_XX
  diff_y_w[is.na(diff_y_w)] <- 0

  dsq_levels <- sort(unique(d$d_sq_XX))
  theta <- list()
  for (l in dsq_levels) {
    mask <- d$ever_change_d_XX == 0L &
            d$d_sq_XX == l &
            !is.na(d$diff_y_XX) &
            d$fd_X_all_non_missing_XX == 1L &
            d$N_gt_XX > 0
    n_obs <- sum(mask)
    if (n_obs <= K + 1L) {
      # Not enough rows; theta_l = 0 (no adjustment).
      theta[[as.character(l)]] <- rep(0, K)
      next
    }
    X_mat <- do.call(cbind, lapply(seq_len(K),
                                    function(k) d[[paste0("resid_X_", k, "_XX")]][mask]))
    y_vec <- diff_y_w[mask]
    XtX <- crossprod(X_mat)
    XtY <- crossprod(X_mat, y_vec)
    theta_l <- tryCatch(
      as.numeric(MASS::ginv(XtX) %*% XtY),
      error = function(e) rep(0, K)
    )
    theta[[as.character(l)]] <- theta_l
  }

  list(theta = theta, diff_X = diff_X, controls = controls,
       dsq_levels = dsq_levels)
}


#' Apply the controls FWL adjustment to diff_y_XX in-place
#'
#' Modifies `d$diff_y_XX` by subtracting `sum_k theta_d[l, k] * diff_X_k`
#' for each row with `d_sq == l`. The kernel then uses the adjusted
#' diff_y as if no controls existed.
#'
#' @keywords internal
#' @noRd
.apply_controls_adjustment <- function(d, prefit) {
  K <- length(prefit$controls)
  if (K == 0L) return(d)
  d[, diff_y_XX_orig := diff_y_XX]
  for (l in prefit$dsq_levels) {
    key <- as.character(l)
    theta_l <- prefit$theta[[key]]
    if (is.null(theta_l) || all(theta_l == 0)) next
    mask <- d$d_sq_XX == l
    adj <- rep(0, nrow(d))
    for (k in seq_len(K)) {
      dx <- d[[paste0("diff_X_", k, "_XX")]]
      adj <- adj + theta_l[k] * ifelse(is.na(dx), 0, dx)
    }
    d[mask, diff_y_XX := diff_y_XX - adj[mask]]
  }
  d
}


# The note fixest prints for each per-baseline control regression that
# DIDmultiplegtDYN fits (did_multiplegt_main.R:540-560): one `feols` per
# baseline level with residualised controls, on the rows before the
# group's switch, dropping rows whose outcome or control first
# difference is missing. didgpu does not call fixest, so it says what
# fixest says, in fixest's words. Needs the columns .prefit_controls
# leaves on `prepped`.
#' @keywords internal
#' @noRd
.controls_na_notes <- function(prepped, controls) {
  K <- length(controls)
  if (K == 0L) return(character(0))
  d <- prepped
  dx_cols <- paste0("diff_X_", seq_len(K), "_XX")
  if (!all(c(dx_cols, "fd_X_all_non_missing_XX", "ever_change_d_XX") %in% names(d))) {
    return(character(0))
  }
  fmt <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
  notes <- character(0)
  for (l in sort(unique(d$d_sq_XX[!is.na(d$d_sq_XX)]))) {
    in_l <- !is.na(d$d_sq_XX) & d$d_sq_XX == l
    # Levels the reference residualises: more than one switch date, and
    # a non-empty control sample (main.R:455-468).
    if (length(unique(d$F_g_XX[in_l])) <= 1L) next
    ctl <- in_l & d$ever_change_d_XX == 0L & !is.na(d$diff_y_XX) &
           d$fd_X_all_non_missing_XX == 1L
    if (!any(ctl)) next
    reg <- in_l & d$F_g_XX > d$time_XX
    if (!any(reg)) next
    lhs <- is.na(d$diff_y_XX[reg])
    rhs <- Reduce(`|`, lapply(dx_cols, function(cn) is.na(d[[cn]][reg])))
    # The balancing fill-in rows have no weight in the reference.
    wts <- if ("no_wt_XX" %in% names(d)) d$no_wt_XX[reg] else logical(sum(reg))
    n <- sum(lhs | rhs | wts)
    if (n == 0L) next
    parts <- c(if (any(lhs)) paste0("LHS: ", fmt(sum(lhs))),
               if (any(rhs)) paste0("RHS: ", fmt(sum(rhs))),
               if (any(wts)) paste0("Weights: ", fmt(sum(wts))))
    notes <- c(notes, sprintf("NOTE: %s observation%s removed because of NA values (%s).",
                              fmt(n), if (n == 1L) "" else "s",
                              paste(parts, collapse = ", ")))
  }
  notes
}


# -------- analytic SEs with controls --------
#
# With controls the reference subtracts, from each group's influence
# contribution, a term for the estimation error in the per-baseline
# control coefficients (did_multiplegt_dyn_core.R:254-347 and 544-590 for
# effects, compute_placebo_effects_polars.R:86-120 and 250-300 for
# placebos):
#
#   part2_g = sum_l sum_j M_{l,j} * ( sum_c invDenom_l[j, c] * in_sum_{c,l}(g)
#                                       * 1[d_sq_g = l, F_g >= 3]  -  theta_l[j] )
#
# for every baseline level l the reference residualises. in_sum and the
# bracket do not depend on the horizon, so .controls_se_prep builds them
# once per fit; .controls_part2 adds the horizon's M_{l,j} and returns
# part2 for one (horizon, direction).

# Stata's invsym, as the reference ports it (invsym_r): a pivoted
# Cholesky inverse, zero on the rows and columns of a singular pivot.
.invsym <- function(M) {
  n <- nrow(M)
  ch <- tryCatch(suppressWarnings(chol(M, pivot = TRUE)), error = function(e) NULL)
  if (is.null(ch)) {
    return(tryCatch(solve(M), error = function(e) MASS::ginv(M)))
  }
  piv <- attr(ch, "pivot"); rank <- attr(ch, "rank") %||% n
  if (rank == n) {
    oo <- order(piv)
    return(chol2inv(ch)[oo, oo, drop = FALSE])
  }
  inv <- matrix(0, n, n)
  if (rank > 0) {
    idx <- piv[seq_len(rank)]
    inv[idx, idx] <- chol2inv(ch[seq_len(rank), seq_len(rank), drop = FALSE])
  }
  inv
}

# Fitted values of a weighted regression of y on X (no intercept) with
# period fixed effects, as `fixest::feols(y ~ X - 1 | t, weights = w)`
# then `predict()` give them: NA where X is incomplete or the period has
# no estimation rows. Collinear columns get coefficient 0, which is how
# fixest's dropped variables enter its predictions.
.time_fe_predict <- function(y, X, t, w) {
  X <- as.matrix(X)
  cx <- stats::complete.cases(X)
  est <- cx & !is.na(y) & !is.na(w) & w > 0
  pred <- rep(NA_real_, length(y))
  if (!any(est)) return(pred)
  te <- t[est]; we <- w[est]
  dm <- function(v) v - stats::ave(v * we, te, FUN = sum) / stats::ave(we, te, FUN = sum)
  yd <- dm(y[est])
  Xd <- apply(X[est, , drop = FALSE], 2L, dm)
  if (is.null(dim(Xd))) Xd <- matrix(Xd, ncol = ncol(X))
  fit <- stats::lm.wfit(Xd, yd, we)
  b <- fit$coefficients; b[is.na(b)] <- 0
  res <- y[est] - as.numeric(X[est, , drop = FALSE] %*% b)
  fe <- tapply(res * we, te, sum) / tapply(we, te, sum)
  ok <- cx & as.character(t) %in% names(fe)
  pred[ok] <- as.numeric(X[ok, , drop = FALSE] %*% b) + fe[as.character(t[ok])]
  pred
}

#' @keywords internal
#' @noRd
.controls_se_prep <- function(prepped, prefit) {
  d <- prepped
  K <- length(prefit$controls)
  if (K == 0L) return(list())
  dx_cols <- paste0("diff_X_", seq_len(K), "_XX")
  rx_cols <- paste0("resid_X_", seq_len(K), "_XX")
  if (!all(c(dx_cols, rx_cols, "ever_change_d_XX", "fd_X_all_non_missing_XX") %in% names(d))) {
    return(list())
  }
  G <- length(unique(d$group_XX))
  gl <- d[, list(dsq = d_sq_XX[1L], Fg = F_g_XX[1L]), by = group_XX]
  data.table::setorder(gl, group_XX)
  w_reg <- if ("no_wt_XX" %in% names(d)) ifelse(d$no_wt_XX, NA_real_, d$N_gt_XX) else d$N_gt_XX
  out <- list()
  for (l in sort(unique(d$d_sq_XX[!is.na(d$d_sq_XX)]))) {
    in_l <- !is.na(d$d_sq_XX) & d$d_sq_XX == l
    # The levels the reference residualises (main.R:455-497).
    if (length(unique(d$F_g_XX[in_l])) <= 1L) next
    ctl <- in_l & d$ever_change_d_XX == 0L & !is.na(d$diff_y_XX) &
           d$fd_X_all_non_missing_XX == 1L
    if (!any(ctl)) next
    X <- as.matrix(d[ctl, rx_cols, with = FALSE])
    yw <- sqrt(d$N_gt_XX[ctl]) * d$diff_y_XX[ctl]
    overall <- crossprod(cbind(yw, X, 1))
    if (is.na(sum(overall))) next
    Mx <- overall[1L + seq_len(K), 1L + seq_len(K), drop = FALSE]
    inv_M <- .invsym(Mx)
    theta <- as.numeric(inv_M %*% overall[1L + seq_len(K), 1L, drop = FALSE])
    rmax <- max(d$F_g_XX[ctl])
    rs <- ctl & d$time_XX >= 2L & d$time_XX <= rmax - 1L &
          d$time_XX < d$F_g_XX & !is.na(d$diff_y_XX)
    inv_Denom <- inv_M * sum(d$N_gt_XX[rs]) * G

    # E_y_hat: the period-FE regression of diff_y on the controls' first
    # differences among the level's not-yet-switched rows (main.R:540-560),
    # kept where at least two such rows share the period.
    reg <- in_l & d$F_g_XX > d$time_XX
    eyi <- rep(NA_real_, nrow(d))
    eyi[reg] <- .time_fe_predict(d$diff_y_XX[reg],
                                 as.matrix(d[reg, dx_cols, with = FALSE]),
                                 d$time_XX[reg], w_reg[reg])
    dummy <- as.numeric(reg & !is.na(d$diff_y_XX))
    denom <- stats::ave(dummy, d$time_XX, as.numeric(in_l), FUN = sum)
    denom[!in_l] <- NA_real_
    E_y_hat <- eyi * as.numeric(denom >= 2)
    T_d <- max(d$F_g_XX[in_l]) - 1L
    N_c <- sum(d$N_gt_XX[in_l & d$time_XX >= 2L & d$time_XX <= T_d &
                         d$time_XX < d$F_g_XX & !is.na(d$diff_y_XX)])
    adj <- ifelse(!is.na(denom) & denom > 1, sqrt(denom / (denom - 1)) - 1, 0)
    win <- as.numeric(d$time_XX >= 2L & d$time_XX <= d$F_g_XX - 1L)
    in_sum <- matrix(0, G, K)
    for (c in seq_len(K)) {
      prod_c <- sqrt(d$N_gt_XX) * d[[rx_cols[c]]]
      prod_c[is.na(prod_c)] <- 0
      tmp <- prod_c * (1 + as.numeric(denom >= 2) * adj) *
             (d$diff_y_XX - E_y_hat) * win / N_c
      s <- rowsum(tmp, d$group_XX, na.rm = TRUE, reorder = TRUE)
      in_sum[match(as.integer(rownames(s)), gl$group_XX), c] <- s[, 1L]
    }
    sel <- as.numeric(!is.na(gl$dsq) & gl$dsq == l & gl$Fg >= 3L)
    inB <- (in_sum %*% t(inv_Denom)) * sel
    inB <- sweep(inB, 2L, theta, "-")
    out[[length(out) + 1L]] <- list(l = l, inB = inB)
  }
  out
}

# part2 for one horizon k and one switching direction, from the switcher
# mask the kernel just built. `cols` names that mask's columns; placebo
# horizons difference the controls between F_g-1-k and F_g-1-2k.
#' @keywords internal
#' @noRd
.controls_part2 <- function(d, se_prep, controls, k, G, N_inc,
                            dist, ratio, never, placebo = FALSE) {
  part2 <- numeric(G)
  if (!length(se_prep) || !is.finite(N_inc) || N_inc == 0) return(part2)
  k <- as.integer(k)
  base <- (G / N_inc) * (d[[dist]] - d[[ratio]] * d[[never]]) *
          as.numeric(k <= d$T_g_XX - 2L) *
          as.numeric(d$time_XX >= k + 1L & d$time_XX <= d$T_g_XX) * d$N_gt_XX
  dx <- lapply(controls, function(cn) {
    x <- d[[cn]]
    g <- d$group_XX
    lag <- function(v, n) {
      out <- data.table::shift(v, n, type = "lag")
      out[data.table::shift(g, n, type = "lag") != g] <- NA
      out
    }
    if (placebo) lag(x, 2L * k) - lag(x, k) else x - lag(x, k)
  })
  for (s in se_prep) {
    in_l <- !is.na(d$d_sq_XX) & d$d_sq_XX == s$l
    for (j in seq_along(controls)) {
      M_lj <- sum(base * in_l * dx[[j]], na.rm = TRUE) / G
      part2 <- part2 + M_lj * s$inB[, j]
    }
  }
  part2
}
