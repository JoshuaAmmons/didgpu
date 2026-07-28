# ============================================================================
# Aggregate committed cells into a result object.
#
# Cell 0 is the point estimate (no resampling). Cells 1..bootstrap_reps are
# the bootstrap iters. The aggregator:
#
#   - Pulls the per-iter effects/placebos vectors into matrices.
#   - SE per coefficient = sd across bootstrap iters.
#   - CI per coefficient = point +/- z * SE  (normal approx).
#   - Joint chi-square p-values for the effects vector and the placebos
#     vector (using the bootstrap empirical covariance).
#
# Output structure mirrors DIDmultiplegtDYN where it can:
#   $coef     : list(b = named numeric, vcov = matrix)
#   $results  : list(N_Effects, N_Placebos, Effects, ATE, Placebos,
#                    p_jointeffects, p_jointplacebo)
#   $args     : the canonical args bundle
#   $manifest : the cell manifest data.frame
#
# Differences from the reference's output:
#   - Standard errors are bootstrap-derived (not the reference's
#     analytical SE). With enough bootstrap reps this converges, but
#     for small bootstrap_reps the SEs WILL differ from the reference.
#     That is intentional: we are checkpointing the bootstrap, which
#     means SEs come out of that bootstrap rather than out of analytic
#     formulas. The Effects[, 1] column (the point estimate) IS expected
#     to match the reference to floating-point precision under the
#     'reference' backend.
#   - Coefficient names are NOT padded with trailing spaces (the
#     reference's "Effect_1    " naming is preserved as a tolerance in
#     compat.R for downstream code that depends on it; the canonical
#     name in didgpu is "Effect_1" with no padding).
# ============================================================================


#' @keywords internal
#' @noRd
.aggregate_to_result <- function(cells, args, panel_hash) {
  if (length(cells) == 0L) {
    stop("No cells available to aggregate. ",
         "Either bootstrap_reps == 0 with no point estimate yet, or the ",
         "checkpoint manifest is empty.")
  }

  # Pull cells in iter order. Convert names from "0","1",... to integers
  # to be sure the point estimate ends up at the right position.
  iters <- sort(as.integer(names(cells)))
  cells <- cells[as.character(iters)]
  has_point <- 0L %in% iters
  if (!has_point) {
    stop("Cannot aggregate: cell b=0 (point estimate) is missing from the ",
         "checkpoint manifest.")
  }
  boot_iters <- setdiff(iters, 0L)

  # Effects matrix : boot_iters x n_effects (point goes in $coef$b).
  e0 <- cells[["0"]]$effects
  n_e <- length(e0)
  p0 <- cells[["0"]]$placebos
  n_p <- length(p0)

  # Drop degenerate bootstrap iterations instead of crashing. Under
  # sparse-switching treatments a resample can contain no valid switcher
  # cell at some horizon, yielding an effects/placebos vector of the wrong
  # length (typically length 0). The vapply() calls below hard-require
  # exact lengths, so a single such iteration used to abort the whole
  # aggregation ("values must be length K ... result is length 0") -- and
  # the failure probability grows with bootstrap_reps, so exactly the
  # large-rep runs users want for final inference were the ones crashing.
  # Dropping failed resamples is standard bootstrap practice; SEs and the
  # bootstrap covariance are computed from the surviving iterations, and a
  # warning reports how many were dropped.
  n_boot_dropped <- 0L
  if (length(boot_iters) > 0L) {
    ok <- vapply(as.character(boot_iters), function(i) {
      ce <- cells[[i]]$effects
      cp <- cells[[i]]$placebos
      length(ce) == n_e && (n_p == 0L || length(cp) == n_p) && !all(is.na(ce))
    }, logical(1))
    if (any(!ok)) {
      n_boot_dropped <- sum(!ok)
      warning(sprintf(
        paste0("didgpu: dropped %d of %d bootstrap iteration(s) whose resample ",
               "produced degenerate cells (no valid switchers at some horizon); ",
               "SEs/covariance use the remaining %d iterations."),
        n_boot_dropped, length(boot_iters), sum(ok)), call. = FALSE)
      boot_iters <- boot_iters[ok]
    }
  }

  e_mat <- if (length(boot_iters) > 0L) {
    m <- vapply(as.character(boot_iters),
                function(i) cells[[i]]$effects,
                numeric(n_e))
    if (is.null(dim(m))) matrix(m, nrow = n_e, ncol = length(boot_iters))
    else m
  } else matrix(numeric(0), nrow = n_e, ncol = 0L)
  e_mat <- t(e_mat)   # iter x effect

  p_mat <- if (length(boot_iters) > 0L && n_p > 0L) {
    m <- vapply(as.character(boot_iters),
                function(i) cells[[i]]$placebos,
                numeric(n_p))
    if (is.null(dim(m))) matrix(m, nrow = n_p, ncol = length(boot_iters))
    else m
  } else matrix(numeric(0), nrow = n_p, ncol = 0L)
  p_mat <- t(p_mat)

  # ---- Match the reference's reported horizon count ----
  # didgpu's data-availability horizon clamp (max L_g per group) can be one
  # step more permissive than the reference's cohort-level T_g clamp, so it may
  # attempt one extra horizon the reference deems infeasible. When that horizon
  # has no switcher reaching it, the point estimate is NA (N_inc == 0) and the
  # reference simply omits the row. Trim the trailing contiguous block of NA
  # point estimates from effects and placebos so the reported horizon count
  # matches did_multiplegt_dyn. Only TRAILING NAs are dropped (a single clamp
  # cutoff, as the reference does); any interior NA is preserved. An all-NA
  # vector is left intact (handled downstream as a zero-effect result).
  .last_estimable <- function(v) { w <- which(!is.na(v)); if (length(w)) max(w) else length(v) }
  ke <- .last_estimable(e0)
  if (ke < n_e) { e0 <- e0[seq_len(ke)]; e_mat <- e_mat[, seq_len(ke), drop = FALSE]; n_e <- ke }
  if (n_p > 0L) {
    # If EVERY placebo is unestimable (all NA point estimates), drop the entire
    # block to match did_multiplegt_dyn (which returns NULL Placebos in that
    # case). Otherwise trim only the trailing contiguous NA block, as for
    # effects. This handles weighted trends_lin / short-panel cases where no
    # group has the F_g - q - 1 pre-period any placebo needs.
    kp <- if (all(is.na(p0))) 0L else .last_estimable(p0)
    if (kp < n_p) { p0 <- p0[seq_len(kp)]; p_mat <- p_mat[, seq_len(kp), drop = FALSE]; n_p <- kp }
  }

  ate_vec <- if (length(boot_iters) > 0L) {
    vapply(as.character(boot_iters),
           function(i) cells[[i]]$ate %||% NA_real_,
           numeric(1))
  } else numeric(0)

  # SEs and CIs. A kept bootstrap iteration can still carry NA at SOME
  # horizons (its resample has switchers overall but none reaching horizon
  # j). Plain sd()/cov() would then return NA for every affected column --
  # with few switchers this wiped out ALL SEs. Compute per-horizon SEs
  # from the finite draws, requiring a minimum bootstrap support of
  # MIN_BOOT_SUPPORT finite draws per horizon (below that the SE is
  # genuinely not estimable and stays NA).
  MIN_BOOT_SUPPORT <- 30L
  .col_sd <- function(m) {
    if (nrow(m) < 2L) return(rep(NA_real_, ncol(m)))
    apply(m, 2L, function(col) {
      v <- col[is.finite(col)]
      if (length(v) >= max(2L, MIN_BOOT_SUPPORT)) stats::sd(v) else NA_real_
    })
  }
  z <- stats::qnorm(0.5 + (args$ci_level %||% 95) / 200)
  e_se <- if (nrow(e_mat) >= 2L) .col_sd(e_mat) else rep(NA_real_, n_e)
  p_se <- if (nrow(p_mat) >= 2L) .col_sd(p_mat) else rep(NA_real_, n_p)
  ate_ok <- ate_vec[is.finite(ate_vec)]
  ate_se <- if (length(ate_ok) >= max(2L, MIN_BOOT_SUPPORT)) stats::sd(ate_ok) else NA_real_

  e_ci_lo <- e0 - z * e_se;  e_ci_hi <- e0 + z * e_se
  p_ci_lo <- p0 - z * p_se;  p_ci_hi <- p0 + z * p_se
  ate0 <- cells[["0"]]$ate %||% NA_real_
  ate_ci_lo <- ate0 - z * ate_se
  ate_ci_hi <- ate0 + z * ate_se

  # NB: paste0("Effect_", seq_len(0)) returns the length-1 string "Effect_"
  # (zero-length recycling), which then mismatches a 0-row Effects matrix and
  # crashes rownames<-. Guard for n_e == 0 exactly as the placebo path does.
  # n_e == 0 arises e.g. with trends_lin on panels where no group has the
  # required F_g-2 pre-period, so no event-study effect is estimable.
  effect_names <- if (n_e > 0L) paste0("Effect_", seq_len(n_e)) else character(0)
  placebo_names <- if (n_p > 0L) paste0("Placebo_", seq_len(n_p)) else character(0)

  # Count columns come from cell b=0 (the point estimate). The reference
  # reports FOUR distinct count columns:
  #   N           = unweighted contributing observations
  #   Switchers   = unweighted switcher cells
  #   N.w         = weighted   contributing observations (sum of weights)
  #   Switchers.w = weighted   switcher cells (sum of switcher weights)
  # The r/reference backends populate dedicated fields for each; older
  # backends (cpu/cuda/fect_*) only run on UNWEIGHTED panels where N_gt is
  # 0/1, so the weighted and unweighted columns coincide -- we fall back to
  # the unweighted field (and to n_inc for the switcher count, which equals
  # the unweighted switcher count there). This keeps unweighted output
  # bit-identical to the pre-weight-fix behavior.
  # Slice each count vector to the (possibly trimmed) horizon count seq_len(n_e)
  # / seq_len(n_p); cell fields still carry the pre-trim length.
  c0 <- cells[["0"]]
  n_eff_e    <- (c0$n_eff_effects %||% rep(NA_integer_, n_e))[seq_len(n_e)]              # N
  n_eff_p    <- (c0$n_eff_placebos %||% rep(NA_integer_, n_p))[seq_len(n_p)]
  n_sw_unw_e <- (c0$n_sw_unw_effects %||% c0$n_inc_effects %||% rep(NA_integer_, n_e))[seq_len(n_e)]   # Switchers
  n_sw_unw_p <- (c0$n_sw_unw_placebos %||% c0$n_inc_placebos %||% rep(NA_integer_, n_p))[seq_len(n_p)]
  n_eff_w_e  <- (c0$n_eff_w_effects %||% n_eff_e)[seq_len(n_e)]                          # N.w
  n_eff_w_p  <- (c0$n_eff_w_placebos %||% n_eff_p)[seq_len(n_p)]
  n_sw_w_e   <- (c0$n_sw_w_effects %||% n_sw_unw_e)[seq_len(n_e)]                        # Switchers.w
  n_sw_w_p   <- (c0$n_sw_w_placebos %||% n_sw_unw_p)[seq_len(n_p)]

  # Effects matrix, shape (n_e x 8) matching DIDmultiplegtDYN.
  Effects <- cbind(
    Estimate = e0, SE = e_se, LB.CI = e_ci_lo, UB.CI = e_ci_hi,
    N = n_eff_e, Switchers = n_sw_unw_e,
    N.w = n_eff_w_e, Switchers.w = n_sw_w_e
  )
  rownames(Effects) <- effect_names

  Placebos <- cbind(
    Estimate = p0, SE = p_se, LB.CI = p_ci_lo, UB.CI = p_ci_hi,
    N = n_eff_p, Switchers = n_sw_unw_p,
    N.w = n_eff_w_p, Switchers.w = n_sw_w_p
  )
  if (n_p > 0L) rownames(Placebos) <- placebo_names

  ate_n   <- if (length(n_sw_unw_e) > 0L) sum(n_sw_unw_e, na.rm = TRUE) else NA_integer_
  ate_n_w <- if (length(n_sw_w_e)   > 0L) sum(n_sw_w_e,   na.rm = TRUE) else NA_integer_
  ATE <- matrix(c(ate0, ate_se, ate_ci_lo, ate_ci_hi,
                  NA_integer_, ate_n, NA_integer_, ate_n_w),
                nrow = 1L,
                dimnames = list("ATE",
                                c("Estimate", "SE", "LB.CI", "UB.CI",
                                  "N", "Switchers", "N.w", "Switchers.w")))

  # Joint chi-square p-values via the bootstrap empirical covariance.
  p_joint_e <- .joint_pvalue(e0, e_mat)
  p_joint_p <- if (n_p > 0L) .joint_pvalue(p0, p_mat) else NA_real_

  # Coefficient vector + bootstrap vcov for the full b parameter
  # (effects then placebos).
  b <- c(e0, p0)
  names(b) <- c(effect_names, placebo_names)
  V <- if (nrow(e_mat) >= 2L) {
    full <- if (n_p > 0L) cbind(e_mat, p_mat) else e_mat
    stats::cov(full, use = "pairwise.complete.obs")
  } else matrix(NA_real_, nrow = length(b), ncol = length(b))
  dimnames(V) <- list(names(b), names(b))

  # predict_het: carry the iter-0 cell's block (a data.frame) through
  # into results$predict_het. If absent, omit the field.
  het_block <- cells[["0"]]$predict_het

  results_list <- list(
    N_Effects      = as.integer(n_e),
    N_Placebos     = as.integer(n_p),
    Effects        = Effects,
    ATE            = ATE,
    Placebos       = Placebos,
    p_jointeffects = p_joint_e,
    p_jointplacebo = p_joint_p,
    n_boot         = length(boot_iters),
    n_boot_dropped = n_boot_dropped
  )
  if (!is.null(het_block)) results_list$predict_het <- het_block

  list(
    coef = list(b = b, vcov = V),
    results = results_list,
    args = c(args, list(panel_hash = panel_hash)),
    cells_used = length(cells)
  )
}


# Joint chi-square p-value: theta0' V^-1 theta0 ~ chi2(k) under H0:
# theta = 0. V is the bootstrap covariance.
#
# Two robustness rules (0.1.2):
#  - Horizons with fewer than 30 finite bootstrap draws are excluded from
#    the joint test (their variance is not estimable), and NA draws at
#    kept horizons are handled with a pairwise-complete covariance.
#  - When the covariance of a high-dimensional coefficient block is
#    near-singular (rcond < 1e-10, common with many horizons and few
#    switchers), the chi-square statistic is numerically unstable: a
#    warning tells the user to prefer a low-dimensional prespecified
#    test (e.g. the leads nearest treatment) over this omnibus p.
.joint_pvalue <- function(theta0, boot_mat) {
  if (nrow(boot_mat) < 2L) return(NA_real_)
  support <- colSums(is.finite(boot_mat))
  keep <- is.finite(theta0) & support >= 30L
  if (!any(keep)) return(NA_real_)
  th <- theta0[keep]
  bm <- boot_mat[, keep, drop = FALSE]
  if (nrow(bm) < length(th) + 1L) return(NA_real_)
  V <- stats::cov(bm, use = "pairwise.complete.obs")
  if (any(!is.finite(V))) return(NA_real_)
  rc <- suppressWarnings(tryCatch(rcond(V), error = function(e) NA_real_))
  if (is.finite(rc) && rc < 1e-10) {
    warning(sprintf(
      paste0("didgpu: the %d-dimensional bootstrap covariance behind a joint ",
             "test is near-singular (rcond = %.1e); the omnibus chi-square ",
             "p-value is numerically unreliable. Prefer a low-dimensional ",
             "prespecified test (e.g. the leads nearest treatment)."),
      length(th), rc), call. = FALSE)
  }
  inv <- try(solve(V), silent = TRUE)
  if (inherits(inv, "try-error")) {
    inv <- MASS::ginv(V)
  }
  q <- as.numeric(t(th) %*% inv %*% th)
  stats::pchisq(q, df = length(th), lower.tail = FALSE)
}


# -------- pretty-print helpers --------

.sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) "***" else
  if (p < 0.01)  "**"  else
  if (p < 0.05)  "*"   else
  if (p < 0.1)   "."   else ""
}

.print_coef_block <- function(m) {
  est <- as.numeric(m[, "Estimate"])
  se  <- as.numeric(m[, "SE"])
  lo  <- as.numeric(m[, "LB.CI"])
  hi  <- as.numeric(m[, "UB.CI"])
  z   <- est / se
  p   <- 2 * stats::pnorm(-abs(z))
  stars <- vapply(p, .sig_stars, character(1))

  df <- data.frame(
    Estimate  = sprintf("%9.4f", est),
    SE        = ifelse(is.na(se), "      NA", sprintf("%8.4f", se)),
    z         = ifelse(is.na(z),  "     NA", sprintf("%7.2f", z)),
    p         = ifelse(is.na(p),  "     NA", sprintf("%7.4f", p)),
    `      CI` = ifelse(is.na(lo) | is.na(hi), "    [NA, NA]",
                        sprintf("[%6.3f, %6.3f]", lo, hi)),
    sig       = stars,
    row.names = rownames(m),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  print(df, right = FALSE)
}


# -------- S3 methods --------

#' Print method for didgpu_result
#'
#' @param x A `didgpu_result` object.
#' @param ... Unused (for S3 method compatibility).
#' @return The input invisibly.
#' @export
print.didgpu_result <- function(x, ...) {
  cat("didgpu result\n")
  cat(sprintf("  backend         : %s\n", x$args$backend %||% "n/a"))
  cat(sprintf("  effects         : %d   placebos: %d\n",
              x$results$N_Effects, x$results$N_Placebos))
  cat(sprintf("  bootstrap reps  : %d (used %d cells)\n",
              x$args$bootstrap_reps, x$cells_used))
  if (!is.null(x$checkpoint_dir) && !is.na(x$checkpoint_dir)) {
    cat(sprintf("  checkpoint_dir  : %s\n", x$checkpoint_dir))
  }
  cat("\nEffects:\n")
  .print_coef_block(x$results$Effects)
  if (x$results$N_Placebos > 0L) {
    cat("\nPlacebos:\n")
    .print_coef_block(x$results$Placebos)
  }
  if (!is.null(x$results$ATE) && !is.na(x$results$ATE[1, "Estimate"])) {
    cat("\nATE (cumulative across event-times):\n")
    .print_coef_block(x$results$ATE)
  }
  cat("\n---\n")
  cat(sprintf("Joint test of effects:  chi2 p = %.4g %s\n",
              x$results$p_jointeffects,
              .sig_stars(x$results$p_jointeffects)))
  if (x$results$N_Placebos > 0L) {
    cat(sprintf("Joint test of placebos: chi2 p = %.4g %s\n",
                x$results$p_jointplacebo,
                .sig_stars(x$results$p_jointplacebo)))
  }
  cat("Signif: *** p<0.001  ** p<0.01  * p<0.05  . p<0.1\n")
  invisible(x)
}

#' Summary method for didgpu_result
#'
#' @param object A `didgpu_result` object.
#' @param ... Unused (for S3 method compatibility).
#' @return The input invisibly.
#' @export
summary.didgpu_result <- function(object, ...) {
  print.didgpu_result(object, ...)
}


#' Coefficient extractor for didgpu_result
#'
#' Returns a named numeric vector of point estimates. Names are
#' "Effect_1", ..., "Effect_n_effects" optionally followed by
#' "Placebo_1", ..., "Placebo_n_placebos" and (if available) "ATE".
#'
#' @param object A `didgpu_result` object.
#' @param which One of `"effects"`, `"placebos"`, `"ate"`, or `"all"`
#'   (default). Controls which subset is returned.
#' @param ... Unused.
#' @return Named numeric vector.
#' @export
coef.didgpu_result <- function(object, which = "all", ...) {
  which <- match.arg(which, c("all", "effects", "placebos", "ate"))
  pieces <- list()
  if (which %in% c("all", "effects")) {
    e <- object$results$Effects
    if (!is.null(e) && nrow(e) > 0L) {
      v <- as.numeric(e[, "Estimate"])
      names(v) <- rownames(e)
      pieces$effects <- v
    }
  }
  if (which %in% c("all", "placebos") && object$results$N_Placebos > 0L) {
    p <- object$results$Placebos
    if (!is.null(p) && nrow(p) > 0L) {
      v <- as.numeric(p[, "Estimate"])
      names(v) <- rownames(p)
      pieces$placebos <- v
    }
  }
  if (which %in% c("all", "ate") && !is.null(object$results$ATE)) {
    a <- object$results$ATE
    if (!is.na(a[1, "Estimate"])) {
      v <- as.numeric(a[, "Estimate"])
      names(v) <- "ATE"
      pieces$ate <- v
    }
  }
  v <- unlist(pieces, use.names = TRUE)
  # Strip the leading list-name prefix ("effects.", "placebos.", "ate.")
  # that unlist() injects when the parent list is named.
  names(v) <- sub("^(effects|placebos|ate)\\.", "", names(v))
  v
}


#' Confidence intervals for didgpu_result
#'
#' Returns the percentile-based CI matrix already computed during
#' aggregation (from the bootstrap distribution). Re-computing at a
#' different level requires a fresh run because the cell-level
#' percentiles are not stored on disk.
#'
#' @param object A `didgpu_result` object.
#' @param parm Optional character vector of coefficient names to subset.
#' @param level Confidence level. Must match the level used at fit time;
#'   otherwise a warning is issued and the stored CI is returned anyway.
#' @param ... Unused.
#' @return A 2-column matrix with the lower and upper bounds.
#' @export
confint.didgpu_result <- function(object, parm = NULL, level = NULL, ...) {
  fit_level <- object$args$ci_level %||% 95
  if (!is.null(level) && abs(level * 100 - fit_level) > 1e-9 &&
      abs(level     - fit_level) > 1e-9) {
    warning(sprintf("Stored CIs are at level %.1f%% (from the fit); ",
                    fit_level),
            "requested level=", level, " is ignored. Re-run didgpu() ",
            "with ci_level = ", level, " to change.")
  }
  # Pick lower/upper-bound columns from the printed tables.
  lo_hi_from <- function(m) {
    if (is.null(m) || nrow(m) == 0L) return(NULL)
    cn <- colnames(m)
    lo <- grep("^LB", cn)[1L]; hi <- grep("^UB", cn)[1L]
    out <- cbind(as.numeric(m[, lo]), as.numeric(m[, hi]))
    rownames(out) <- rownames(m)
    out
  }
  pieces <- list()
  pieces$effects <- lo_hi_from(object$results$Effects)
  if (object$results$N_Placebos > 0L) pieces$placebos <- lo_hi_from(object$results$Placebos)
  if (!is.null(object$results$ATE) &&
      !is.na(object$results$ATE[1, "Estimate"])) {
    pieces$ate <- lo_hi_from(object$results$ATE)
    if (!is.null(pieces$ate)) rownames(pieces$ate) <- "ATE"
  }
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (length(pieces) == 0L) return(matrix(numeric(0), nrow = 0L, ncol = 2L,
                                           dimnames = list(NULL,
                                                            c("LB", "UB"))))
  out <- do.call(rbind, pieces)
  colnames(out) <- c(sprintf("%g%% LB", fit_level),
                     sprintf("%g%% UB", fit_level))
  if (!is.null(parm)) {
    miss <- setdiff(parm, rownames(out))
    if (length(miss)) {
      warning("parm not in result: ", paste(miss, collapse = ", "))
    }
    out <- out[intersect(parm, rownames(out)), , drop = FALSE]
  }
  out
}


#' Event-study plot of a didgpu result
#'
#' Plots estimates against event-time horizon: pre-treatment placebos
#' at negative horizons, post-treatment effects at positive horizons,
#' with the stored CIs as vertical error bars and a horizontal dashed
#' line at 0 for reference. Uses base R graphics — no ggplot2
#' dependency. Returns the input invisibly so calls can be chained.
#'
#' @param x A `didgpu_result` object.
#' @param ... Extra graphical parameters passed to the underlying
#'   `plot()` call (e.g. `main`, `xlab`, `ylab`, `col`, `pch`, `lwd`,
#'   `xlim`, `ylim`).
#' @param show_zero_line Logical. Draw a dashed line at y = 0
#'   (default TRUE).
#' @param show_zero_horizon Logical. Mark the boundary between placebo
#'   and effect horizons with a vertical dashed line (default TRUE).
#' @param ci Logical. Draw error bars at the stored CI level
#'   (default TRUE; suppressed when bootstrap_reps = 0 because the
#'   CIs are NA).
#' @return The input invisibly.
#' @export
plot.didgpu_result <- function(x, ...,
                                show_zero_line = TRUE,
                                show_zero_horizon = TRUE,
                                ci = TRUE) {
  e  <- x$results$Effects
  pl <- if (x$results$N_Placebos > 0L) x$results$Placebos else NULL
  if ((is.null(e) || nrow(e) == 0L) && is.null(pl)) {
    stop("Nothing to plot: result has no effects or placebos.")
  }

  # Build the (horizon, estimate, LB, UB) frame.
  rows_eff <- if (!is.null(e) && nrow(e) > 0L) {
    data.frame(h  = seq_len(nrow(e)),
               y  = as.numeric(e[, "Estimate"]),
               lb = as.numeric(e[, grep("^LB", colnames(e))[1L]]),
               ub = as.numeric(e[, grep("^UB", colnames(e))[1L]]))
  } else NULL
  rows_pl <- if (!is.null(pl) && nrow(pl) > 0L) {
    data.frame(h  = -seq_len(nrow(pl)),
               y  = as.numeric(pl[, "Estimate"]),
               lb = as.numeric(pl[, grep("^LB", colnames(pl))[1L]]),
               ub = as.numeric(pl[, grep("^UB", colnames(pl))[1L]]))
  } else NULL
  df <- rbind(rows_pl, rows_eff)
  df <- df[order(df$h), , drop = FALSE]

  # Suppress CIs if no bootstrap was run.
  if (isTRUE(ci) && (is.null(x$args$bootstrap_reps) ||
                      x$args$bootstrap_reps == 0L)) ci <- FALSE

  # Defaults for graphical params, overridable via ...
  dots <- list(...)
  if (is.null(dots$xlab)) dots$xlab <- "Event-time horizon"
  if (is.null(dots$ylab)) dots$ylab <- "Estimate"
  if (is.null(dots$main)) dots$main <- "didgpu event-study"
  if (is.null(dots$pch))  dots$pch  <- 19L
  if (is.null(dots$col))  dots$col  <- "black"
  if (is.null(dots$xlim)) dots$xlim <- range(df$h) + c(-0.5, 0.5)
  if (is.null(dots$ylim)) {
    yvals <- if (isTRUE(ci)) c(df$y, df$lb, df$ub) else df$y
    yvals <- yvals[is.finite(yvals)]
    if (length(yvals) == 0L) yvals <- c(-1, 1)
    pad <- 0.05 * diff(range(yvals))
    if (pad == 0) pad <- 0.1 * max(1, abs(yvals[1L]))
    dots$ylim <- range(yvals) + c(-pad, pad)
  }

  do.call(plot,
          c(list(df$h, df$y, type = "p"), dots))

  if (isTRUE(show_zero_line)) {
    graphics::abline(h = 0, lty = 2L, col = "grey50")
  }
  if (isTRUE(show_zero_horizon)) {
    graphics::abline(v = 0.5, lty = 2L, col = "grey50")
  }
  if (isTRUE(ci)) {
    graphics::segments(df$h, df$lb, df$h, df$ub,
                       col = dots$col, lwd = max(1L, dots$lwd %||% 1L))
    # Endcaps.
    w <- 0.1
    graphics::segments(df$h - w, df$lb, df$h + w, df$lb,
                       col = dots$col, lwd = max(1L, dots$lwd %||% 1L))
    graphics::segments(df$h - w, df$ub, df$h + w, df$ub,
                       col = dots$col, lwd = max(1L, dots$lwd %||% 1L))
  }

  invisible(x)
}


#' Variance-covariance matrix of estimates
#'
#' Returns the empirical covariance matrix of the bootstrap replicate
#' distribution over `(Effects, Placebos)`, computed at fit time and
#' stored on the result. When `bootstrap_reps = 0` (or only one rep),
#' returns a square NA matrix because the covariance is undefined.
#'
#' The ordering matches `coef(object)`'s default (effects first, then
#' placebos). The ATE row/column is NOT included — it is a linear
#' combination of the per-event-time effects, so its variance can be
#' recovered as `t(w) %*% vcov(object) %*% w` where `w` is the
#' incidence-weighted vector.
#'
#' @param object A `didgpu_result` object.
#' @param ... Unused.
#' @return A square matrix.
#' @export
vcov.didgpu_result <- function(object, ...) {
  V <- object$coef$vcov
  if (is.null(V)) {
    nm <- names(object$coef$b)
    return(matrix(NA_real_, length(nm), length(nm),
                  dimnames = list(nm, nm)))
  }
  V
}
