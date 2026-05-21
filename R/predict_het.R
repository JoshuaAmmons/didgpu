# ============================================================================
# predict_het: regress per-group ATE contribution on time-invariant
# covariates. For each requested event-time `i`, fits
#
#   prod_het_i = S_g_het * (Y_{F_g + i - 1} - Y_{F_g - 1})
#                ~ covariate(s) + interaction(F_g, d_sq, S_g [, trends_nonparam])
#
# weighted by the user's weight column, with HC1 (heteroskedasticity-
# robust) standard errors on the covariate coefficients and a joint
# F-test across all covariates.
#
# One row in the returned data frame per (event-time x covariate).
# Reference: did_multiplegt_main.R:1641-1745.
# ============================================================================


#' Compute predict_het block for one fit
#'
#' @param prepped The data.table coming out of .prep_panel.
#' @param het_vars Character vector of time-invariant covariate column names.
#' @param het_effects Integer vector of event-times to do het regression for.
#' @param l_eff Maximum feasible event-time.
#' @param weight_col NULL or character; user's weight column (post-prep
#'   it lives as `N_gt_XX`).
#' @param trends_nonparam_col NULL or character; cohort-extension column.
#' @param ci_level Numeric in (0, 100); confidence level for CIs.
#'
#' @return A data.frame with columns
#'   `effect, covariate, Estimate, SE, t, LB, UB, N, pF`.
#'
#' @keywords internal
#' @noRd
.compute_predict_het <- function(prepped, het_vars, het_effects,
                                   l_eff,
                                   trends_nonparam_col = NULL,
                                   ci_level = 95) {
  if (length(het_vars) == 0L) {
    return(data.frame(effect = integer(0), covariate = character(0),
                      Estimate = numeric(0), SE = numeric(0),
                      t = numeric(0), LB = numeric(0), UB = numeric(0),
                      N = integer(0), pF = numeric(0),
                      stringsAsFactors = FALSE))
  }
  d <- prepped

  # Resolve het_effects: -1 means "all 1..l_eff".
  if (any(het_effects == -1L)) {
    het_effects <- seq_len(l_eff)
  } else {
    het_effects <- sort(unique(as.integer(het_effects)))
    bad <- het_effects[het_effects < 1L | het_effects > l_eff]
    if (length(bad)) {
      stop("predict_het: requested event-times out of range 1..", l_eff,
           ": ", paste(bad, collapse = ", "))
    }
  }

  # Validate that het_vars are time-invariant per group (reference checks
  # this at main.R:99-110 and warns otherwise).
  for (v in het_vars) {
    if (!v %in% names(d)) {
      stop("predict_het: covariate '", v, "' not in prepped data.")
    }
    sd_by_g <- d[, list(s = stats::sd(get(v), na.rm = TRUE)),
                 by = group_XX]
    sd_by_g$s[is.na(sd_by_g$s)] <- 0
    if (mean(sd_by_g$s) > 0) {
      warning(sprintf(
        "predict_het: variable '%s' is time-varying within group; ",
        v),
        "using its per-group mean. Reference would drop the variable; ",
        "we keep it for visibility.")
    }
  }

  # Y at F_g - 1 (the immediately pre-switch outcome). NA if F_g == 1
  # (no pre-switch obs) or if the unit is unobserved at F_g - 1.
  d[, Yg_Fg_min1_XX := ifelse(time_XX == F_g_XX - 1L, outcome_XX,
                                NA_real_)]
  d[, Yg_Fg_min1_XX := mean(Yg_Fg_min1_XX, na.rm = TRUE),
    by = group_XX]
  d[, Yg_Fg_min1_XX := ifelse(is.nan(Yg_Fg_min1_XX), NA_real_,
                                Yg_Fg_min1_XX)]
  d[, feasible_het_XX := !is.na(Yg_Fg_min1_XX)]

  # S_g_het: +1 for in-switchers, -1 for out-switchers, NA for never.
  d[, S_g_het_XX := ifelse(is.na(S_g_XX), NA_integer_,
                             ifelse(S_g_XX == 0L, -1L, 1L))]

  # Per-row weight (default 1; weight column was folded into N_gt_XX).
  # Use the raw N_gt_XX (which is 0 for rows that should be excluded).
  d[, weight_XX := N_gt_XX]

  results <- vector("list", length(het_effects))
  ci_z <- stats::qnorm(0.5 + ci_level / 200)  # for normal-based CI fallback

  for (idx in seq_along(het_effects)) {
    i <- het_effects[idx]

    # Per-group outcome at F_g + i - 1 (the i-th post-switch period).
    d[, Yg_Fg_i_XX := ifelse(time_XX == F_g_XX - 1L + i, outcome_XX,
                               NA_real_)]
    d[, Yg_Fg_i_XX := mean(Yg_Fg_i_XX, na.rm = TRUE),
      by = group_XX]
    d[, Yg_Fg_i_XX := ifelse(is.nan(Yg_Fg_i_XX), NA_real_, Yg_Fg_i_XX)]

    # Signed outcome difference: positive contribution if Y went up for
    # in-switchers, also positive if Y went down for out-switchers.
    d[, prod_het_XX := S_g_het_XX * (Yg_Fg_i_XX - Yg_Fg_min1_XX)]

    # Feasible sample: groups that reach event-time i and have valid
    # Y at F_g - 1. Take one row per group (the first one in time order).
    data.table::setorder(d, group_XX, time_XX)
    d[, gr_id_XX := seq_len(.N), by = group_XX]
    sample <- d[gr_id_XX == 1L & feasible_het_XX &
                  (F_g_XX - 1L + i) <= T_g_XX &
                  !is.na(prod_het_XX), ]

    # Drop rows whose covariates are NA so lm doesn't complain.
    keep <- Reduce("&",
                   lapply(het_vars,
                          function(v) !is.na(sample[[v]])))
    sample <- sample[keep, ]
    N_sample <- nrow(sample)

    if (N_sample < length(het_vars) + 2L) {
      # Too few observations to fit anything meaningful.
      for (v in het_vars) {
        results[[idx]] <- rbind(results[[idx]], data.frame(
          effect = i, covariate = v,
          Estimate = NA_real_, SE = NA_real_, t = NA_real_,
          LB = NA_real_, UB = NA_real_,
          N = as.integer(N_sample), pF = NA_real_,
          stringsAsFactors = FALSE))
      }
      next
    }

    # Build the formula. Include cohort-interaction factors that have
    # >1 level, parametrised as full main + interaction (lm drops
    # collinear columns automatically).
    factor_vars <- c("F_g_XX", "d_sq_XX", "S_g_XX")
    if (!is.null(trends_nonparam_col)) {
      factor_vars <- c(factor_vars, trends_nonparam_col)
    }
    factor_vars <- factor_vars[vapply(factor_vars, function(fv) {
      length(unique(sample[[fv]])) > 1L
    }, logical(1))]
    if (length(factor_vars) > 0L) {
      cohort_term <- paste(sprintf("factor(%s)", factor_vars),
                            collapse = ":")
    } else {
      cohort_term <- NULL
    }
    rhs <- paste(c(het_vars, cohort_term), collapse = " + ")
    fml <- stats::as.formula(paste("prod_het_XX ~", rhs))

    fit <- tryCatch(
      stats::lm(fml, data = sample, weights = sample$weight_XX),
      error = function(e) NULL)
    if (is.null(fit)) {
      for (v in het_vars) {
        results[[idx]] <- rbind(results[[idx]], data.frame(
          effect = i, covariate = v,
          Estimate = NA_real_, SE = NA_real_, t = NA_real_,
          LB = NA_real_, UB = NA_real_,
          N = as.integer(N_sample), pF = NA_real_,
          stringsAsFactors = FALSE))
      }
      next
    }

    # HC1 robust covariance: V = (n / (n - k)) * (X'X)^-1 * sum(e_i^2 * X_i X_i') * (X'X)^-1
    # for weighted regression we use the weighted X and residuals.
    bread <- .hc1_vcov(fit)
    coefs <- stats::coef(fit)
    # Identify rows in beta corresponding to het_vars.
    var_pos <- match(het_vars, names(coefs))
    if (any(is.na(var_pos))) {
      # A het_var got dropped by lm (collinear with cohort interactions).
      # Report NA for it.
      for (j in seq_along(het_vars)) {
        v <- het_vars[j]
        pos <- var_pos[j]
        if (is.na(pos)) {
          results[[idx]] <- rbind(results[[idx]], data.frame(
            effect = i, covariate = v,
            Estimate = NA_real_, SE = NA_real_, t = NA_real_,
            LB = NA_real_, UB = NA_real_,
            N = as.integer(N_sample), pF = NA_real_,
            stringsAsFactors = FALSE))
        }
      }
      var_pos <- var_pos[!is.na(var_pos)]
      het_vars_kept <- het_vars[!is.na(match(het_vars, names(coefs)))]
    } else {
      het_vars_kept <- het_vars
    }
    if (length(var_pos) == 0L) next

    est <- unname(coefs[var_pos])
    se  <- sqrt(diag(bread)[var_pos])
    tval <- est / se
    df_resid <- stats::df.residual(fit)
    # t-distribution critical value at ci_level (matches reference: qt at
    # 0.975 for 95% CI).
    tcrit <- stats::qt(0.5 + ci_level / 200, df = df_resid)
    lb <- est - tcrit * se
    ub <- est + tcrit * se

    # Joint F-test on the het_vars block: H0: beta[var_pos] == 0.
    # Wald = (b_sub)' V_sub^{-1} b_sub, F = Wald / q.
    b_sub <- est  # = unname(coefs[var_pos]) (length q)
    V_sub <- bread[var_pos, var_pos, drop = FALSE]
    q <- length(var_pos)
    wald <- tryCatch(
      as.numeric(matrix(b_sub, nrow = 1) %*% solve(V_sub) %*% matrix(b_sub, ncol = 1)),
      error = function(e) NA_real_)
    f_stat <- if (!is.na(wald)) wald / q else NA_real_
    pF <- if (is.na(f_stat)) NA_real_
          else stats::pf(f_stat, q, df_resid, lower.tail = FALSE)

    results[[idx]] <- rbind(results[[idx]], data.frame(
      effect    = i,
      covariate = het_vars_kept,
      Estimate  = est,
      SE        = se,
      t         = tval,
      LB        = lb,
      UB        = ub,
      N         = as.integer(N_sample),
      pF        = pF,
      stringsAsFactors = FALSE))
  }

  # Clean up scratch columns.
  d[, c("Yg_Fg_min1_XX", "feasible_het_XX", "S_g_het_XX",
        "Yg_Fg_i_XX", "prod_het_XX", "gr_id_XX",
        "weight_XX") := NULL]

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}


# HC1 robust variance for a weighted lm fit.
# Reference formula: V = (n / (n - k)) * (X'WX)^-1 * X' diag(w^2 e^2) X * (X'WX)^-1
# This matches sandwich::vcovHC(., type = "HC1") for weighted lm.
#'
#' Extracts X, residuals, and weights FROM THE FITTED MODEL so they all
#' have the same length (post-NA-omission). Drops X columns whose
#' coefficient is NA (collinear / dropped by lm).
#'
#' @keywords internal
#' @noRd
.hc1_vcov <- function(fit) {
  X <- stats::model.matrix(fit)
  # Drop columns lm dropped via aliasing (their coef is NA).
  coef_full <- stats::coef(fit)
  alive <- !is.na(coef_full)
  X <- X[, alive, drop = FALSE]
  n <- nrow(X); k <- ncol(X)
  resid <- stats::residuals(fit)
  w <- stats::weights(fit)
  if (is.null(w)) w <- rep(1, n)
  if (length(w) != n) {
    # Defensive — shouldn't happen because residuals() and weights()
    # both come from the same fitted object.
    stop(".hc1_vcov: weights length (", length(w),
         ") != model.matrix rows (", n, ")")
  }
  # Weighted X'X inverse (the "bread"). With lm weights, this is solve(X' W X).
  XtWX <- crossprod(X, w * X)
  XtWX_inv <- tryCatch(solve(XtWX),
                        error = function(e) MASS::ginv(XtWX))
  # Meat: sum_i w_i^2 e_i^2 X_i X_i' = X' diag(w^2 e^2) X.
  meat <- crossprod(X, (w^2 * resid^2) * X)
  V <- XtWX_inv %*% meat %*% XtWX_inv
  # HC1 small-sample correction: n / (n - k).
  V <- V * (n / max(1, n - k))
  # Re-expand to the full coef vector with NA for dropped columns.
  V_full <- matrix(NA_real_, length(coef_full), length(coef_full),
                    dimnames = list(names(coef_full),
                                    names(coef_full)))
  V_full[alive, alive] <- V
  V_full
}
