# ============================================================================
# didgpu_bacon(): Goodman-Bacon (2021) decomposition of the static TWFE DiD.
#
# The two-way fixed-effects estimator beta^DD on a binary, staggered,
# absorbing treatment equals a weighted average of all 2x2 DiD comparisons
# between "timing groups". This returns that decomposition: every 2x2's
# estimate and weight, grouped into
#   - "Treated vs Untreated"        : a timing group vs the never-treated   (clean)
#   - "Earlier vs Later Treated"    : earlier group treated, later group is
#                                     control while still untreated          (clean)
#   - "Later vs Earlier Treated"    : later group treated, EARLIER group is
#                                     an already-treated control            (forbidden)
# The total weight on the "forbidden" already-treated comparisons is the
# Goodman-Bacon diagnostic for how heavily TWFE leans on bad comparisons —
# the bias channel under dynamic / heterogeneous treatment effects.
#
# Cheap: it's group-level cell means + closed-form weights, no bootstrap,
# no iteration beyond the one demeaning used to report the TWFE coefficient.
#
# Reference: Goodman-Bacon (2021), "Difference-in-differences with variation
# in treatment timing", J. Econometrics 225(2). Canonical balanced-panel,
# binary, staggered-adoption case (the same scope as the bacondecomp pkg).
# ============================================================================


#' Goodman-Bacon decomposition of the TWFE DiD estimator
#'
#' Decomposes the static two-way fixed-effects DiD coefficient into the
#' weighted average of all 2x2 timing-group comparisons (Goodman-Bacon
#' 2021). The headline diagnostic is the total weight placed on "forbidden"
#' comparisons that use already-treated units as controls — the source of
#' TWFE bias under heterogeneous, dynamic treatment effects. Pair it with
#' [didgpu_twfe()] (the estimate) and [didgpu()] / [didgpu_cs()] (robust
#' alternatives).
#'
#' Scope: a **balanced** panel with a **binary, absorbing** (staggered-
#' adoption) treatment — the canonical Goodman-Bacon case. Always-treated
#' units (treated in the first period, no pre-period) are dropped with a
#' note, as in the reference implementation. For unbalanced panels,
#' treatment that switches off, or continuous treatment, use [didgpu()].
#'
#' @param df A data.frame / data.table panel.
#' @param outcome,group,time,treatment Character column names: outcome,
#'   unit id, time id, and the binary 0/1 absorbing treatment indicator.
#' @return A `didgpu_bacon` object: a list with `beta_twfe` (the static
#'   TWFE DiD coefficient on the decomposition sample), `comparisons` (a
#'   data.frame: `type`, `treated`, `control` timing values, `estimate`,
#'   `weight`), `summary` (per-type total weight and weighted-average 2x2
#'   estimate), `forbidden_weight` (total weight on already-treated-control
#'   comparisons), and `n_always_treated` (units dropped).
#' @seealso [didgpu_twfe()], [didgpu()], [didgpu_cs()].
#' @examples
#' \donttest{
#' p <- didgpu_simulate_panel(n_units = 80L, n_periods = 12L,
#'                            tau_profile = c(0.5, 1.0), seed = 7L)
#' p$D <- as.integer(p$D >= 0.5)            # binary staggered treatment
#' bd <- didgpu_bacon(p, "Y", "unit", "period", "D")
#' bd
#' }
#' @export
didgpu_bacon <- function(df, outcome, group, time, treatment) {
  stopifnot(is.data.frame(df) || data.table::is.data.table(df))
  for (nm in c("outcome", "group", "time", "treatment")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L || !nzchar(v))
      stop("`", nm, "` must be a single non-empty column name.")
    if (!v %in% names(df)) stop("column not in df: ", v)
  }
  d <- data.table::as.data.table(df)
  d <- d[!is.na(get(outcome)) & !is.na(get(group)) &
         !is.na(get(time)) & !is.na(get(treatment))]
  Y <- as.numeric(d[[outcome]]); G <- d[[group]]; Tt <- d[[time]]
  D <- as.numeric(d[[treatment]])

  if (!all(D %in% c(0, 1)))
    stop("didgpu_bacon requires a binary 0/1 treatment; got values outside {0,1}. ",
         "For continuous treatment use didgpu().")

  units <- sort(unique(G)); periods <- sort(unique(Tt))
  nU_all <- length(units); Tn <- length(periods)
  if (nrow(d) != nU_all * Tn)
    stop("didgpu_bacon requires a BALANCED panel (every unit observed in ",
         "every period). Use didgpu() for unbalanced data.")

  ri <- match(G, units); ci <- match(Tt, periods)
  Ymat <- matrix(NA_real_, nU_all, Tn); Dmat <- matrix(NA_real_, nU_all, Tn)
  Ymat[cbind(ri, ci)] <- Y; Dmat[cbind(ri, ci)] <- D

  if (any(apply(Dmat, 1L, function(r) any(diff(r) < 0))))
    stop("didgpu_bacon requires an ABSORBING (staggered) treatment: once a ",
         "unit is treated it stays treated. For treatment that turns off, ",
         "use didgpu() (it handles in/out switchers).")

  # First treated column index per unit; Inf = never treated.
  first_tr <- apply(Dmat, 1L, function(r) { w <- which(r == 1)[1L]
                                            if (is.na(w)) Inf else w })
  always <- which(first_tr == 1L)            # treated from period 1: no pre-period
  n_always <- length(always)
  if (n_always == nU_all)
    stop("every unit is always-treated; nothing to decompose.")
  if (n_always > 0L) {
    keep <- setdiff(seq_len(nU_all), always)
    Ymat <- Ymat[keep, , drop = FALSE]
    Dmat <- Dmat[keep, , drop = FALSE]
    first_tr <- first_tr[keep]
  }
  N <- nrow(Ymat)

  # --- Static TWFE coefficient on the decomposition sample (FWL) ---------
  # beta = sum(D~ * Y~) / sum(D~ * D~), ~ = two-way demeaned. Reuses the
  # same iterative demeaning as didgpu_twfe() so it matches lm() exactly.
  Yv <- as.numeric(t(Ymat)); Dv <- as.numeric(t(Dmat))   # row-major (unit, time)
  g_idx <- rep(seq_len(N), each = Tn)
  t_idx <- rep(seq_len(Tn), times = N)
  Md <- .twfe_demean(cbind(Yv, Dv), g_idx, t_idx, tol = 1e-12, max_iter = 5000L)
  Yd <- Md[, 1L]; Dd <- Md[, 2L]
  vDD <- sum(Dd * Dd)
  beta_twfe <- if (vDD > 0) sum(Dd * Yd) / vDD else NA_real_

  # --- Timing groups -----------------------------------------------------
  fin <- is.finite(first_tr)
  g_times <- sort(unique(first_tr[fin]))     # column indices of adoption
  has_U <- any(!fin)
  units_of <- function(g) which(first_tr == g)
  n_of  <- function(g) length(units_of(g)) / N
  Dbar  <- function(g) (Tn - g + 1) / Tn     # share of periods treated
  cellY <- function(uset, cols) {
    if (length(uset) == 0L || length(cols) == 0L) return(NA_real_)
    mean(Ymat[uset, cols, drop = FALSE])
  }
  # Map column index -> the actual time value, for human-readable output.
  tval <- function(ci) periods[ci]

  comps <- list()
  add <- function(type, treated, control, estimate, w_num) {
    comps[[length(comps) + 1L]] <<- data.frame(
      type = type, treated = treated, control = control,
      estimate = estimate, w_num = w_num, stringsAsFactors = FALSE)
  }

  U_units <- if (has_U) which(!fin) else integer(0)
  n_U <- length(U_units) / N

  for (gi in seq_along(g_times)) {
    k <- g_times[gi]
    uk <- units_of(k); nk <- length(uk) / N; Dk <- Dbar(k)
    pre_k  <- seq_len(k - 1L)                 # t < adoption_k
    post_k <- k:Tn                            # t >= adoption_k

    # Treated vs Never-treated.
    if (has_U) {
      est <- (cellY(uk, post_k) - cellY(uk, pre_k)) -
             (cellY(U_units, post_k) - cellY(U_units, pre_k))
      add("Treated vs Untreated", tval(k), NA, est, nk * n_U * Dk * (1 - Dk))
    }

    # Pairs with a later-adopting group l (k earlier, l later).
    for (gj in seq_along(g_times)) {
      l <- g_times[gj]
      if (l <= k) next
      ul <- units_of(l); nl <- length(ul) / N; Dl <- Dbar(l)
      MID  <- k:(l - 1L)                      # k treated, l not yet
      POST <- l:Tn                            # both treated

      # Earlier (k) treated vs later (l) still-untreated control — CLEAN.
      est_good <- (cellY(uk, MID) - cellY(uk, pre_k)) -
                  (cellY(ul, MID) - cellY(ul, pre_k))
      add("Earlier vs Later Treated", tval(k), tval(l), est_good,
          nk * nl * (Dk - Dl) * (1 - Dk))

      # Later (l) treated vs earlier (k) ALREADY-treated control — FORBIDDEN.
      est_bad <- (cellY(ul, POST) - cellY(ul, MID)) -
                 (cellY(uk, POST) - cellY(uk, MID))
      add("Later vs Earlier Treated", tval(l), tval(k), est_bad,
          nk * nl * Dl * (Dk - Dl))
    }
  }

  comp <- do.call(rbind, comps)
  total_w <- sum(comp$w_num)
  comp$weight <- comp$w_num / total_w
  comp$w_num <- NULL
  rownames(comp) <- NULL

  # Per-type summary: total weight + weighted-average 2x2 estimate.
  types <- unique(comp$type)
  summ <- do.call(rbind, lapply(types, function(ty) {
    sub <- comp[comp$type == ty, , drop = FALSE]
    w <- sum(sub$weight)
    avg <- if (w > 0) sum(sub$weight * sub$estimate) / w else NA_real_
    data.frame(type = ty, weight = w, avg_estimate = avg,
               n_comparisons = nrow(sub), stringsAsFactors = FALSE)
  }))
  rownames(summ) <- NULL

  beta_check <- sum(comp$weight * comp$estimate)
  forbidden_w <- sum(comp$weight[comp$type == "Later vs Earlier Treated"])

  out <- list(
    beta_twfe        = beta_twfe,
    beta_check       = beta_check,     # == beta_twfe (Goodman-Bacon identity)
    comparisons      = comp,
    summary          = summ,
    forbidden_weight = forbidden_w,
    n_always_treated = n_always,
    n_units          = N,
    n_periods        = Tn,
    has_never_treated = has_U
  )
  class(out) <- "didgpu_bacon"
  out
}


#' Print method for didgpu_bacon
#' @param x A `didgpu_bacon` object from [didgpu_bacon()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.didgpu_bacon <- function(x, ...) {
  cat("Goodman-Bacon decomposition of the TWFE DiD estimator\n")
  cat(sprintf("  %d units, %d periods%s%s\n", x$n_units, x$n_periods,
              if (x$has_never_treated) ", with never-treated controls" else "",
              if (x$n_always_treated > 0L)
                sprintf("; dropped %d always-treated unit(s)", x$n_always_treated)
              else ""))
  cat(sprintf("  TWFE DiD estimate: %.6f\n", x$beta_twfe))
  cat("\nDecomposition by comparison type (weight * estimate):\n")
  s <- x$summary
  for (i in seq_len(nrow(s))) {
    cat(sprintf("  %-26s  weight %6.3f   avg 2x2 %9.5f\n",
                s$type[i], s$weight[i], s$avg_estimate[i]))
  }
  cat(sprintf("\n  Weight on FORBIDDEN (already-treated control): %.3f\n",
              x$forbidden_weight))
  cat("  Higher forbidden weight => TWFE leans more on bad comparisons and\n")
  cat("  is more exposed to bias under heterogeneous/dynamic effects.\n")
  cat("  Robust alternatives: didgpu(), didgpu_cs().\n")
  # Sanity: the weighted 2x2 sum reproduces the TWFE coefficient.
  if (is.finite(x$beta_twfe) &&
      abs(x$beta_check - x$beta_twfe) > 1e-6 * (1 + abs(x$beta_twfe)))
    cat(sprintf("\n  [warning] decomposition sum %.6f != TWFE %.6f\n",
                x$beta_check, x$beta_twfe))
  invisible(x)
}
