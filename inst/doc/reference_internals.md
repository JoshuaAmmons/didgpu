# DIDmultiplegtDYN 2.2.0 — Internal Structure

This document is the working spec for didgpu's clean-room reimplementation.
It is the result of a deep read of the reference package source (dumped to
`reference_source/` by `reference_source/_dump.R`). All line numbers refer to
those dumped files; they are byte-identical to the installed package source
modulo the `deparse()` reformatting.

Last refreshed against DIDmultiplegtDYN 2.2.0.

---

## 1. Call graph

`did_multiplegt_dyn` (`did_multiplegt_dyn.R:3`) validates args, then per
by-group calls:

```
did_multiplegt_dyn
  -> did_multiplegt_dyn_by_check       (validation only, by_check.R:3)
  -> [optional] did_multiplegt_by_path
  -> did_multiplegt_main
       -> [if controls] feols + MASS::ginv     (main.R:411-476)
       -> did_multiplegt_dyn_core              (1-2x by switcher direction)
            [pure compute: weighted sums, no regressions]
       -> [optional] feols/lm for predict_het  (main.R:1307-1386)
       -> [optional] write.csv                 (main.R:1502)
  -> [optional] did_multiplegt_bootstrap
       -> did_multiplegt_main (per bootstrap iter, full re-run)
  -> _design / _dfs / _normweights / _graph / _save_sample (post-hoc)
```

Key fact: under `trends_lin = FALSE` (the default), the core is called
**twice** per estimation — once per switcher direction `in`/`out` —
**not** once per event-time. One call computes all event-times together.

---

## 2. Pre-core data preparation (in `did_multiplegt_main.R`)

By the time the core is invoked, `df` has been:

- Column-renamed with `_XX` suffix (`main.R:65-67, 117-119, 139-147`).
- Balanced with `plm::pdata.frame` + `make.pbalanced(...,
  balance.type="fill")` (`main.R:262-263, 341`).
- Group/time recoded to consecutive integers (`main.R:145-146`).

Added derived columns include:

| Column | What | Built at |
|---|---|---|
| `d_sq_XX` | baseline (period-1) treatment per group | 174-181 |
| `F_g_XX` | first-switch period per group (`T_max+1` for never-switchers) | 187-200, 224 |
| `S_g_XX` | switcher direction: 1 = up ("in"), 0 = down ("out"), NA = never | 294-302 |
| `L_g_XX` | post-switch horizon = `T_g_XX - F_g_XX + 1` | 322 |
| `N_gt_XX` | effective per-cell weight (0 if outcome/treatment missing) | 270-271 |
| `diff_y_XX` | first difference of outcome | 344-345 |
| `diff_d_XX` | first difference of treatment | 345 |
| `L_u_XX`, `L_a_XX` | global max horizons for "in" and "out" | 496-525 |
| `l_XX`, `l_placebo_XX` | clamped effects/placebo lengths | 526-557 |

For each event-time `k = 1..l_XX`, eight scratch columns are pre-allocated
in main (`U_Gg{k}_plus_XX`, `count{k}_plus_XX`, etc.; lines 586-611) so
the core can write into them in place via `data.table` `:=`.

---

## 3. The inner hot loop

`did_multiplegt_dyn_core` runs `for (i in 1:l_u_a_XX)` at `core.R:60`.
Per event-time `i`, every operation is one of:

- **Long difference**: `diff_y_i_XX = outcome - lag(outcome, i)` by group
  (line 108-109). Fused gather+subtract.
- **Mask construction**: `distance_to_switch_i_XX` is 1 iff
  `t == F_g - 1 + i AND i <= L_g AND S_g == direction AND N_gt_control > 0`
  (lines 197/207). Element-wise on `(G*T)`.
- **Grouped sums**: `N_gt_control_i_XX = sum(never_change * N_gt) by (time, d_sq, trends_nonparam)`
  (line 137), `N1_t_i_XX = sum(distance_to_switch_i_wXX) by time` (line 211),
  `N1_t_i_g_XX` by `(time, d_sq, trends_nonparam)` (line 223).
- **Cohort means**: `mean_cohort_i_ns_t_XX` by `(d_sq, trends_nonparam, time)`
  (line 302), `mean_cohort_i_s_t_XX` by `(d_sq, F_g, d_fg, trends_nonparam)`
  (line 329).
- **The U-statistic** (lines 484-491):

  ```
  U_Gg_i_temp = dummy * (G / N_inc_i) * 1[t in [i+1, T_g]] * N_gt
                * (distance - (N_t_g / N_t_control) * never_change)
                * diff_y_i
  U_Gg_i      = sum_t U_Gg_i_temp  by group
  ```

  This is *not* a regression. It is a weighted DiD construction. Each row
  contributes O(1) work; a per-group reduction sums them.
- **Variance kernel** (lines 498-533): same shape as the U-statistic,
  multiplied by `DOF_gt_i * (diff_y_i - E_hat_gt_i)`.

Aggregation across event-times (lines 893-935):

```
w_i = N_inc_i / sum_j N_inc_j
U_Gg_num     += w_i * U_Gg_i
U_Gg_num_var += w_i * U_Gg_i_var
U_Gg_den     += w_i * delta_D_i      # cumulative "dose"
U_Gg     = U_Gg_num / U_Gg_den
U_Gg_var = U_Gg_num_var / U_Gg_den
```

Placebos (lines 576-891) mirror this structure with
`diff_y_pl_i = lag(outcome, 2i) - lag(outcome, i)`.

**The only true matrix solve in the pipeline lives outside the core**, in
`main.R:411-433`: a small `(K+2) x (K+2)` Gram inverted via `MASS::ginv`
once per baseline level, for the controls residualization. Dropped under
the "no controls" simplification.

---

## 4. Bootstrap

`did_multiplegt_bootstrap.R:26-57`:

```r
for (j in 1:bootstrap) {
  df_boot <- df[ list_to_vec(xtset[sample(1:length(xtset), replace = TRUE)]), ]
  df_boot <- df_boot[order(group, time), ]
  df_est  <- did_multiplegt_main(df_boot, ...)   # SAME function as point
  bresults_effects[j, ] <- df_est$Effects[, 1]
  bresults_ATE[j, 1]    <- df_est$ATE[1]
  bresults_placebo[j, ] <- df_est$Placebos[, 1]
}
```

Cluster bootstrap: `bs_group <- ifelse(!is.null(cluster), cluster, group)`
(line 19). `xtset` is a list of row-indices keyed by `bs_group` level; one
iter samples `length(xtset)` clusters with replacement and concatenates
their rows. Each iter is a full `did_multiplegt_main` re-run — there is
no shortcut. SEs = `sd(bresults_*)`, normal CIs.

**Natural checkpoint unit: one bootstrap iter.** Iters are independent;
per-iter output is small (effects vector + ATE scalar + placebos
vector); the core call cannot be split by event-time without re-running
the controls pre-fit each subset. Per-baseline-level `l` is a separate
parallelisable dimension inside `main` but it's gated on `controls`.

---

## 5. Purity / side effects

- `did_multiplegt_dyn_core`: pure function of inputs. Uses `assign()`
  into its own frame; no global writes, no I/O. Mutates input `df` in
  place via `data.table` `:=` — callers reassign `df <- data$df`.
- `did_multiplegt_bootstrap`: no global state, no I/O. Calls
  `progressBar()` (just `cat`). Non-deterministic without external
  `set.seed()`.
- `did_multiplegt_main`: pure except optional `write.csv` at line 1502
  when `save_results` is non-NULL.

---

## 6. Switchers pre-processing in code

The conceptual "switchers/stayers/yet-to-switch" machinery reduces to:

1. `d_sq_XX` — baseline treatment (cohort key).
2. `diff_from_sq_XX = treatment - d_sq`.
3. Monotonicity tags: `ever_strict_increase_XX`, `ever_strict_decrease_XX`.
   Non-monotone groups dropped (unless `dont_drop_larger_lower`).
4. `F_g_XX` — first switch date; `T_max+1` for never-switchers.
5. Drop cohorts (by `d_sq, trends_nonparam`) with no `F_g` variation.
6. Missing-treatment linear fill, panel truncation.
7. `avg_post_switch_treat_XX`, then `S_g_XX = 1[avg_post > d_sq]`.
8. `L_g_XX, L_u_XX, L_a_XX, L_placebo_u_XX, L_placebo_a_XX`.

That's it. A per-group integer trio `(d_sq, F_g, S_g)` plus per-group
horizon `L_g`. The two-direction core dispatch slices on `S_g`.

---

## 7. Memory hotspots

Per core call, the data.table `df` grows by:

- ~25 columns per event-time `i` (the `U_Gg_i_*`, mask, count, DOF cols).
- ~25 per placebo `i` (mirrored).
- ~15 * `l_XX * |d_sq levels| * K` for controls (zero in our case).

Peak is at end of the main `i` loop. For G=1000, T=20, l=10, placebo=5:
~20k rows * ~600 columns = ~100 MB per core call. Per bootstrap iter,
double that (resampled `df_boot` + main's enlarged copy).

The largest dense intermediates outside the data.table are `data_XX`
(per-baseline Gram subset, `main.R:401`) and the `feols` design matrix
at `main.R:476` — both controls-only.

---

## 8. Binary non-absorbing reduction — branches to drop

For our target (binary on/off treatment, no controls, no continuous,
no normalized, no by_path, no predict_het, no trends_lin), the dead
branches we ignore in the reimplementation:

**`did_multiplegt_dyn.R`**: by_path (117-129, 138-142, 199-205,
214-221); design/dfs/normweights/graph post-processors (168-187); XLSX
writers (207-212).

**`did_multiplegt_main.R`**: predict_het (90-112, 247-251, 1252-1403);
continuous (61-62, 202-208, 297-318, 541-544, 855-858); controls
machinery including ginv solve and feols (77-81, 348-491); trends_lin
(252-260 and inline conditionals); normalized / normalized_weights
(621-625, 636-641, 754-757, 809-811, 894-896, 1054-1056, 1089-1091);
effects_equal test (1404-1457).

**`did_multiplegt_dyn_core.R`**: less_conservative_se (111-128 inline +
340-390); controls (225-273, 504-533, 659-688, 818-848); normalized
(536-551, 850-865, 942-950); trends_lin (553-575, 868-891);
continuous in `delta_temp` (912-922); same_switchers enforcement of
`still_switcher_i_XX` (139-199 — take the simpler `else` at 201-209).

**Keep both** `switchers_core="in"` (line 685) and `switchers_core="out"`
(line 773) calls — binary non-absorbing can have downward switchers.

What remains is precisely the U-statistic accumulation: long
differences + two grouped reductions per event-time + per-group sum,
all element-wise on length-`G*T` vectors. **This is the GPU target.**

---

## Items needing further investigation

- Exact behaviour of `S_g_XX = NA` in some count expressions
  (`main.R:937-940`).
- Whether `feols`'s out-of-sample `predict()` is bit-reproducible across
  versions (matters only if/when we add controls support).
- `main.R:920` references `first_obs_by_gp` (no `_XX` suffix) which
  appears to be a typo in the reference — the unsuffixed column is
  never defined, so the line is a no-op. Verify this conclusion in a
  separate pass before depending on it.
- Placebo `same_switchers_pl` cross-terms are large enough to deserve a
  dedicated read pass before reimplementing.
