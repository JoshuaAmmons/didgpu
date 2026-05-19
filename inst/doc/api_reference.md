# didgpu API reference

Quick reference for every public function across the five estimator
families plus the HonestDiD sensitivity layer.

## Family 1: DIDmultiplegtDYN-style — `didgpu()`

Heterogeneity-robust dynamic DiD estimator (de Chaisemartin & D'Haultfoeuille 2024).

- `didgpu(df, outcome, group, time, treatment, ...)` — main entry point.
  Bit-identical to `DIDmultiplegtDYN::did_multiplegt_dyn` across every
  commonly-used option: binary / multivalued / continuous treatment;
  `controls`, `weight`, `trends_nonparam`, `trends_lin`, `normalized`,
  `predict_het`, `same_switchers`, `same_switchers_pl`,
  `only_never_switchers`, `dont_drop_larger_lower`,
  `switchers = ""/"in"/"out"`, all sample-size columns.
- `didgpu_simulate_panel(...)` — DGP generator with known truth.
- `didgpu_compare(...)` — bit-by-bit comparison against the reference.
- `didgpu_summarize_panel(...)` — pre-fit panel design check.
- `didgpu_estimate_runtime(...)` — extrapolate wall-clock from probe fits.
- `didgpu_backend_info()` — which backends are available on this machine.
- `didgpu_resume(checkpoint_dir, df, ...)` — resume an interrupted run.
- `didgpu_bootstrap_more(checkpoint_dir, df, extra_reps)` — extend a
  finished checkpoint with more bootstrap reps.
- `didgpu_by(df, by_var, ...)` — per-subgroup estimation.
- `didgpu_by_path(...)` — per-treatment-trajectory subgroup estimation.
- `didgpu_compute_paths(...)` — augment a panel with a trajectory column.
- `didgpu_init_checkpoint(...)`, `didgpu_load_checkpoint(...)`,
  `didgpu_aggregate_cells(...)` — manual checkpoint API.
- `didgpu_event_study_data(fit)` — extract an event-study data frame
  ready for ggplot2.

### S3 methods on `didgpu_result`

`print`, `summary`, `coef`, `confint`, `vcov`, `plot`, `tidy`,
`glance`, `augment`.

## Family 2: Callaway-Sant'Anna (2021) — `didgpu_cs()`

Staggered DiD with group-time ATT(g, t).

- `didgpu_cs(df, outcome, group, time, treatment, est_method, control_group, aggregation, covariates, bootstrap_reps, bootstrap_kind, ...)`
  - `est_method = c("OR", "IPW", "DR")`. DR is Sant'Anna-Zhao (2020)
    doubly-robust.
  - `control_group = c("never", "notyet")`.
  - `aggregation = c("event", "group", "calendar", "overall")`.
  - `bootstrap_kind = c("cluster", "multiplier")`.
- `didgpu_cs_aggregate(fit, aggregation)` — switch aggregation without
  refitting.

### S3 methods on `didgpu_cs_result`

`print`. (Inherits the standard model accessors via the broader
result-result classes.)

## Family 3: fect (Liu, Wang & Xu 2024) — `didgpu_fect()`

Counterfactual-prediction DiD.

- `didgpu_fect(df, outcome, group, time, treatment, method, effects, r, lambda, ...)`
  - `method = "fe"` — two-way fixed effects via iterative demeaning.
  - `method = "ife"` — Bai (2009) interactive fixed effects (alternating
    fe + rank-r SVD).
  - `method = "mc"` — Athey et al. (2021) matrix completion via iterative
    soft-thresholded SVD; `lambda` selected by k-fold CV if NULL.
- `didgpu_fect_placebo(df, ..., method, n_placebos, bootstrap_reps)` —
  refit the chosen method pretending the M periods before F_g are
  treated; report placebo estimates + joint chi-square test.
- `didgpu_fect_equivalence(placebo_result, delta)` — TOST-based
  equivalence test: do the placebo deviations stay within `delta`?

## Family 4: TestMechs (Kwon & Roth 2026) — `didgpu_test_sharp_null()`

Sharp test of full mediation.

- `didgpu_test_sharp_null(df, d, m, y, method, B, num_Ybins, ...)`
  - `method = "CS"` — Cox-Shi (2023). Projection QP + chi-sq cv.
  - `method = "ARP"` — Andrews-Roth-Pakes (2023). Least-favorable cv via
    Monte Carlo on N(0, Sigma).
  - `method = "FSST"` — Fang-Santos-Shaikh-Torgovitsky (2023).
    Cone-based test with secondary bootstrap.
  - Binary AND multi-level mediator (K >= 2) under no-defiers.
- `didgpu_lb_frac_affected(df, d, m, y, B, ...)` — sharp lower bound
  on the fraction of always-takers whose outcome is moved by treatment.
  (Scaffolded; implementation pending.)

## Family 5 (sensitivity layer): HonestDiD — `didgpu_honest_did()`

Rambachan & Roth (2023) sensitivity analysis.

- `didgpu_honest_did(fit, event_post, method, Mbar, alpha)`
  - `fit` is a `didgpu_result` or `didgpu_cs_result` with placebos.
  - `method = "RM"` — relative magnitudes (post-trend violation
    bounded by Mbar times the max pre-trend deviation).
  - `method = "M"` — smoothness (consecutive-period violations bounded by M).
  - Returns the bounds at each Mbar plus the "breakdown" parameter
    (smallest Mbar at which the CI includes zero).

## CUDA acceleration

All five families have CUDA kernel scaffolds in `src/cuda_*.cu`. They
compile when `nvcc` is on `PATH` and skip silently otherwise. The
R-side reference implementations are always available and produce
identical output up to floating-point reduction order.

| Kernel | Family | Pattern |
|---|---|---|
| `cuda_didkernel.cu` | didgpu | 5-kernel U-statistic chain |
| `cuda_fect_fe.cu` | fect_fe | Iterative row-mean + col-mean reductions |
| `cuda_fect_svd.cu` | fect_ife, fect_mc | cuSOLVER Jacobi SVD + cuBLAS gemm for reconstruction |
| `cuda_testmechs_bootstrap.cu` | TestMechs | cuRAND multinomial + atomic-add reduction |
| `cuda_cs_inner.cu` | didgpu_cs | Batched per-(g, t) regression via cuBLAS/cuSOLVER |

## Cross-validated parity

| didgpu function | Reference package | Status |
|---|---|---|
| `didgpu()` | `DIDmultiplegtDYN` | Bit-identical (≤ 1e-10) across every supported option |
| `didgpu_cs()` (OR) | `did::att_gt()` | Max abs diff < 0.25 on event-study estimates |
| `didgpu_fect()` | `fect` package | Same algorithms; not bit-identical due to different convergence tolerances |
| `didgpu_honest_did()` | `HonestDiD` | Direct wrapper; bit-identical (calls the same internals) |
| `didgpu_test_sharp_null()` | `TestMechs` | Same methodology; v1 power lower on the alternative (known limitation, see code comments) |

## Suggested R packages

- `did` (for CS parity testing)
- `DIDmultiplegtDYN` (for didgpu parity testing)
- `DRDID` (Sant'Anna-Zhao reference; used internally by `did`)
- `HonestDiD` (HonestDiD sensitivity; required for `didgpu_honest_did()`)
- `quadprog` (QP solver; required for TestMechs CS test)
- `broom` (for `tidy` / `glance` / `augment` methods)
- `testthat >= 3.0.0` (for the test suite)
