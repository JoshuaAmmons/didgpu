# didgpu 0.1.0

First public release. Five estimator families plus a sensitivity layer,
each with CUDA kernels for the hot paths.

## Callaway-Sant'Anna (2021) — `didgpu_cs()`

- `est_method = c("OR", "IPW", "DR")` — all three inner estimators. DR is
  Sant'Anna-Zhao (2020) doubly-robust.
- `control_group = c("never", "notyet")` — never-treated OR not-yet-treated.
- `covariates =` for OR/IPW/DR adjustment.
- Pre-treatment placebos computed automatically; joint chi-square test
  on the placebo block via `fit$placebo`.
- Four aggregations (`event` / `group` / `calendar` / `overall`) with
  `didgpu_cs_aggregate()` for switching post-fit.
- `bootstrap_kind = c("cluster", "multiplier")` — cluster bootstrap on
  units, or multiplier wild bootstrap on per-unit influence functions
  (much faster for large B).
- CUDA: the OR per-(g, t) inner regression runs on the GPU
  (`src/cuda_cs_inner.cu`, in-thread Cholesky per cell) with per-row
  influence functions; cluster + multiplier bootstrap SEs run on the
  GPU too. IPW / DR inner regressions still use the R path (batched
  logistic-regression kernel is future work).
- Cross-validated against the reference `did` package on simulated
  panels (max abs diff < 0.25 on event-study estimates).

## TestMechs (Kwon & Roth 2026) — `didgpu_test_sharp_null()`

- All three test methods: `"CS"` (Cox-Shi 2023), `"ARP"`
  (Andrews-Roth-Pakes 2023), `"FSST"` (Fang-Santos-Shaikh-Torgovitsky 2023).
- Both binary mediator (K = 2) and multi-level (K >= 2) under
  no-defiers.
- Generic CS engine `.testmechs_cs_test(theta_hat, Sigma, A, A_eq, b_eq)`
  is reusable for any moment-inequality test.
- Nonparametric and Bayesian (Dirichlet) bootstrap of the partial-density
  vector beta.obs.
- CUDA bootstrap kernel `src/cuda_testmechs_bootstrap.cu` (cuRAND
  multinomial; the main acceleration target). Live on Linux/WSL and
  wired through `.testmechs_bootstrap_cuda`; the "nonparametric"
  method runs on the GPU, "bayes" uses the R path. cuRAND vs R's
  MT19937 differ per-replicate, so bootstrap moments match within
  Monte-Carlo error rather than bit-for-bit.

## Leave-one-out robustness — `didgpu_loo()`

- Drops one entity at a time (cohort / unit / cluster / arbitrary
  column level) and re-fits the estimator.
- Works on all three estimator families: `didgpu_result`,
  `didgpu_cs_result`, `didgpu_fect_result`.
- Returns a `didgpu_loo_result` data.frame with `leave_out`,
  `estimate`, `delta`, and `delta_pct` columns, sorted by
  `abs(delta)` descending so the most-influential drop is on top.
- `print()` shows the top-N most influential rows with an
  interpretation hint; `plot()` draws a tornado plot of deltas.
- Default `by = "cohort"` (leave-one-cohort-out, the standard DiD
  diagnostic). Pass `by = "unit"`, `"cluster"`, or any column name
  to drop on a different key.

## HonestDiD (Rambachan & Roth 2023) — `didgpu_honest_did()`

- Sensitivity analysis on event-study DiD estimates.
- `method = c("RM", "M")` — relative-magnitudes OR smoothness bounds.
- Reports the breakdown parameter (smallest Mbar at which the CI
  includes zero) so users can read off how robust their conclusion is.
- Works on both `didgpu_result` (DIDmultiplegtDYN-style) and
  `didgpu_cs_result` (Callaway-Sant'Anna) fits.

## fect family (counterfactual-prediction estimators) — `didgpu_fect()`

## Estimator (bit-identical to `DIDmultiplegtDYN::did_multiplegt_dyn`)

- Binary, multivalued, and continuous treatment.
- `effects`, `placebo`, `switchers = ""/"in"/"out"`, ATE.
- `weight`, `controls`, `trends_nonparam` cohort extension.
- `only_never_switchers`, `same_switchers`, `dont_drop_larger_lower`.
- `normalized = TRUE` (per-unit-of-treatment), `trends_lin = TRUE`
  (linear cohort trends with cumulative-recovery placebos).
- `same_switchers_pl` (placebo-side same-switchers gate; mirrors the
  reference's constraint that it must be paired with `same_switchers`).
- `predict_het` (heterogeneity regression with HC1 robust SEs and joint
  F-test).
- `didgpu_by_path()` for treatment-trajectory subgroup analysis (the
  equivalent of the reference's `by_path` argument).
- Sample-size columns (`N`, `Switchers`, `N.w`, `Switchers.w`) match
  the reference exactly.

## Long-running workflow

- Per-cell checkpointing to disk with atomic writes (`saveRDS`
  tmp + rename) and append-only `manifest.csv`. Resumable on crash.
- `didgpu_resume(checkpoint_dir, df, ...)` — re-invokes with every
  stored arg restored from `meta.json`.
- `didgpu_bootstrap_more(checkpoint_dir, df, extra_reps)` — extend a
  finished run with more bootstrap reps without rework.
- `didgpu_by(df, by_var, ...)` — per-subgroup fits, each with its own
  checkpoint subdirectory.
- `n_workers > 1L` parallelises the bootstrap loop via
  `parallel::makeCluster`; bit-identical to sequential at the same seed.

## Backends

- `"r"` — pure R via `data.table`. 60× faster than the reference at
  200 K rows.
- `"reference"` — delegate to `DIDmultiplegtDYN::did_multiplegt_dyn`,
  used as the parity oracle.
- `"cuda"` — **live on Linux/WSL2** (built + verified end-to-end on an
  NVIDIA RTX 4000 Ada, CUDA 12.6). Live GPU paths: the CS cluster
  bootstrap (**179–228× faster** than R via the influence-function
  shortcut), the CS multiplier bootstrap, the CS OR point estimate
  (bit-exact vs R), and the TestMechs nonparametric bootstrap. The
  fect SVD path is size-gated — it only engages for very large
  balanced panels, since cuSOLVER loses to CPU LAPACK on the small
  matrices typical of fect. Every GPU path falls back transparently
  to R when CUDA is unavailable or would be slower, so `backend =
  "cuda"` is always safe. See `BENCHMARKS.md` and
  `tests/testthat/test-cuda-equivalence-grid.R` (142 lock-step
  assertions). Tests skip GPU paths when `nvcc` / a device is absent.
- `"cpu"` — Rcpp+Eigen, scaffolded only.

## R interface

- S3 methods: `print`, `summary`, `coef`, `confint`, `vcov`, `plot`,
  plus `tidy`, `glance`, `augment` via `broom`.
- Plot is base-R (no `ggplot2` dependency); event-study with stored CIs.
- Diagnostic helpers: `didgpu_summarize_panel`, `didgpu_estimate_runtime`,
  `didgpu_compare` (compare against the reference).

## fect family (counterfactual-prediction estimators)

- `didgpu_fect(method = "fe")` — two-way fixed effects, iterative
  demeaning of the controls-only outcome matrix.
- `didgpu_fect(method = "ife")` — Bai (2009) interactive fixed effects.
  Alternating fe-step + rank-r SVD of the residual matrix until
  convergence.
- `didgpu_fect(method = "mc")` — Athey et al. (2021) matrix completion.
  Iterative soft-thresholded SVD on the controls-only matrix.
- All three reuse `didgpu()`'s checkpoint / resume / parallel bootstrap
  infrastructure. Results are returned as `didgpu_fect_result` (extends
  `didgpu_result`) so all the standard accessors (`coef`, `confint`,
  `vcov`, `plot`, `tidy`, `glance`) work the same way.
- CUDA kernels for fect live in `src/cuda_fect_fe.cu` and
  `src/cuda_fect_svd.cu` (the latter uses cuSOLVER's
  `cusolverDnDgesvdj` for the SVD primitive shared by `ife` and `mc`)
  and are wired through R. **However**, they are size-gated: on the
  small, tall-skinny matrices typical of fect panels the per-iteration
  cuSOLVER SVD is 100–300× slower than R's LAPACK (cuSOLVER handle +
  H2D/D2H overhead dwarfs the tiny SVD). `.fect_cuda_svd_worthwhile()`
  only routes to the GPU for very large balanced panels
  (`n_units ≥ 2000` and `n_units·n_periods ≥ 2e5`); below that
  `backend = "cuda"` transparently uses R's `svd()`. See `BENCHMARKS.md`.

## Testing

- 300+ tests across 24 test files; `R CMD check` passes with only
  pre-existing intentional WARN (CUDA `.cu` files in `src/`) and the
  declared GNU make `SystemRequirements` NOTE.
- Adversarial fuzz harness (`tests/testthat/test-fuzz.R`) covers 21
  scenarios: vanilla / weight / controls / switchers / normalized /
  placebos / trends_lin / only_never + same_switchers / kitchen sink /
  trends_lin sink / multivalued / bootstrap-stability /
  parallel-equals-sequential / checkpoint round-trip / degenerate /
  very-small. Default `DIDGPU_FUZZ_N = 8` for fast CI; bump via env
  var for deep local runs (validated at N = 200, no failures).
