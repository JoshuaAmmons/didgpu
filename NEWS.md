# didgpu 0.1.2

## Bug fixes

- **`didgpu_fect()` kept units whose counterfactual is not identified,
  which made `method = "ife"` badly wrong -- sometimes the wrong sign.**
  A unit's counterfactual comes only from its UNTREATED observations:
  the unit fixed effect needs at least one, an r-factor loading needs
  several. fect drops units below a minimum (`fect.default`): `min.T0`
  is `1` for `method = "fe"` and `5` for `"ife"`, `"mc"`, `"both"`,
  `"gsynth"` and `"cfe"`. didgpu kept every unit and extrapolated.

  With a weak factor structure that is a nuisance; with a strong one it
  is fatal. On a 100-unit panel with two latent factors and a TRUE ATT
  of +1.0:

      factors    truth   fect       didgpu before   didgpu after
      L ~ 0.5     1.0    +1.00142   +0.97446        +1.00206
      L ~ 1.5     1.0    +0.99880   +0.72146        +0.98872
      L ~ 3.0     1.0    +1.00064   -0.17445        +0.94980

  The clue was that `fect`'s returned `lambda` was 84 x 2 on a 100-unit
  panel: it had silently dropped 16 units, exactly those with fewer
  than 5 untreated periods.

  `didgpu_fect()` gains a documented `min_T0` argument. `NULL` (the
  default) follows fect's rule; `$n_units_dropped` and `$min_T0` are
  returned and dropped units are warned about. This subsumes the
  earlier always-treated fix, since those units have zero untreated
  periods.

  Note a consequence: because the defaults differ by method,
  `method = "ife", r = 0` does NOT generally equal `method = "fe"` --
  they fit different samples. `fect` behaves the same way, and didgpu
  now reproduces both of its numbers on a panel where 22 of 60 units
  have short histories (`fe` +0.302820, `ife(r = 0)` +0.236833). Pass
  `min_T0 = 1` to make them agree.

- **`didgpu_fect(method = "fe")` stopped before converging.** The fit
  had a second stopping rule, `abs(loss - prev_loss) < tol`, where
  `loss` is a SUM of squared residuals. Its absolute change falls below
  a tolerance meant for parameter units long before the parameters
  settle, so it fired first: on a 60x10 panel the fit exited after 6
  iterations with `delta = 2.0e-04` against the requested `1e-05` and
  reported itself finished. Convergence is now judged on the
  parameters alone. `fe` converges in 9 iterations and its agreement
  with `fect::fect` improved from 2.81e-05 to 7.80e-07.

## New features

- **`didgpu_fect()` reports convergence diagnostics and warns when a
  fit does not converge.** The solver already produced `iter`, `delta`
  and `lambda` per cell but nothing surfaced them, so a caller could
  not tell a converged fit from one that had exhausted `max_iter` --
  both returned a number and looked identical. Results now carry
  `$diagnostics` (`iter`, `delta`, `converged`, `tol`, `max_iter`,
  `lambda`, `n_nonzero_singular`) and a non-converged fit warns.

  This immediately exposed two problems that had been invisible: `fe`
  was exiting after 6 iterations (fixed above), and `method = "ife"`
  exhausts `max_iter` without converging at default settings. The
  latter is a known outstanding defect -- `ife` is materially wrong as
  the factor structure strengthens -- and is not yet fixed.

- **The CUDA backend silently ignored most estimation options.** Its
  compatibility guard tested only `controls`, `weight` and
  `trends_nonparam`, so every other option passed through to a kernel
  that does not implement it and a wrong number came back with no
  warning and no fallback.

  `normalized` was the damaging case. With a multivalued treatment,
  `backend = "cuda"` returned the UNnormalised effects, which do not
  vary with dose -- so a dose-response analysis looked exactly as
  though the treatment had been binarised. Against
  `DIDmultiplegtDYN` with `normalized = TRUE` on a three-level dose:

      backend   |diff| vs reference
      r         5.551e-17
      cpu       5.551e-17
      cuda      8.465e-01

  `backend = "auto"` resolves to `"cuda"` whenever a GPU is present, so
  this was the default path on a CUDA machine. `dont_drop_larger_lower`
  was separately dropped, because the CUDA path called `.prep_panel()`
  without forwarding it.

  The guard is now identical to the CPU backend's and names the
  offending option when it falls back. All backends now agree with the
  reference to 5.551e-17 on a normalised multivalued treatment.

  To be explicit, since this was reported as binarisation: `didgpu()`
  does NOT binarise the treatment. Under the default
  `normalized = FALSE` the dCDH dynamic effect is the average outcome
  change among switchers, keyed on switching TIMES rather than dose
  magnitude, so rescaling a treatment that keeps the same switching
  pattern leaves the estimate unchanged -- in `DIDmultiplegtDYN` too.
  Dose enters under `normalized = TRUE`, where it moved the estimate by
  0.56 and matched the reference exactly.

- **`didgpu_fect(method = "mc")` ignored fixed effects and was not
  invariant to the level of the outcome.** Athey et al. (2021) estimate
  `Y = L + unit FE + time FE`, penalising the nuclear norm of `L`
  alone. The implementation soft-thresholded the RAW outcome matrix, so
  the penalty shrank the level itself and the residual `Y - Y_hat`
  absorbed it. `lambda` compounded this by being scaled to the singular
  values of the raw matrix, so the penalty also grew with the level.

  On a known-zero DGP the reported ATT moved with a pure location
  shift -- `Y + 0` gave `+0.269`, `Y + 100` gave `+2.551` -- while
  `fe`, `ife` and `fect::fect` all returned `-0.024` at every level. On
  a positive, trending outcome it manufactured large, monotonically
  rising, significant effects where `fe` and `ife` both found a null.

  The fit now removes two-way fixed effects, soft-thresholds the SVD of
  the RESIDUAL, and adds the fixed effects back; `lambda` is scaled to
  that residual. The estimator is now exactly level-invariant, recovers
  the known-zero null (`-0.031`), and tracks `fect::fect` to 6.9e-03
  (was 2.9e-01). It is not yet bit-for-bit: `lambda` selection still
  differs from `fect`'s own cross-validation.

- **`didgpu_loo()` reported a pre-treatment placebo instead of the
  ATT.** The headline for a `didgpu_cs_result` was
  `fit$aggregation$estimate[1]`, documented as working "for all four
  aggregations". It does not: with the default `aggregation = "event"`
  row 1 is the MOST NEGATIVE event time, i.e. the longest
  pre-treatment horizon.

  It also concealed itself, because dropping a single entity seldom
  changes which cells populate the earliest lead -- so nearly every
  entity returned an IDENTICAL estimate, which reads as "no entity is
  influential" rather than as a bug. On a 56-year panel with a 1985
  cohort it returned the `e = -42` cell (`-0.006317`) instead of the
  ATT (`-0.032831`), for both `by = "cohort"` and `by = "unit"`.

  `didgpu_loo()` now always reports the OVERALL ATT -- the
  `n_treated`-weighted mean over post-treatment cells -- regardless of
  which aggregation the fit requested, so the answer no longer depends
  on an unrelated display choice. Both the fast re-aggregation path and
  the refit path were affected and both are fixed.

- **`didgpu_cs()` influence functions were wrong; multiplier-bootstrap
  standard errors were two to eight times too narrow.** The per-cell
  influence function was taken to be the treated units' demeaned
  residual. It is not. It needs (a) the treated arm normalised by
  `E[D]`, (b) the comparison arm normalised by `E[p(X)(1-D)/(1-p(X))]`
  -- not `E[D]` -- and (c) estimation-effect terms for the nuisance
  parameters, which load onto CONTROL units. `OR` set every control
  unit's influence to zero outright; `IPW`/`DR` omitted the
  normalisers. Measured against `did::att_gt()` on a 200-unit panel,
  `bootstrap_kind = "multiplier"` returned SEs at ~0.13x (`OR`) and
  ~0.50x (`IPW`/`DR`) of the correct width. Point estimates were
  correct throughout, so nothing looked wrong. The default
  `bootstrap_kind = "cluster"` never touches these and was correct.

  The three per-cell estimators now mirror DRDID -- the package `did`
  itself calls -- function for function (`reg_did_panel`,
  `std_ipw_did_panel`, `drdid_panel`). Verified against `did::att_gt()`
  to machine precision for all three methods, with and without
  covariates: max |diff| ~4e-16 on ATT(g, t) and ~1e-16 on its SE.

  Anyone who used `bootstrap_kind = "multiplier"` should re-run:
  confidence intervals were far too narrow and p-values far too small.

- **The batched CUDA inner kernel for `didgpu_cs()` is disabled pending
  a matching rewrite.** Its ATT agrees with the CPU path to ~4e-16, but
  it computes the OLD influence functions, and those now feed both the
  multiplier bootstrap and the aggregation SEs. After the CPU rewrite,
  per-cell `max |dSE|` against `did::att_gt()` was 1.11e-16 for
  backends `"r"` and `"cpu"` but 1.96e-01 for `"cuda"`. CS therefore
  routes through the validated CPU path, so all backends -- `"r"`,
  `"cpu"` and `"cuda"` -- now return identical numbers; verified at
  cell `max |dSE|` = 1.11e-16 and aggregate `max |dSE|` = 0 for every
  backend. Re-enable for kernel development with
  `options(didgpu.cs_cuda_inner = TRUE)`.

## New features

- **`didgpu_cs()` aggregations now carry standard errors and confidence
  intervals.** `.cs_aggregate()` previously propagated point estimates
  only, so `$aggregation` had no `se` column (while its empty-result
  stub declared one) and `didgpu_tidy()` could only report `NA`.
  Aggregate SEs are now derived from the influence functions the same
  way `did::aggte()` derives them, including the correction for having
  ESTIMATED the aggregation weights (`did:::wif`) -- without which SEs
  degrade badly at long event times, where few cohorts contribute
  (0.98x of correct at event 0, falling to 0.33x at event 9).
  Verified identical to `did::aggte(type = "dynamic")` for all three
  estimators, with and without covariates: max |diff| ~5e-16 on the
  aggregate estimate and ~3e-17 on its SE.

  Note that pre-treatment ATT(g, t) still differ from `did` by
  construction: didgpu uses a universal base period, `did` defaults to
  a varying one. Parity is asserted on post-treatment cells and
  non-negative event times.

- **`didgpu_fect(method = "ife", r = 0)` now runs.** An
  interactive-fixed-effects model with zero factors is the plain
  two-way FE model, and sweeping `r = 0..k` is the standard way to ask
  whether a result depends on the factor structure. It previously
  errored with `requires numeric/complex matrix/vector arguments`:
  `svd(M, nu = 0, nv = 0)` omits the `u` component entirely (it is
  present only when `nu > 0`), so the truncated-SVD helper multiplied
  `NULL`. `r = 0` now returns an empty factor term and the fit reduces
  to two-way FE -- verified identical to `method = "fe"`.

- **`didgpu_tidy()` accepts a `didgpu_cs_result`.** The CS result class
  does not carry the `didgpu_result` parent, so tidy failed its
  `stopifnot()` with the opaque `inherits(x, "didgpu_result") is not
  TRUE`, even though `didgpu_loo()` and `didgpu_honest_did()` both
  document and accept that class. Tidy now emits one row per
  `aggregation` entry (`kind = "aggregate"`) followed by one row per
  ATT(g, t) cell (`kind = "att_gt"`). `.cs_aggregate()` propagates
  point estimates only, so `std.error` and the CI columns are `NA` on
  the aggregate rows rather than fabricated; the cell rows carry real
  SEs and CIs when `bootstrap_reps > 0`. Passing anything else now
  names both accepted classes instead of printing a bare `inherits()`
  assertion.

- **`didgpu_cs()` cluster bootstrap no longer crashes on panels whose
  columns collide with internal variable names.** `.cs_bootstrap_se()`
  held the panel as a `data.table` and subset it with
  `d[d[[args$group]] == u, ]`. `[.data.table` evaluates its `i`
  expression with the table's COLUMNS in scope, so a panel carrying a
  column literally named `d` shadowed the local `d` with the treatment
  VECTOR: `d[[args$group]]` became `treatment[["g"]]` and failed with
  `subscript out of bounds`. Since `d` is an ordinary name for a
  treatment indicator, real panels hit this routinely, and because it
  fired only when `bootstrap_reps > 0` it looked data-dependent rather
  than name-dependent -- the same panel worked at `reps = 0` and
  crashed at `reps > 0`. didgpu's own simulated panels never tripped it
  because their treatment column is `D`.

  The bootstrap now keeps the panel as a plain data.frame and resolves
  every column lookup before indexing, so no column name can shadow an
  internal. Point estimates and SEs are bit-identical on panels that
  previously worked. As a side effect the per-unit row lookup is
  precomputed once instead of rescanning the whole panel for every
  (replicate x pick), removing an O(B * n_units * nrow) cost.

- **`didgpu_fect()` no longer reports always-treated units' levels as
  treatment effects.** A unit treated in every observed period has no
  control cell, so its unit fixed effect (`fe`) / factor loading
  (`ife`, `mc`) is unidentified. The fitter set an unidentified unit
  effect to `0`, which made the imputed counterfactual the time effect
  alone -- so the unit's entire LEVEL landed in the residual and was
  reported as treatment effect. Always-treated units are selected on
  level (they are exactly the units already treated before the sample
  window opened), so the bias did not average out. On a known-zero DGP
  with 10 such units the reported ATE was **+1.63 against a true effect
  of 0**, and it survived every factor count and `method = "fe"`.
  Units with no untreated period are now dropped before estimation
  with a warning, and the count is returned as
  `$n_always_treated_dropped`. This matches the reference `fect`
  package ("units whose number of untreated periods <1 are dropped
  automatically") and `didgpu_bacon()`, which already did this.
  On the same DGP all three methods now recover the null.

  Anyone who ran `didgpu_fect()` on a panel containing always-treated
  units should re-run: the ATT was biased upward by their level.

- **Note on `fect` parity.** `didgpu_fect(method = "fe")` agrees with
  `fect::fect()` definitionally but NOT bit-for-bit at the default
  `tol = 1e-5`: the alternating-projections fit stops early, leaving a
  gap of ~6e-5 on a 60-unit panel (~2e-6 at `tol = 1e-8`, ~2e-8 at
  `tol = 1e-12`). Pass a tighter `tol` when exact agreement matters.
  `fect` has been added to `Suggests` so the comparison is now covered
  by the test suite.

- **Reference parity restored on UNBALANCED panels.** Baseline treatment
  `d_sq` was taken from the *global* first period rather than each
  group's own first period with non-missing treatment. A group entering
  the panel late therefore had `d_sq = NA` (the internal balancing merge
  creates the row but leaves treatment missing), which (a) made `F_g`
  fall through to `T_max + 1`, reclassifying the group as a
  never-switcher, and (b) propagated the NA into the `(time, d_sq)`
  cohort-key encoding used by the fast backends, forcing a fallback
  branch that grouped cohorts differently from the reference. Backends
  `"r"`, `"cpu"` and CUDA consequently disagreed with `DIDmultiplegtDYN`
  on any panel with unequal group lengths -- max |diff| 1.8e-01 on a
  60-unit reprex -- while `backend = "reference"` happened to agree.
  Baseline treatment now uses each group's own first non-missing period,
  matching the reference definition exactly; the same reprex now agrees
  to 5.6e-17. Balanced-panel results are bit-identical to before.

  This went undetected by the randomized differential suite because
  `didgpu_simulate_panel()` could only emit balanced panels, on which
  the two definitions coincide. The simulator gained a
  `late_entry_frac` argument (default `0`, RNG-stream preserving) and
  `tests/testthat/test-unbalanced-parity.R` now checks every CPU backend
  against the reference on an unbalanced panel.

  Users who ran didgpu on an unbalanced panel with any backend other
  than `"reference"` should re-run: point estimates were affected.

- **Bootstrap SEs survive partial-NA iterations.** A kept bootstrap
  iteration can carry NA at some horizons (its resample has switchers
  overall but none reaching horizon j). Plain `sd()`/`cov()` then
  returned NA — with few switchers this silently wiped out EVERY SE, CI
  and joint p-value while the point estimates looked fine. Per-horizon
  SEs are now computed from the finite draws with a minimum bootstrap
  support of 30 finite draws per horizon (below that the SE stays NA,
  honestly), and the stored `$coef$vcov` uses a pairwise-complete
  covariance.

- **Joint tests: near-singular covariance now warns.** With many horizons
  and few switchers, the bootstrap covariance behind the omnibus
  joint-effects/placebo chi-square can be near-singular; the statistic is
  then numerically unstable (tiny eigenvalues amplify arbitrary linear
  combinations). `.joint_pvalue()` now excludes horizons with inadequate
  bootstrap support and emits a warning when `rcond(V) < 1e-10`, advising
  a low-dimensional prespecified test (e.g. the leads nearest treatment)
  instead of the omnibus p.

- **CUDA: unbalanced panels with late-entrant groups no longer crash with
  error 700.** Groups unobserved at the global first period have NA
  baseline treatment (`d_sq`); on the CPU path they fall out of every
  cohort mask and contribute zero, but on the CUDA path the NA flowed
  through `match()` into `cohort_key` as `NA_integer_`, which reaches the
  kernels as `INT_MIN` and caused an illegal memory access
  (`CUDA DID kernel failed with code 700`) in
  `k_finalize_dist_and_kernel` — and error 700 poisons the CUDA context,
  so every subsequent didgpu call in the process failed too. Any
  real-world unbalanced panel (units entering the sample over time) hit
  this immediately. NA-key rows are now parked in a padding cohort whose
  kernel contribution is identically zero, matching CPU semantics
  bit-for-bit. Diagnosed with `compute-sanitizer` (invalid 8-byte global
  read at `base + INT_MIN * 8`); regression-tested in
  `test-cuda-late-entrants.R`, including a context-not-poisoned check.

- **Bootstrap aggregation no longer crashes on degenerate resamples.**
  Under sparse-switching treatments (few clean switchers), a bootstrap
  resample can contain no valid switcher cell at some horizon, yielding a
  zero-length effects vector. `.aggregate_to_result()`'s `vapply()` calls
  hard-required full-length vectors, so a single such iteration aborted
  the entire estimation with `values must be length K ... result is
  length 0` — and because the failure probability grows with
  `bootstrap_reps`, exactly the large-rep runs users want for final
  inference were the ones crashing. Degenerate iterations are now dropped
  with a warning that reports the count; SEs and the bootstrap covariance
  use the surviving iterations (standard failed-resample practice).
  Found while re-estimating thin subsamples at 2,000 reps; e.g. a
  78-DAO/4-switcher spec had 103/2,000 degenerate resamples.

## CRAN resubmission fixes

- **`test-parallel.R` now skips on CRAN.** The previous submission's
  pretest exceeded CRAN's 2-core cap (`_R_CHECK_LIMIT_CORES_`) despite
  the test honouring the env var, producing the only test failure.
  Both `test_that` blocks in `test-parallel.R` now call
  `skip_on_cran()`. The parallel path is fully exercised in our GitHub
  Actions CI matrix on every push.
- **README links rewritten to absolute GitHub URLs.** Two references
  (`WINDOWS_BUILD_STATUS.md`, `BENCHMARKS.md`) are intentionally excluded
  from the source tarball via `.Rbuildignore`; the README now points to
  their canonical GitHub URLs so they resolve from the rendered README
  on CRAN. One stale link (`../NOTES_did_gpu_checkpointed.Rmd`, an
  out-of-tree file that no longer exists) was removed.
- **DESCRIPTION typography.** Single-quoted the software-name references
  `'Rcpp'` and `'CUDA'` per CRAN convention.

# didgpu 0.1.1

## Bug fixes

- **CUDA: support Blackwell (`sm_120`) GPUs.** The GPU build previously
  targeted only Turing–Hopper (`sm_75`–`sm_90`) with no PTX fallback, so on
  Blackwell parts (RTX PRO Blackwell, RTX 50xx) the kernels had no device
  image and effect estimates silently came back as **all zeros** while the
  CPU/R backends were correct. Added `compute_120,sm_120` plus a
  `compute_120` PTX target to both `src/Makevars` and `src/Makevars.win`.
  GPU effects are again bit-identical to the CPU/R backends on Blackwell
  (verified on an RTX PRO 5000 Blackwell, CUDA 13.2).
- **Windows build: ship `src/didgpu_cuda.def`.** The DLL export list was
  git-ignored (`/src/*.def`), so clean checkouts — and the r-universe /
  `install_github` source build — could not link `didgpu_cuda.dll`. It is
  now tracked.
- **Windows build: locate CUDA runtime DLLs under `bin/x64`.** CUDA 13
  moved the redistributable DLLs from `bin/` to `bin/x64/`; the bundling
  step now searches both so `didgpu_cuda.dll` loads at runtime.

## Reference parity with `DIDmultiplegtDYN` 2.3.x

didgpu was originally validated bit-for-bit against an older
`DIDmultiplegtDYN`. Two of its outputs were deliberately changed upstream;
didgpu now tracks the current (fixed) behavior:

- **`predict_het` standard errors now use HC2.** The reference switched the
  heterogeneity-regression variance from HC1 to HC2
  (`sandwich::vcovHC(type = "HC2")`) in v2.3.1 ("explicit CI formulas"
  fix). didgpu now does the same (new `sandwich` dependency), so the
  predict_het `SE`/`t`/`LB`/`UB`/`pF` columns match again.
- **Placebo `N` counts each contributing cell once.** For bidirectional
  panels the reported placebo sample size was double-counting controls
  shared between the switcher-in and switcher-out comparisons. It now uses
  the in-direction count, matching the reference's per-row
  `coalesce(in, out)` combiner across both `in>out` and `out>in` panels.
  Point estimates were never affected.

## Correctness fixes (found by randomized differential testing vs the reference)

- **`same_switchers`: the placebo now uses the same restricted switcher set
  as the effects.** Under `same_switchers = TRUE` the placebo block was
  computed on the full switcher set rather than the consistent-switchers
  subset, producing a biased placebo estimate (and inflated placebo `N`)
  relative to `did_multiplegt_dyn`. The placebo distance now honours the
  effects-based `still_switcher` restriction (the reference derives the
  placebo distribution from the same_switchers-gated effect distance), so
  placebo estimates and counts match again.
- **`trends_lin`: no longer crashes on panels with zero estimable effects.**
  On short panels where no group has the `F_g-2` pre-period that
  `trends_lin` requires, result aggregation crashed with
  "length of 'dimnames' [1] not equal to array extent" (an unguarded
  row-name build on a 0-row table). It now returns an empty, no-estimable-
  effects result cleanly.
- **Weighted `N` / `Switchers` columns now separate unweighted counts from
  weighted sums.** On weighted panels (`weight =`) the four reported count
  columns were all populated from the same (weighted, truncated) switcher
  mass, so `N` / `Switchers` reported weighted sums instead of observation
  counts, and `N.w` lost the fractional weight (each cell's weight was
  floored before summing). The estimator now reports `N` and `Switchers` as
  the true unweighted observation / switcher counts and `N.w` /
  `Switchers.w` as the (unfloored) weighted sums, matching
  `did_multiplegt_dyn` exactly. Point estimates were never affected — the
  weighted switcher mass that drives the Neyman pooling and ATE weights is
  unchanged; only the reported count columns moved. On unweighted panels all
  four columns coincide as before (bit-identical output).
- **Weighted estimates: Neyman direction-pooling no longer truncates the
  switcher mass.** On weighted panels with switchers in *both* directions, the
  per-direction weighted switcher mass was floored to an integer before being
  used as the Neyman pooling weight (`w_in = N_in / (N_in + N_out)`) and as the
  across-horizon ATE weight, biasing the pooled event-study estimates by
  ~1e-3. The exact (unfloored) mass is now used throughout the estimate path,
  matching `did_multiplegt_dyn`. Single-direction (`switchers = "in"/"out"`)
  and unweighted panels were never affected (the weight is 0/1 or the mass is
  already integer). Found by randomized weighted×flag differential testing.
- **Reported effect/placebo count: drop trailing unestimable horizons.**
  didgpu's horizon clamp uses each group's own data availability (`max L_g`),
  which can be one step more permissive than `did_multiplegt_dyn`'s
  cohort-level `T_g` clamp. When that extra horizon has no switcher reaching it
  (NA estimate, zero switchers), the reference omits the row; didgpu now trims
  the trailing block of such unestimable effect/placebo rows so the reported
  horizon count matches. Estimable horizons, estimates, and all four count
  columns are unchanged. Verified across 359 boundary-stress panels
  (switchers in/out/both × effects 4–5 × weighted/unweighted): zero horizon
  count mismatches and, importantly, no case where didgpu reported an extra
  horizon with positive switchers. Found by overnight differential testing.
- **Placebos: drop the whole block when every horizon is unestimable.** On
  weighted `trends_lin` (and other short-panel cases) the placebo block can
  be entirely unestimable — no group has the `F_g - q - 1` pre-period any
  placebo horizon needs. `did_multiplegt_dyn` returns `NULL` placebos in that
  case; didgpu used to emit a placebo matrix of NAs. Aggregation now drops
  the placebo block when every placebo point estimate is NA, matching the
  reference. Found by extended weighted×flag differential testing
  (4 of 4 affected panels now match; all other scenarios untouched —
  verified 30 seeds × 14 scenarios = 420 weighted comparisons, 0 fails).

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
- CUDA: all three inner regressions (OR / IPW / DR) run on the GPU
  (`src/cuda_cs_inner.cu`) with per-row influence functions. OR uses
  an in-thread Cholesky per cell; IPW/DR add a per-cell IRLS logistic
  propensity model replicating `stats::glm.fit` (ATT agrees with R to
  ~1e-8), and DR layers on the outcome-regression augmentation. With
  per-cell IFs, the cluster + multiplier bootstrap SEs all run on the
  GPU — the DR cluster bootstrap is ~192x faster than R.
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
- `didgpu_equivalence(fit, delta)` — pre-trends equivalence (TOST) test on
  the placebo estimates. Instead of "failed to reject a zero pre-trend"
  (weak, and worst exactly when underpowered), it tests
  H0: |pre-trend| >= delta and REJECTING is positive evidence the
  pre-trend is within +/- delta. Reports per-horizon and joint
  (intersection-union) verdicts plus the smallest defensible margin
  (`breakdown_delta`). Mirrors `didgpu_fect_equivalence()`.
- `didgpu_joint_placebo(fit, horizons)` — the joint chi-square placebo
  test (`p_jointplacebo`) restricted to a chosen pre-treatment window,
  reusing the stored bootstrap covariance. Test parallel trends only over
  the leads you care about; the full-window call reproduces the headline
  `p_jointplacebo` exactly.
- `didgpu_bacon()` — Goodman-Bacon (2021) decomposition of the static TWFE
  DiD into its 2x2 timing-group comparisons, with the total weight on
  "forbidden" already-treated-control comparisons as the bias diagnostic.
  Validated by the exact identity (weighted 2x2 sum == the TWFE
  coefficient from `didgpu_twfe()`). Balanced, binary, absorbing panels.
- `didgpu_did_static()` — de Chaisemartin & D'Haultfoeuille (2020) DID_M
  instantaneous estimator. Unlike the staggered-adoption methods it allows
  treatment to turn on AND off (non-absorbing): it compares each switcher's
  period-over-period outcome change to same-baseline stayers and averages
  over all switch events, with a cluster bootstrap SE. Native
  reimplementation; cross-checked against `DIDmultiplegt::did_multiplegt`.
- `didgpu_freyaldenhoven()` — Freyaldenhoven, Hansen & Shapiro (2019)
  pre-event panel event study. `estimator = "OLS"` is the two-way FE
  event study; `estimator = "FHS"` adds an auxiliary proxy covariate as an
  endogenous regressor and 2SLS-instruments it with a far policy lead to
  purge a confound that generates pre-trends. Native reimplementation of
  the first-difference parameterization; coefficients match
  `eventstudyr::EventStudy` (OLS and FHS) to machine precision.
- `didgpu_cs_continuous()` — Callaway, Goodman-Bacon & Sant'Anna (2024)
  difference-in-differences with a CONTINUOUS treatment (dose). Estimates
  the dose-response curve: the level effect ATT(d) and the causal response
  ACRT(d) = ATT'(d), via a B-spline regression of the within-unit outcome
  change on the dose, vs a never-treated comparison; multiplier-bootstrap
  SEs. Native reimplementation (spline basis via splines2); ATT(d)/ACRT(d)
  match `contdid::cont_did` exactly.
- `didgpu_did_continuous()` — de Chaisemartin & D'Haultfoeuille (2024)
  continuous treatment with NO STAYERS. When the dose changes for (almost)
  every unit there is no pure control group, so identification is in first
  differences: with dY, dD the within-unit changes, the common trend
  E[dY|dD=0] is recovered from quasi-stayers (dD near 0), giving the level
  effect effect(d) = E[dY|dD=d] - E[dY|dD=0] and the average causal response
  ACR(d). `estimator = "parametric"` fits a polynomial in dD (sqrt(n));
  `estimator = "nonparametric"` is a local-linear (kernel) fit (n^2/5; flagged
  EXPERIMENTAL — no maintained R reference exists to bit-validate it).
  Multiplier-bootstrap SEs. Both estimators validated by simulation against a
  known dose-response.
- All eight auxiliary estimators above were benchmarked to verify they belong
  on the CPU (none has a GPU-amenable hot path; see `BENCHMARKS.md`). The audit
  also caught and fixed two quadratic bootstraps: `didgpu_did_static`'s cluster
  bootstrap was O(n_units^2) per replicate (pre-splitting by cluster makes it
  O(n); 12-23x faster, bit-identical SEs), and `didgpu_did_continuous` no longer
  recomputes the O(n^2) overall-ACR on every bootstrap replicate (nonparametric
  bootstrap ~675x faster; reported effect(d)/ACR(d) unchanged).

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
- `"cuda"` — **live on Windows and Linux/WSL2** (built + verified
  end-to-end on an NVIDIA RTX 4000 Ada, CUDA 12.6; bit-identical results
  on both). On Windows it needs **no admin rights**: a user-local CUDA
  toolkit plus a two-DLL split (`didgpu_cuda.dll` built by nvcc/MSVC,
  the R-facing `didgpu.dll` built by Rtools/MinGW, bridged by a pure-C
  ABI) sidesteps the MinGW↔MSVC link barrier — see
  `WINDOWS_BUILD_STATUS.md`. Live GPU paths: the CS cluster
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

## Validation

- **Monte Carlo CI coverage by switcher count** (`tools/validation/`): with
  the cluster bootstrap, 95% CIs achieve ~94% coverage at 8+ switchers,
  ~87% at 4, but only ~66% at 2 — the classic few-treated-clusters
  failure, not specific to didgpu. Point estimates are unbiased at every
  count. Practical rule: treat estimates identified off fewer than ~5
  clean switchers as diagnostics; their bootstrap CIs materially
  undercover. The degenerate-resample drop does not distort coverage
  where support is adequate (8-13 switchers: nominal).
