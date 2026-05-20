# didgpu

A GPU-capable, checkpointed hub for causal inference on panel and cross-sectional data.
Currently includes **five estimator families plus a sensitivity layer**:

- **DIDmultiplegtDYN-style** dynamic difference-in-differences (de Chaisemartin & D'Haultfoeuille 2024) — `didgpu()`. Bit-for-bit numerical equivalence with the reference [DIDmultiplegtDYN](https://cran.r-project.org/package=DIDmultiplegtDYN) package across every commonly-used option.
- **Callaway & Sant'Anna (2021)** staggered DiD with group-time ATT(g, t) — `didgpu_cs()`. All three inner estimators (OR / IPW / DR-doubly-robust), both control groups (never- / not-yet-treated), covariate adjustment, pre-treatment placebos with joint test, four aggregations (event-study / group / calendar / overall), cluster or multiplier wild bootstrap. Cross-validated against the reference [did](https://cran.r-project.org/package=did) package.
- **fect** counterfactual-prediction estimators (Liu, Wang & Xu 2024) — `didgpu_fect()` with `method = "fe"` (two-way FE), `"ife"` (interactive fixed effects, Bai 2009), and `"mc"` (matrix completion, Athey et al. 2021). Placebo + equivalence tests + CV-based lambda for MC.
- **TestMechs** (Kwon & Roth 2026) sharp mediation testing — `didgpu_test_sharp_null()`. All three test methods (CS / ARP / FSST), binary AND multi-level mediators (K >= 2).
- **HonestDiD** (Rambachan & Roth 2023) sensitivity analysis — `didgpu_honest_did()`. Bounds post-treatment effects under user-specified restrictions on pre-trend violations; reports the "breakdown" parameter at which the conclusion flips. Both `"M"` (smoothness) and `"RM"` (relative magnitudes) methods. Wraps the reference [HonestDiD](https://cran.r-project.org/package=HonestDiD) package.

Plus the **naive TWFE baseline** every applied paper reports for comparison — `didgpu_twfe()`, a two-way fixed-effects dynamic event study (distributed-lag form) with cluster-robust SEs. Its `Effect_k` / `Placebo_j` output mirrors `didgpu()` so the bias-prone TWFE estimate sits right next to the robust one. (Point estimates are bit-exact to `lm()`/`fixest`; intended as a baseline — TWFE is biased under heterogeneous effects.)

Designed for long-running econometric work: per-cell checkpointing to disk, resumable runs after crash or OOM, configurable backends (pure R, optional CUDA), and bit-for-bit numerical equivalence with the reference packages.

## Status

| Feature                                   | Status   |
|-------------------------------------------|----------|
| Binary treatment, on/off                  | ✅ done  |
| Multiple post-treatment event-times       | ✅ done  |
| Pre-treatment placebos                    | ✅ done  |
| Both switcher directions (in + out)       | ✅ done  |
| `switchers = "in"` / `"out"` restriction  | ✅ done  |
| ATE (cumulative effect)                   | ✅ done  |
| Cluster bootstrap (group or custom)       | ✅ done  |
| `weight=` column (weighted DiD)           | ✅ done (bit-identical to reference) |
| Per-cell checkpoint + resume              | ✅ done  |
| Reference-compatible output structure     | ✅ done  |
| Multivalued discrete treatment (3+ levels) | ✅ done (bit-identical, no extra arg needed) |
| Continuous treatment (`continuous=K`)     | ✅ done (bit-identical, k = 1, 2 tested) |
| `controls=` (covariates) — point estimates + bootstrap SEs | ✅ done (bit-identical to reference) |
| `trends_nonparam=` (per-cohort trend column) | ✅ done (bit-identical) |
| `only_never_switchers=` (strict control set) | ✅ done (bit-identical) |
| `same_switchers=` (require valid controls at every event-time) | ✅ done (bit-identical) |
| `dont_drop_larger_lower=` (non-monotone group handling) | ✅ done (bit-identical) |
| `weight=` (weighted DiD)                  | ✅ done (bit-identical) |
| `normalized=TRUE` (per-unit-of-treatment effects) | ✅ done (bit-identical; binary, multivalued, continuous; placebos too) |
| `trends_lin=TRUE` (group-specific linear trends) | ✅ done (bit-identical to reference; effects + placebos + weight + controls + normalized + sample-size columns) |
| `predict_het` (heterogeneity regression) | ✅ done (bit-identical to reference; HC1 robust SEs + joint F-test) |
| `same_switchers_pl` (placebo-side same-switchers gate) | ✅ done (bit-identical to reference) |
| `didgpu_by()` (subgroup-by-subgroup estimation) | ✅ done (wraps didgpu(); per-subgroup checkpoints) |
| `didgpu_by_path()` (treatment-trajectory subgroup analysis) | ✅ done (mirrors reference's by_path argument) |
| `n_workers=` (parallel bootstrap)         | ✅ done (bit-identical to sequential) |
| CUDA backend                              | ✅ live on Linux/WSL (built + verified end-to-end on an RTX 4000 Ada; see GPU acceleration below) |
| Rcpp+Eigen CPU backend                    | 🟡 scaffolded (smoke .cpp compiles; real port TBD) |

For the supported subset (binary, no controls), the r-backend's output matches the reference bit-for-bit on point estimates, SEs, ATE, and the four sample-size columns. See `tests/testthat/test-r-backend.R`, `test-bidirectional.R`, and `test-reference-parity.R` (100+ assertions, all green).

## GPU acceleration

The CUDA backend is built and verified end-to-end on Linux/WSL2 (NVIDIA
RTX 4000 Ada, CUDA 12.6). Set `backend = "cuda"` on a supported call to
use it; every GPU path falls back transparently to the R implementation
when CUDA is unavailable **or when the GPU would be slower** (see the
fect note below), so `backend = "cuda"` is always safe.

### Where the GPU helps — and by how much

| Path | GPU status | Speedup vs R | Notes |
|------|-----------|--------------|-------|
| `didgpu_cs()` **cluster bootstrap** | ✅ live | **179–228×** | Influence-function shortcut; the headline win. R re-runs the full estimator per replicate (~25–50 s for B=200); CUDA does one `(B × n_units) @ (n_units × n_cells)` product (~0.1 s). |
| `didgpu_cs()` **IPW / DR cluster bootstrap** | ✅ live | **173–196×** | DR (the doubly-robust gold standard) at ~192×. Per-cell IRLS logistic propensity kernel matches R's `glm.fit`; ATT agrees to ~1e-8. |
| `didgpu_cs()` multiplier (wild) bootstrap | ✅ live | 1.2–1.7× | R is already IF-based; GPU win is bounded by the IF-matrix copy. |
| `didgpu_cs(est_method = "OR"/"IPW"/"DR")` point estimate | ✅ live | ~1× | Matches R (1e-12 no-cov OR, 1e-6 with covariates / IRLS). CS inner regressions are small, so H2D/D2H roughly cancels the compute win — the bootstrap is where the GPU pays off. |
| TestMechs bootstrap | ✅ live | 4–18× | Nonparametric partial-density bootstrap on GPU; bootstrap moments match R within Monte-Carlo error. |
| `didgpu_fect(method = "mc")` at scale | ✅ live | 3.6–7.9× | Full-SVD matrix completion; engages for large balanced panels (`n_units ≥ 2000`). Verified correct to 2.2e-10. |
| `didgpu_fect()` (fe / ife / small mc) | 🔵 size-gated | ~1× (small panels) | GPU SVD only helps very large panels; below the gate it transparently uses R's LAPACK (far faster for small matrices). ife uses R at all sizes. |

Full numbers and methodology in [`BENCHMARKS.md`](BENCHMARKS.md).

### The headline

For an applied researcher running `didgpu_cs(bootstrap_reps = 1000,
bootstrap_kind = "cluster")` on a typical panel, the CUDA path turns a
~3–4 minute job into well under a second — fast enough that re-running
after a spec tweak is interactive instead of a coffee break.

Every CUDA path is pinned against its R counterpart in
`tests/testthat/test-cuda-equivalence-grid.R` (142 assertions) so the
GPU and CPU results stay in lock-step.

## Install

```r
# From the package source directory
install.packages("didgpu", repos = NULL, type = "source",
                 INSTALL_opts = "--no-multiarch")
```

For optional CUDA support, install the NVIDIA CUDA Toolkit so `nvcc` is on PATH, then reinstall. See [`inst/doc/cuda_setup.md`](inst/doc/cuda_setup.md) for Windows details.

## Quick start

```r
library(didgpu)

# A simulated panel with known event-time profile.
p <- didgpu_simulate_panel(n_units = 100L, n_periods = 20L,
                            frac_treated = 0.6,
                            tau_profile = c(0.5, 1.0, 1.2, 1.0),
                            seed = 17L)

# Estimate ATT at horizons 1..4 with one placebo, 100 bootstrap reps,
# saving each iter to disk so we can resume if anything blows up.
fit <- didgpu(
  df = p,
  outcome   = "Y", group = "unit", time = "period", treatment = "D",
  effects   = 4L, placebo = 1L,
  bootstrap_reps = 100L,
  checkpoint_dir = "checkpoints/example",
  backend   = "r"
)

print(fit)
```

```
didgpu result
  backend         : r
  effects         : 4
  placebos        : 1
  bootstrap reps  : 100 (used 101 cells)
  checkpoint_dir  : checkpoints/example

Effects:
         Estimate     SE  LB.CI  UB.CI    N  Switchers   N.w  Switchers.w
Effect_1   0.5079 0.0671 0.3764 0.6394  364         44  364           44
Effect_2   0.8701 0.0699 0.7332 1.0070  321         40  321           40
Effect_3   1.0478 0.0758 0.8993 1.1963  281         36  281           36
Effect_4   ...

Placebos:
          Estimate     SE  ...
Placebo_1  -0.1778 0.1110 ...

Joint test of all effects:  p = 0.0000
Joint test of all placebos: p = 0.4530
```

## Checkpoint & resume

A long bootstrap on a real panel might take an hour or more. If it crashes (OOM, power loss, kernel panic), the checkpoint directory holds every completed iter:

```
checkpoints/example/
├── meta.json         # the run's config + panel hash
├── manifest.csv      # one row per committed cell
└── cells/
    ├── b0000.rds     # point estimate
    ├── b0001.rds     # bootstrap rep 1
    └── ...
```

Re-invoking with the same call resumes:

```r
fit <- didgpu(p, ..., checkpoint_dir = "checkpoints/example", resume = TRUE)
# [didgpu] resuming checkpoints/example: 47/101 cells already done
# [didgpu] cell b=47    1.71s   (48/101 total)
# ...
```

Re-running with the same seed and config produces identical aggregated output (locked in by `test-reference-parity.R::resume yields identical aggregate`).

You can also extend an existing checkpoint with more bootstrap reps:

```r
# Started with 100 reps; want 200 without re-running the first 100.
didgpu_bootstrap_more("checkpoints/example", df = p, extra_reps = 100L)
```

Or simply re-invoke without re-stating every argument:

```r
# Reads outcome/group/time/treatment/effects/placebo/... from meta.json,
# validates that `p` has the same panel hash, resumes from where it left off.
didgpu_resume("checkpoints/example", df = p)
```

## Working with results

`didgpu()` returns a `didgpu_result` object that plays well with the
standard R model-accessor methods:

```r
fit <- didgpu(p, "Y", "unit", "period", "D",
               effects = 4L, placebo = 2L, bootstrap_reps = 100L)

coef(fit)                       # named vector of estimates
confint(fit)                    # 2-column CI matrix at the stored ci_level
vcov(fit)                       # bootstrap covariance of (effects, placebos)
tidy(fit); glance(fit)          # broom-compatible
plot(fit)                       # base-R event-study plot with error bars
```

## Backends

```r
didgpu_backend_info()
#     backend available                                                                notes
# 1 reference      TRUE                                            DIDmultiplegtDYN 2.2.0
# 2         r      TRUE binary; effects + placebos + controls + weight + trends_nonparam
# 3       cpu     FALSE                                              stub; needs Rcpp port
# 4      cuda     FALSE          stub; install CUDA Toolkit with full headers + nvcc
```

- `"reference"` — delegates to `DIDmultiplegtDYN::did_multiplegt_dyn`. Always available if the reference is installed. Used as the parity oracle.
- `"r"` — standalone pure-R port using `data.table` primitives. 14–60× faster than the reference depending on panel size (see benchmark below). **Bit-identical** to the reference across every commonly-used DIDmultiplegtDYN option (binary / multivalued / continuous treatment; `controls`, `weight`, `trends_nonparam`, `trends_lin`, `normalized`, `predict_het`, switcher restrictions, `only_never_switchers`, `same_switchers`, `same_switchers_pl`, `dont_drop_larger_lower`, all sample-size columns, cluster-bootstrap SEs at the same seed).
- `"cpu"` — Rcpp port of the per-event-time inner kernel for the binary-no-controls case. **Bit-identical** to the `"r"` backend. ~1.5-3× faster on top of `"r"` (≈ 100-190× vs reference). Falls back to `"r"` transparently for unsupported feature combinations (controls, weight, normalized, trends_lin, same_switchers, predict_het, continuous, trends_nonparam).
- `"cuda"` — scaffolded but blocked on a CUDA Toolkit install with full headers. The kernel, host launcher, R glue, and build infrastructure are all in place; it will compile and run as soon as `nvcc` is on `PATH` and the package is reinstalled from source.

`backend = "auto"` picks the best available, preferring cuda > cpu > r > reference.

## Performance

`backend = "r"` vs. reference (`DIDmultiplegtDYN::did_multiplegt_dyn`) per single point-estimate fit, median of 3 runs after warmup:

| Panel (units × periods) | Rows    | cpu-backend | r-backend  | reference   | cpu vs ref |
|---|---|---|---|---|---|
| 100 × 20                | 2 K     | 0.02 s      | 0.07 s     | 0.98 s      | 58×        |
| 200 × 50                | 10 K    | 0.02 s      | 0.10 s     | 2.48 s      | 103×       |
| 500 × 100               | 50 K    | 0.06 s      | 0.28 s     | 13.84 s     | 244×       |
| 1000 × 200              | 200 K   | 0.26 s      | 0.81 s     | 49.07 s     | **188×**   |

For a 100-rep cluster-bootstrap on a 200 K-row panel:
- reference: ~70 minutes
- r-backend, sequential: ~95 seconds
- r-backend, `n_workers = 8`: ~30 seconds

`didgpu(..., n_workers = N)` distributes bootstrap iters across `N` worker processes via `parallel::makeCluster`. Each cell is saved atomically so resume works the same as in sequential mode. Bit-identical to sequential at the same seed.

The speedup widens with panel size because the per-row work in the r-backend is in vectorised `data.table` primitives that are already at C speed; the reference's per-fit overhead (arg validation, multiple internal data copies, optional features) scales worse.

Reproduce with `inst/scripts/bench_scaling.R`.

## Numerical equivalence with the reference

For every commonly-used DIDmultiplegtDYN option, the r-backend agrees with `DIDmultiplegtDYN::did_multiplegt_dyn` on:

- Per-event-time DID estimates (`Effects[, 1]`) — exact match (max abs diff `< 1e-10` across all tested seeds, panel shapes, and option combinations)
- Per-event-time placebo estimates (`Placebos[, 1]`) — exact match
- ATE — match to within a few machine epsilon
- Sample size columns (`N`, `Switchers`, `N.w`, `Switchers.w`) — exact match
- `predict_het` regression block (`Estimate`, `SE`, `t`, `LB`, `UB`, `N`, `pF`) — exact match within 1e-10
- Cluster bootstrap SEs — exact match when both backends are invoked through the same didgpu orchestrator (because both then use the same cluster-resampled panels with the same seed)

This is enforced by a comprehensive test suite (300+ tests across 24 test files) plus an adversarial fuzz harness (`tests/testthat/test-fuzz.R`) that compares both backends on randomized panels across many shapes and option combinations. Run with `Sys.setenv(DIDGPU_FUZZ_N = "100")` for a deep local validation pass.

## Internals

The package design (call graph, kernel formulas, where bugs would hide) is documented in [`inst/doc/reference_internals.md`](inst/doc/reference_internals.md). That document was derived from a careful read of every line of the reference package source; it doubles as the spec for any future backend port.

The high-level architecture is in [`../NOTES_did_gpu_checkpointed.Rmd`](../NOTES_did_gpu_checkpointed.Rmd) ("as-built status" section).

## License

MIT.
