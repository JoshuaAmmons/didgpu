# Overnight session status — 2026-05-19

Autonomous session run while Joshua was away. Summary of what landed,
what was deliberately skipped, and what should happen next.

## What landed tonight (all pushed to `main`)

| Commit  | What |
|---------|------|
| `88fd2c2` | **fect CUDA size gate** — perf fix. Per-iter cuSOLVER SVD was 100–300× *slower* than R for fect's small matrices; gated so `backend="cuda"` falls back to R `svd()` below a size threshold. fect is now ~1× (was 0.01×). |
| `9407e7c` | **Equivalence test grid (#90)** — pins every deterministic GPU path == R path within documented tolerance. |
| `4a97d60` | **README + NEWS GPU docs (#91)** — coverage matrix, benchmark headline, accurate status. |
| `8b414fd` | **testmechs benchmark (#89)** — 4–18× speedup; rounds out the benchmark suite. First STATUS draft. |
| `d828785` | **IPW/DR backend="cuda" safety contract** — +10 assertions: IPW/DR fall back to R bit-identically; cluster bootstrap SE finite/positive for all methods. |
| `264cee9` | **fect at-scale characterization + ife hard-wired to R.** Large-panel benchmark: ife loses even at 8000×50 (0.26×) → ife never uses CUDA; mc wins 3.6–7.9× → keeps the gate. |
| `f48b8c5` | **Bootstrap edge-case tests** — single cluster (exact `colSums(IF)`), n_dims=1, B=1, N(0,1) variance scaling, single cohort. |
| `e8b2ea4` | **mc-at-scale correctness verified** — GPU mc fit agrees with R to **2.2e-10** on a 2500×80 panel, 3.71× faster. |

Earlier in the day (same session, already pushed): Phase 0b (first
end-to-end CUDA), Phase 1 wiring (#79–81), Phase 2 cluster bootstrap
(#82), multiplier bootstrap (#83), CS OR kernel + per-row IF (#84
part 1), and the CS benchmark.

## Benchmark headlines (RTX 4000 Ada, CUDA 12.6, WSL2) — see BENCHMARKS.md

- **CS cluster bootstrap: 179–228×** ← the marquee win. B=200 went from
  ~25–50 s (R) to ~0.1 s (CUDA), via the influence-function shortcut.
- **TestMechs bootstrap: 4–18×** (grows with B; 18× at B=1000).
- **fect_mc at scale: 3.6–7.9×** (≥2000 units; verified correct to 2.2e-10).
- CS multiplier bootstrap: 1.2–1.7× (R already IF-based).
- CS OR point estimate: ~1× (bit-exact; small matrices, copy-bound).
- fect_ife / small fect: ~1× (uses R; GPU has no favourable regime).

## Test status

**611 CI assertions + 17 edge-case assertions, 0 failures, 0
warnings**, 82 skipped (all `DIDmultiplegtDYN` not installed in WSL).
Verified before every commit. Plus two `tools/` verification scripts
(mc-at-scale correctness; large-panel benchmark) that are too slow
for CI.

## What I deliberately did NOT do, and why

1. **#84 part 2 — IPW + DR batched logistic regression on GPU.**
   The no-covariate IPW/DR paths reduce to mean-differences (easy but
   low value — users pick OR for no-cov). The valuable with-covariate
   path needs a per-cell logistic regression (IRLS / Newton-Raphson)
   that must match R's `glm.fit` (with step-halving + deviance
   convergence) closely enough to pass the 1e-6 equivalence bar.
   That's a real numerical-validation task — a near-miss would pass
   some equivalence tests and fail others, leaving murky state. Too
   risky to land unsupervised. **This is the #1 candidate for the next
   supervised session** — the cluster-bootstrap kernel already exists,
   so finishing IPW/DR inner-with-IF would extend the 200× bootstrap
   win to the doubly-robust estimator (the CS gold standard).

2. **#86 / #87 — fused fect kernels.** Ruled OUT by the benchmark
   (documented in BENCHMARKS.md). Even a perfectly device-resident
   loop runs gesvdj on a hundreds-by-tens matrix where the GPU has no
   edge over CPU LAPACK. The size gate is the correct resolution.

3. **#88 — batched LOO refit kernel** and a **didgpu() IF-shortcut
   cluster bootstrap.** Both need a per-unit influence function for
   the de Chaisemartin & D'Haultfoeuille estimator, which isn't
   currently computed anywhere in the codebase. Deriving it is a
   research task, not an overnight one. Until then didgpu()'s bootstrap
   stays on the per-rep R path.

4. **#95–#97 — Windows two-DLL build, r-universe distribution,
   fresh-VM UX validation.** These require Windows build tools (MSVC +
   Rtools) and a Windows VM. Can't be done from WSL. They are the
   critical path to the user's actual goal (Windows + RStudio + GPU,
   `install.packages(...)` and go) and should be the focus once the
   Linux kernel work is at a good stopping point — which it now is.

## Recommended next steps, in priority order

1. **Phase 5 (#95): Windows two-DLL build.** This is what unblocks the
   end-user vision. The C-ABI contract (`inst/include/didgpu_cuda_api.h`)
   is already designed for it. Needs a Windows session with MSVC +
   Rtools + CUDA 12.x. The MinGW↔MSVC ABI bridge is the known hard part
   (documented in ROADMAP.md).
2. **Phase 6 (#96): r-universe + GitHub Releases** so the win-binary is
   installable. Depends on #95.
3. **#84 part 2 (IPW/DR logistic)** — supervised, with the equivalence
   grid as the acceptance test. Extends the 200× cluster-bootstrap win
   to DR.
4. Optionally: derive the de Chaisemartin IF to unlock #88 + didgpu()
   GPU bootstrap (biggest remaining algorithmic win, but research-y).

## How to resume

```
cd "C:\Users\ammonsj\DID GPU\didgpu"
wsl.exe -d Ubuntu -u root -- bash "/mnt/c/Users/ammonsj/DID GPU/didgpu/tools/build-on-wsl.sh"
wsl.exe -d Ubuntu -u root -- Rscript -e "library(didgpu); library(testthat); test_dir('/mnt/c/Users/ammonsj/DID GPU/didgpu/tests/testthat', reporter='minimal')"
```

All WSL invocations use `-u root` to bypass the sudo password prompt.
