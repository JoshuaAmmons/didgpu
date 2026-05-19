# didgpu Timeline Estimates

> Living document. Updated when assumptions change or actuals come in.
> Companion to `ROADMAP.md` (which describes *what* we're building; this
> describes *how long*).
>
> **Last updated:** 2026-05-19 — initial baseline.

## How to read this

Each item has three numbers in **focused work hours** (mine, the AI's
implementation time):

- **O** = optimistic — everything goes right, no surprises
- **M** = most likely — my honest best guess
- **P** = pessimistic — meaningful complications stack up, but not catastrophic

The expected value uses the PERT formula: **E = (O + 4M + P) / 6**. That
weights the most-likely estimate while still accounting for the tails.

A few important notes:

1. **Hours are focused implementation time**, not calendar time. Calendar
   time depends on how many of those hours per week you can put in
   (running my commands, reviewing output, deciding on direction).
   Translation table at the bottom.
2. **Tail risk is real.** P is "things go wrong but we recover" — not
   catastrophic. For tasks where P is much larger than M, the cause is
   usually identified in the *risk factors* section below.
3. **Estimates compound.** Sum-of-P is not the actual worst case (that
   would assume every task hits its tail simultaneously). The realistic
   upper bound for the whole project is closer to the PERT total + 1
   standard deviation, NOT the sum of P's.

---

## Phase 0a — WSL2 dev env setup

| Sub-task | O | M | P | E |
|---|---|---|---|---|
| Run `tools/setup-wsl-env.sh` (your hands-on time) | 0.3 | 0.5 | 2.0 | 0.7 |

**Phase 0a total: O=0.3, M=0.5, P=2.0, E≈0.7 hours.**

Almost entirely your hands-on time, not mine. Variance is network speed
(2-3 GB CUDA Toolkit download) and whether sudo / apt repos behave.

---

## Phase 0b — Linux build with CUDA

| Sub-task | O | M | P | E |
|---|---|---|---|---|
| `R CMD INSTALL .` first attempt | 0.2 | 0.5 | 2 | 0.7 |
| Smoke test existing CUDA paths (saxpy, did, fect_fe) | 0.5 | 1 | 3 | 1.3 |
| Debug any Linux-specific Makevars / path / arch issues | 0 | 2 | 8 | 2.7 |

**Phase 0b total: O=0.7, M=3.5, P=13, E≈5 hours.**

Should be smooth — the .cu files are standard CUDA, Linux is the friendly
platform. P=13h is the "header path mismatch, wrong gcc version, libdevice
not found, debug it" scenario.

---

## Phase 1 — Wire 3 orphan CUDA scaffolds

Per scaffold (×3 — `cs_inner`, `fect_svd`, `testmechs_bootstrap`):

| Sub-task | O | M | P | E |
|---|---|---|---|---|
| Rcpp wrapper for the extern "C" kernel | 1 | 2 | 4 | 2.2 |
| R-side `backend = "cuda"` dispatch in `R/<family>.R` | 0.5 | 1.5 | 4 | 1.8 |
| Parity test (CUDA == R within tolerance) | 1 | 2 | 6 | 2.5 |
| Debug numerical mismatches / kernel bugs | 0 | 2 | 8 | 2.7 |
| **Per-scaffold subtotal** | **2.5** | **7.5** | **22** | **9.2** |

**Phase 1 total (×3): O=7.5, M=22.5, P=66, E≈28 hours.**

Risk: numerical parity is the wildcard. Kernels were scaffolded but never
verified against R. Could be 30-min fixes or could be days of bisecting.

---

## Phase 2 — Fill kernel gaps (the big one)

| Kernel | O | M | P | E |
|---|---|---|---|---|
| Cluster bootstrap (biggest user-facing win) | 8 | 20 | 60 | 23 |
| Multiplier (wild) bootstrap | 4 | 12 | 30 | 14 |
| CS IPW + DR augmentation | 16 | 40 | 100 | 46 |
| CS aggregation + influence accumulation | 8 | 20 | 50 | 23 |
| fect_ife fused FE + SVD alternation | 16 | 40 | 100 | 46 |
| fect_mc CV inner loop | 8 | 20 | 60 | 24 |

**Phase 2 total: O=60, M=152, P=400, E≈176 hours.**

This is the largest phase by far — ~4-5 work weeks at 40h/wk.

Risk factors specific to Phase 2:

- **DR (doubly-robust) CS estimator** is the single highest-variance kernel.
  Combines IPW propensity + outcome regression + augmentation; correctness
  is sensitive to indexing, weighting, missing-cell handling. P=100h
  reflects "rewrite from scratch after realizing the math was off."
- **fect_ife** alternation has convergence subtleties — what if SVD step
  doesn't decrease the objective? What's the tie-breaking on rank
  selection? Could be 2 weeks if we hit a numerical instability.
- The bootstrap kernels are the lowest-risk because the math is simple
  (resample + GEMM); they're mostly engineering.

---

## Phase 3 — Batched LOO + benchmarks

| Sub-task | O | M | P | E |
|---|---|---|---|---|
| Batched LOO refit kernel | 16 | 32 | 80 | 37 |
| Benchmark suite (all 5 families × multiple N, B, K configs) | 8 | 16 | 40 | 19 |

**Phase 3 total: O=24, M=48, P=120, E≈56 hours (~1.5 work weeks).**

Batched LOO is structurally similar to the per-(g,t) inner regression
batched along a new axis — building on existing kernels. Risk: memory
pressure when K is large; might need streaming.

---

## Phase 4 — Verification + docs

| Sub-task | O | M | P | E |
|---|---|---|---|---|
| Equivalence test grid (CUDA == R across params) | 8 | 24 | 80 | 29 |
| README + NEWS + vignette updates for GPU coverage | 2 | 6 | 16 | 7 |

**Phase 4 total: O=10, M=30, P=96, E≈36 hours (~1 work week).**

P=80h on equivalence tests is "finds a real bug in a kernel and we go
back to Phase 1/2 to fix" scenario. Most likely smaller.

---

## Phase 5 — Windows two-DLL build

| Sub-task | O | M | P | E |
|---|---|---|---|---|
| Refactor `didgpu_init.cpp` to host-pointer C-API | 16 | 32 | 80 | 37 |
| Build pipeline: nvcc → didgpu_cuda.dll, dlltool → libdidgpu_cuda.dll.a, Rtools g++ → didgpu.dll | 16 | 48 | 160 | 61 |
| Runtime DLL bundling + DLL search path setup | 8 | 16 | 40 | 19 |

**Phase 5 total: O=40, M=96, P=280, E≈117 hours (~3 work weeks).**

This is the second-largest phase and the second-highest variance. Windows
DLL pipelines have a long tail of weird failure modes:

- `dlltool` quirks producing slightly-wrong import libraries
- Antivirus quarantining the freshly-built DLL during the build
- DLL search-path ordering bugs (system32 vs package libs/x64)
- ABI subtleties (calling convention, struct alignment)
- `cudart64_X.dll` version mismatches vs the user's driver

P=160h on the build pipeline reflects "stuck on a Windows-specific gotcha
for a week despite a clean spec." E=61h is realistic if we're disciplined.

---

## Phase 6 — Distribution

| Sub-task | O | M | P | E |
|---|---|---|---|---|
| Push to GitHub, set up org/repo conventions | 2 | 4 | 8 | 4.7 |
| r-universe entry + builder config | 2 | 4 | 16 | 6.7 |
| Custom Windows binary URL routing | 4 | 8 | 40 | 12 |

**Phase 6 total: O=8, M=16, P=64, E≈23 hours (~3 work days).**

The custom Windows binary URL is the unknown — r-universe is designed
to build everything themselves, and getting them to point at our
prebuilt Release asset for the win.binary slot may require manual config
or back-and-forth with the r-universe team. P=40h covers "they don't
support this; we host our own miniCRAN-style repo instead."

---

## Phase 7 — Fresh-Windows-VM UX validation

| Sub-task | O | M | P | E |
|---|---|---|---|---|
| Spin up VM / fresh R library, install NVIDIA driver | 1 | 2 | 6 | 2.5 |
| Run `install.packages(...)` + smoke test GPU pipeline | 1 | 2 | 8 | 2.8 |
| Fix discovered issues (loops back to Phase 5 if needed) | 0 | 8 | 40 | 12 |

**Phase 7 total: O=2, M=12, P=54, E≈17 hours (~2 work days).**

P=54h is "fresh-VM test reveals a missing DLL or path issue, back to
Phase 5 for a fix, retest cycle." E=17h is realistic.

---

## Grand totals

| | Hours | Notes |
|---|---|---|
| **Optimistic (O)** | 152 | Everything goes right |
| **Most-likely (M)** | 380 | Honest best guess |
| **Expected (PERT)** | 462 | (O + 4M + P) / 6 weighted average |
| **Pessimistic (P)** | 1,095 | Sum of all tails (unrealistic — variance doesn't compound that way) |
| **Realistic upper bound** | ~600 | PERT + ~1σ for the high-variance phases |

### Calendar-time translation

| Work pace | Optimistic (152h) | Expected (462h) | Realistic upper (600h) |
|---|---|---|---|
| Full-time (40h/wk) | 4 wks | 12 wks (~3 mo) | 15 wks (~4 mo) |
| Half-time (20h/wk) | 8 wks | 23 wks (~5 mo) | 30 wks (~7 mo) |
| Evenings + weekends (10h/wk) | 15 wks | 46 wks (~11 mo) | 60 wks (~14 mo) |
| Casual (5h/wk) | 30 wks | 92 wks (~21 mo) | 120 wks (~28 mo) |

So the realistic answer is: **3-7 months if you can dedicate substantial
time, ~1 year if it's evenings and weekends, ~2+ years if it's casual.**

My earlier "5-7 weeks" claim was wildly off. Apologies.

---

## What would shift these estimates

### Things that would shorten

- **Skip Phase 2's harder kernels.** If we're willing to ship with just
  cluster + multiplier bootstrap (the biggest user wins) and defer CS DR,
  fect_ife alternation, fect_mc CV → shaves ~80h off Phase 2 (M=152 → M=70).
- **Skip the batched LOO** (Phase 3). Sequential LOO still works on GPU
  via the existing per-replicate path; the batched version is a perf
  optimization. Saves ~37h.
- **Defer macOS support entirely.** We're already targeting Linux + Win;
  Mac is just "doesn't break the R backend." Already in the plan, no
  savings.
- **Reuse existing R-package CUDA build patterns** (look at how `torch`,
  `gpuR`, `Rcpp.GPU` do Phase 5). Could save 20-40h on Phase 5 if a
  template exists.
- **AI-paired implementation.** I implement, you review/run. Assumed
  throughout; if you're doing more yourself, hours stretch.

### Things that would lengthen

- **CUDA version churn.** If CUDA 14 ships during the project and we want
  to support it, add ~40h of compatibility work.
- **Need to support CUDA 11.x** for users on older hardware. Adds ~30h
  of conditional compilation.
- **CRAN submission.** If you decide later you want this on CRAN proper
  (not just r-universe), CRAN's strict policies around non-CRAN
  dependencies, binary size, and CUDA `SystemRequirements` add ~60-100h
  of policy work.
- **Antivirus quarantine on the built DLL.** Windows Defender and
  enterprise AVs sometimes flag freshly-compiled DLLs. Workaround is
  code-signing the binaries (~$300/yr cert + 5-10h integration).
- **R 4.5 transition.** R changes APIs sometimes; if we hit a transition,
  add ~20h.

### Things that won't move these estimates

- Phase 0a/0b velocity — they're small line items
- Documentation polish — bounded at Phase 4
- Pure-R backend changes — already done; out of scope here

---

## Risk register

| Risk | Likelihood | Impact (hours) | Mitigation |
|---|---|---|---|
| CS DR estimator math has subtle bug requiring rewrite | medium | +40 | Cross-check against `did`/`DRDID` early in Phase 2 |
| fect_ife alternation has numerical instability | medium | +30 | Compare iteration-by-iteration vs the `fect` reference |
| r-universe won't host custom Windows binary URL | medium | +20 | Host own static mini-repo on GitHub Pages as backup |
| Windows DLL bundling has antivirus issues | medium | +30 | Test on fresh VM early, document workaround |
| CUDA 13.2.1 Windows install is still broken on user box at Phase 5 | high | +10 | Use the shadow-tree workaround from earlier session |
| NVIDIA driver update breaks WSL passthrough mid-project | low | +10 | Pin driver version in development env |
| User-side GPU breaks (hardware failure) | low | +∞ | Backup builds work without GPU; package still ships |

---

## Update protocol

When updating this document:

1. Bump the **Last updated** date at the top.
2. If a sub-task completes, replace its O/M/P/E with the **actual hours
   spent** in **bold**, and add a brief note about what was harder/easier
   than expected.
3. If a sub-task estimate changes mid-project (because we learned
   something), put the **old E → new E** in italics so we can audit
   our calibration later.
4. Add new sub-tasks discovered during work as they appear.
5. Sum-of-actuals + remaining estimates gives a live "where are we" view.

Example after Phase 0b completes:

> **Phase 0b — Linux build with CUDA**
>
> | Sub-task | O | M | P | E | Actual |
> |---|---|---|---|---|---|
> | `R CMD INSTALL .` first attempt | 0.2 | 0.5 | 2 | 0.7 | **0.4** |
> | Smoke test existing CUDA paths | 0.5 | 1 | 3 | 1.3 | **0.8** |
> | Debug Linux-specific issues | 0 | 2 | 8 | 2.7 | **0.5** (compute_75 not supported on this driver, dropped from gencode list — 30 min fix) |
> | **Total** | **0.7** | **3.5** | **13** | **5** | **1.7** |
>
> *Faster than expected because Linux is the friendly platform.*
